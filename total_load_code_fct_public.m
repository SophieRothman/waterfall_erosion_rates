function [En_suspt_mm_yr En_suspt_st_mm_yr En_Sklar_mm_yr] = total_load_code_fct_public(D,S,qb,H,tau_star,ustar,tau_star_c,tstage)
% MATLAB Code to calculate total load erosion following methods of Lamb et
% al (2008). Code written by Michael Lamb, modified and commented by Joel 
% Scheingross, December 2013.

% Note: Use this code at your own risk. This code was written to produce
% the figures presented in Lamb et al (2008), and uses numerical
% approximations in place of analytical solution in select cases. Different 
% applications may require modification of the code, such as increasing the 
% resolution used to calculate turbulent fluctuations or in discritizing 
% the law of the wall calcluations. You are responsible for inspecting the
% code and modifying it as neccessary to suit your own applications.

% All references to equations refer to those from Lamb et al (2008)

% Please report errors in the code to jscheingross at gmail dot com.

%Inputs
% D: grain size (meters)
% S: slope (input as a fraction i.e., for a 10% gradient, use 0.1)
% qb: volumetric sediment flux per unit width, note this is the total flux (bed + suspended load) (m2/s)
% H: flow depth (m)
% tau_star: Shields Stress (dimensionless)
% ustar: shear velocity (m/s)
% tau_star_c: critical Shields Stress (dimensionless)

%Outputs
% En_suspt_mm_yr: Instanteous erosion rate (mm/yr) from Lamb et al (2008) ignorning the effects of viscous dampening
% En_suspt_st_mm_yr: Instanteous erosion rate (mm/yr) from Lamb et al (2008) including the effects of viscous dampening
% En_Sklar_mm_yr: Instanteous erosion rate (mm/yr) from Sklar and Dietrich (2004)
%% Constants
den_sed = 2650; %Density of sediment (kg/m3)
den_w = 1000; %Density of fluid (kg/m3)
R=(den_sed-den_w)/den_w; %submerged specific density of sediment
g = 9.81; %gravitational acceleration (m/s2)
nu = 10^-6; %kinematic viscosity of water
von_Karmen = 0.41; %von Karmen's constant

% Bedrock properties
kv = 10^6; %non-dimensional coefficient from eq (3), using value calibrated by Sklar and Dietrich (2004)
young = 5*10^4*10^6; %Generic Young's Modulus for rock as reported by Sklar and Dietrich (2004)
strength = 7 * 10^6; %Generic tensile strength for rock 

%% Set drag coeficient (following Dietrich 1982)
CSF = 0.8;  %1 is for spheres, 0.8 is for natural
PS = 3.5;  %6 is for spheres, 3.5 is for natural

Dstar = (R.*g.*D.^3)./(nu.^2);
X = log10(Dstar);
R1 = -3.76715+1.92944.*X - 0.09815.*(X.^2) - 0.00575.*(X.^3) + 0.00056.*(X.^4);
R2 = log10(1-((1-CSF)./0.85))-(((1-CSF).^2.3).*tanh(X-4.6)) + 0.3.*(0.5-CSF).*((1-CSF).^2).*(X-4.6);
R3 = (0.65-((CSF./2.83).*tanh(X-4.6))).^(1+((3.5-PS)./2.5));
Wstar = R3.*10.^(R2+R1);
ws = (R.*g.*nu.*Wstar).^(1./3);  %Dietrich settling velocity
cdrag = (4/3).*(R.*g.*D)./(ws.^2);

%% Fluid Stresses
%%% The user must set two of the three of slope, height, and transport
%%% stage, and then calculate the third, do this in a seprate code.
% tau_star = tau./((den_sed-den_w).*g.*D); %Shields Stress  (dimensionless)
% tstage = tau_star./tau_star_c; %Transport Stage
% ustar = sqrt(tau./den_w); %shear velocity (m/s)

%% Compute flow velocity with Law of the wall
z0 = 3.*D./30; %set z0, NOTE: Here we're using the input D to set the local roughness instead 
  % of explicitly specifying roughness value. Users may need to modify this
  % for their own specific conditions.
dz = (H-z0)./1000; %space dz evently between z0 and the top of the water column
z = z0:dz:H; %find velocities from z0 to the top of the water column, (for z < z0, log explodes)
flow_vel = (ustar./von_Karmen) .* log(z./z0); %law of the wall (equation 21)
Uf = sum((ustar./von_Karmen) .* log(z./z0)).*dz./H; %depth averaged velocity (equation 22)

%% Compute bed load height and velocity
hb = D.*1.44.*(tstage-1).^0.5; %height of the bed load layer (equation 25)
Us = (R.*g.*D).^0.5.*1.56.*(tstage-1).^0.56; %bed load velocity (equation 24)
Us(Us>Uf)=Uf; % Don't let particle velocity exceed fluid velocity

%% Compute Suspended Load
if hb < H %for cases when the water depth is greater than the bedload height
    b = hb; %bottom of the suspended load layer is the same as height of bedload
    beta = 1; %coefficient in rouse number (equation 27)
    P = ws./(von_Karmen.*ustar.*beta); %Rouse Parameter
    
    %Use a log scale to clean up the integration
    di5 = (log(H)-log(b))./1000; %we'll calculate at 1000 intervals between the bottom of the suspended load and top of the column
    i5 = log(b):di5:log(H); %space the intervals on a log scale
    z = exp(i5); % Take the log intervals and convert it back to linear
    z(length(z))=H; %Make sure the last element is the top of the water column
    dz = diff(z);
    dz = [dz(1) dz];
    chi=sum((((1-(z(z>z0)./H))./(1-(b./H))).*(b./z(z>z0))).^P .* log(z(z>z0)./z0) .*dz(z>z0))...
    ./ (Uf.*H) .* (ustar./0.41);
    
    cb = qb./(Us.*hb +Uf.*H.*chi); % near bed sediment concentration (equation 18)
    
    %Find the concentration profile
    c=0;
    c(1) = cb;
    c(2:length(z)+1) = cb.*(((1-(z./H))./(1-(b./H))).*(b./z)).^P ; %apply Rouse profile (equation 26)
    z=[0 z];
    c(z==H)=0; %make concentration go to zero at top of the water column
    
    % Calculate the fall distance
    gradc(1:length(c))=0;
    gradc(2:length(c)) = -diff(c);
    Hfall = (1./cb).*sum(z.*gradc); % weighted fall distance (equation 32)
    
else %for cases when the bedload saltates throughout the entire water column
    hb = H; %the bedload layer height is the same as the water depth
    cb = qb./(Us.*hb); %we only have bedload so equation (13) applies
    Hfall = hb;
end

if cb == 0 % if there is no near bed sediment concentration, the fall height is 0
    Hfall = 0;
end

%% Settling Velocity (with turbulence)

% Probability for fluctuations
sig = ustar; %assume turbulent fluctuations scale with u* (end of paragraph 35)
dx = sig./100; %the size of bins to subdivide
X = [-6*sig:dx:6*sig]; %spread distribution for six sigma (this is wi, JS thinks)
f = normpdf(X,0,sig); %normal Guasian distribution, centered at zero
X = X./ws;  % normalize by ws to match with psifall below

% Calculate impact velocity due to gravity, accounting for fact that particles may not reach terminal settling velocity
wst = sqrt((4.*R.*g.*D)./(3.*cdrag)); %terminal settling velocity (equation 31)
Scos = cos(atan(S));  %cosine of the bed angle, account for fact that impacts are not normal to the bed
wfall = wst.*Scos.*sqrt(1-exp((-3.*cdrag.*den_w.*Hfall)./(2.*den_sed.*D.*Scos))); %bed normal component of settling velocity (equation 30)

% Add turbulence into settling veloicty
psifall = wfall./ws;
settlematrix = psifall + X;
settlematrix(settlematrix<0) = 0; %no negative impacts
psi_fall3 = sum((settlematrix.^3).*f).*dx;
E1 = psi_fall3 .* cb;  %we'll use this term to calculate erosion rate including turbulent fluctuations

%% Include Viscous Dampening Effect with Stokes number correction
wi_st = settlematrix;
Stokes_thresh = 30; %set the threshold for Stokes dampening, should be between ~10 - 90
wi_st((D.*wi_st.*ws.*den_sed./(9.*nu.*den_w))<Stokes_thresh) = 0; %Set impact velocity to zero for Stokes # < stokes_threshold
psi_fall3_st = sum((wi_st.^3).*f).*dx;
E1_st = psi_fall3_st .* cb;  %we'll use this term to calculate erosion rate including turbulent fluctuations and viscous dampening

%% Calculate bedload transport capacity
qbt = (Us.*hb.*cb); % (equation 13), this is the bed load flux, use this for the cover term in the model
qt = 5.7.*(R.*g.*D.^3).^0.5.*(tau_star- tau_star_c).^(3/2); %bed load transport capacity following Fernandez Luque and van Beek (1976)


%% Calculate non-dimensional erosion rates
% With turbulence, no viscous dampening
A1 = 0.36; %follow Lamb et al, 2008; A1 < 1 accounts for the fact that some of the particles near the bed are advected upward because of lift forces
En_suspt = A1./kv .* (ws./(g.*D).^0.5).^3 .* E1.*(1-(qbt./qt)); %non dimensional erosion rate
En_suspt(En_suspt<0)=0;

% With turbulence and stokes correction (include viscous dampening)
En_suspt_st = A1./kv .* (ws./(g.*D).^0.5).^3 .* E1_st.*(1-(qbt./qt)); %non dimensional erosion rate
En_suspt_st(En_suspt_st<0)=0;

% Sklar and Dietrich (2004) Erosion
En_Sklar = 0.456.*(R.^1.5) .*tau_star_c.^1.5 ./ kv .* (qb./qt) .* (1 - qb./qt) .* ...
    (tstage - 1).*(1-(ustar./ws).^2).^(3./2);
En_Sklar(En_Sklar<0)=0;


%% Extra guards for special cases
if tstage <=1 %if tau_star < tau_c
    En_suspt_st = 0;
    En_suspt = 0;
    En_Sklar = 0;
end

if H < D %if grains are emerged
    En_suspt = 0;
    En_suspt_st = 0;
    En_Sklar = 0;
end

% This next one isn't in the 2008 paper, but (I think) was used in Mike's
% original code 
% Verify that cb is below entrainment capacity
% Rp = (R.*g.*D).^0.5 .* D ./ nu;
% Xi = (ustar./ws) .*Rp.^0.6 .*S.^0.08;
% B= 7.8*10^-7;
% entrainGP = B.* Xi.^5 ./ (1 + (B./0.3 .* Xi.^5)); %Garcia and Parker fmla - with Wright correction
% 
% if cb > entrainGP
%     En_susp = 0;
%     En_suspt = 0;
%     En_suspt_st = 0;
% end

%% Dimensional (mm/yr) erosion rate
Efactor = 1000.*60*60*24*365.*den_sed.*young.*(g.*D).^(3/2)./(strength.^2); % multiply by En to get mm/yr (this is instanteous erosion, include an intermittency factor if you want to compare to long term averages)
En_suspt_mm_yr = En_suspt.*Efactor; 
En_suspt_st_mm_yr = En_suspt_st.*Efactor; 
En_Sklar_mm_yr = En_Sklar.*Efactor; 

end