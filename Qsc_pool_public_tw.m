function [Qsc_pool, h_TW, c_wall_z, Cb, delta, c_bed_vs_r] = Qsc_pool_public_tw(r_pool, ...
    Q, h_pool, H, W, D, h_brink, dr_spacing, dz_spacing, Htl)

%% Notes
% MATLAB Code to calculate plunge pool sediment transport capacity
% following the method of Scheingross and Lamb (2016)

% Code written by J. Scheingross

% Note: Use this code at your own risk. This code was written to produce
% the figures presented in the accompanying paper. Different 
% applications may require modification of the code. You are responsible 
% for inspecting the code and modifying it as neccessary to suit your own 
% applications.

% All references to equations refer to those from Scheingross and Lamb (2016)

% Please report errors in the code to jscheingross .a.t. gmail dot com.

%% Inputs and outputs
%Inputs
% r_pool: plunge pool radius (meters)
% Q: water discharge (m3/s)
% h_pool: plunge pool depth (z_lip - z_sed)(m)
% H: waterfall drop height (m)
% W: Average channel width (m)
% D: grain size (m)
% S: Slope upstream of the waterfall (enter as dimensionless quantity [i.e., 5% slope should be entered as 0.05])
% Cf_river: Friction coefficient in the river, we use Cf_river = 10^-2, but this varies in natural rivers, see Parker 2008 book chapter for a review

%Outputs
% Qsc_pool: plunge pool sediment transport capacity in m3/s

%% Constants
k1 = 1; %constant in eddy viscosity formula (Eq. 28)
k2 = 0.02; %constant in sediment entrainment function (Eq. 30) 
tau_star_c = 0.045; %critical Shields stress 
Cf_pool = 0.001; % plunge pool friction factor
Cd = 2.6; %waterfall jet diffusion coefficient, we use the value used in Flores-Cervantes et al (2006)
gamma = 0; %virtual origin for jet-half width (Eq. 13)
R = 1.65; %submerged sediment density
g = 9.81; %gravitational acceleration (m/s2)
nu = 1e-6; %kinmatic viscosity of water
rho = 1000; %water density (kg/m3)

%% Calculate flow depth at waterfall brink
q = Q./W; %water discharge per unit width (m2/s)
% Fr_n = sqrt(S./Cf_river); % Fr upstream of the waterfall (Eq. 10)
% if Fr_n > 1;
%     Fr_lip = ((0.4+Fr_n.^2)./Fr_n.^2).^(3/2).*Fr_n; %Fr at the waterfall lip for supercritical flow
% else
%     Fr_lip = (1.4./Fr_n.^(2/3)).^(3/2).*Fr_n; %Fr at the waterfall lip for subcritical flow
% end
% h_brink = (q./(sqrt(g).*Fr_lip)).^(2/3); %flow depth at the waterfall brink

%% Calculate tailwater depth (z_water - z_lip)
h_TW = Htl;%(Q./(W.*sqrt(g))).^(2/3); %tailwater depth (Eq. 11)
h_total = h_pool + h_TW; % this is the total water depth (z_water - z_sed)

%% Plunge pool basal shear stress following Stein et al 1993
V_brink = q./h_brink; % Solve for velocity at brink assuming conservation of mass
beta = atan(sqrt(2.*g.*H)./V_brink); %angle of jet impact with water surface (Eq. 3)
ujet = sqrt(V_brink.^2 + 2.*g.*H);% jet velocity at impact (Eq. 5) following Stein and Julien 1993 conservation of energy equation

jet_area = Q./ujet; %jet area at impact with plunge pool water surface following conservation of mass
r_jet = sqrt(jet_area./pi); % waterfall jet diameter at impact with water surface, Eq.(4)

lambda = 2.*Cd.^2.*r_jet.*sin(beta); %length scale of ZOFE, Eq. 2

if (lambda/h_total) > 1 ; %if the pool-floor is the ZOFE
    u_impact = ujet; % Jet impact velocity on the pool floor (Eq. 6a)
else %if the pool-floor is in the ZOEF
    u_impact = ujet.*sqrt(lambda./h_total);% Jet impact velocity on the pool floor (Eq. 6b)
end

tau = Cf_pool.*rho.*(u_impact.^2); %plunge pool basal shear stress (Eq. 1)
tau_star = tau./(1650*g.*D); % plunge pool basal Shields stress
ustar = sqrt(tau./1000); % plunge pool basal shear velocity (m/s)

%% Particle settling velocity (Ferguson and Church 2004)
c1 = 20; %constant from Ferguson and Church 2004
c2 = 1.1; %constant from Ferguson and Church 2004
ws= (R.*g.*D.^2)./((c1.*nu)+sqrt(0.75.*c2.*R.*g.*D.^3)); % Eq. 16

%% Wnet
w_up = Q./((pi.*r_pool.^2) - (pi.*r_jet.^2)); %upward advective velocity (m/s) (Eq. 14)
w_net = ws-w_up;
if w_up>0 %if w_up is greater than 0, check to make sure w_net is not negative  
    w_net(w_net<0)=0; % if w_net < 0 that means that the upward current is greater than gravitational 
                      % settling, so that particles do not gravitationally settle out, set w_net to 
                      % zero for these cases so that the concentration vertically is always equal to cb 
                      %(i.e. there is no vertical decay of sediment concentration above the bed) 
elseif w_up<0 %this means the jet has a larger diameter than the plunge pool
    w_net = 0; %set w_net to zero for these cases so that the concentration vertically is always equal to cb
end

%% Entrainment relation (van Rijn 1984 style approach)
Cb = k2.*((tau_star./tau_star_c)-1).^1.5; % Eq. 30

if tau_star/tau_star_c < 1; %find cases where tau* < tau*c
    Cb = 0; %if shear stress drops below critical, set entrainment to zero
end

if Cb > 0.2; %find excessively large (and physically unreasonable) values of near-bed sed. conc.
    Cb = 0.2; %Don't let near-bed sediment concentration exceed 0.2.
end

%% Bedload layer height
hb = 1.44.*D.*sqrt((tau_star./tau_star_c)-1); %Use the Sklar and Dietrich empirical relation to calculate bedload layer height

if hb < D; %if the bedload-layer height is less than the grain diameter
    hb = D; %the bedload layer height must be at least as high as the grain size
end

%% Kinematic eddy viscosity and Ld
eddy_visc = k1.*ustar.*lambda; %kinematic eddy viscosity (Eq. 28)
L_d = eddy_visc./w_net; %turbulence diffusion length scale (Eq. 20)

%% Calculate half-width and extent of jet-descending region
if h_total > lambda; %if the pool floor is in ZOFE
    % In ZOEF set jet-descending region to scale with depth:
    delta = 2.*0.1.*(h_total-gamma); %Eq. 13a
else
    % In ZOFE set jet-descending region to a constant:
    delta = 2.*0.1.*(lambda-gamma); %Eq. 13b
end

%% Calculate concentration from lip to top of tailwater
h_tw_vec = linspace(h_pool,h_total,100); %set up a vector for water heights 
    % above z_lip and below the watersurface, we'll integrate over this to
    % calculate the averaged sediment concentration at the pool lip

% Calculate the vertical decay of sediment concentration
vertical_decay = exp(-1.*((h_tw_vec-hb)./L_d)); %calculate Eq. 24 for all depths from z_lip to z_water
i2 = h_tw_vec < hb; %check if points in the tailwater are within the bedload layer 
        % this will only happen if the pool is very shallow.  
vertical_decay(i2) = 1; %within the bed load layer, concentration does not vary with z

% Calcluate the radial decay of sediment concentration 
k1_radial = 1./(besseli(0,delta./L_d)+(besseli(1,r_pool./L_d)...
    .*besselk(0,delta./L_d))./besselk(1,r_pool./L_d));
k2_radial = k1_radial.*besseli(1,r_pool./L_d)./besselk(1,r_pool./L_d);
radial_decay = k1_radial.*besseli(0,r_pool./L_d)+k2_radial.*besselk(0,r_pool./L_d); 

    
c_rpool_z = Cb.*vertical_decay.*radial_decay; % (Eq. 27)
% c_rpool_z is a vector of sediment concentration between the downstream
% plunge pool lip and water surface

if r_pool < delta %if we're within the jet-descending region
    c_rpool_z(:) = Cb.*vertical_decay; %make sure there's no radial decay 
end
if (1/L_d) == 0; % these correspond to cases when w_net = 0
    c_rpool_z(:) = Cb; %keep concentration constant when w_net = 0 (that is, particles have no gravitaitonal settling)
end
    
% Take the depth-averaged concentration to computer Qsc_pool 
dz = mean(diff(h_tw_vec));
Qsc_pool = Q.*sum(c_rpool_z).*dz./h_TW; %numerical integration of Eq. 29

%% Check against reference transport value
qsc_star = Qsc_pool./((2.*r_pool).*sqrt(1.65.*g.*D.^3)); %non-dimensional sediment flux
qsc_star_thresh = 2e-5; % Threshold for zero transport based on Parker et al, 1982 assuming Wr = 0.002 and tau*c = 0.045, see Parker and Klingman 1982, eq. 3 
if qsc_star < qsc_star_thresh; %find fluxes lower than the threshold value
    Qsc_pool = 0; %set Qsc_pool = 0 for qsc_star values below the threshold
end

%% Calculate concentration along wall
% dz = h_pool./100;
% h_pos = linspace(dz,h_pool,100);

dz = h_pool./dz_spacing;
h_pos = linspace(dz,h_pool,dz_spacing);


% Calculate the vertical decay of sediment concentration
vertical_decay_wall = exp(-1.*((h_pos-hb)./L_d)); %calculate Eq. 24 for all depths from z_lip to z_water
i3 = h_pos < hb; %check if points in the tailwater are within the bedload layer 
        % this will only happen if the pool is very shallow.  
vertical_decay_wall(i3) = 1; %within the bed load layer, concentration does not vary with z

c_wall_z = zeros(size(h_pos));
if r_pool < delta %if we're within the jet-descending region
    c_wall_z(:) = Cb.*vertical_decay_wall; %make sure there's no radial decay 
elseif (1/L_d) == 0; % these correspond to cases when w_net = 0
    c_wall_z(:) = Cb; %keep concentration constant when w_net = 0 (that is, particles have no gravitaitonal settling)
else
    c_wall_z = Cb.*vertical_decay_wall.*radial_decay; % 
    % c_wall_z is a vector of sediment concentration between the downstream
    % plunge pool lip and the bottom of the pool
end
if length(c_wall_z)==1
    fdshj
end

%% Calculate concentration along the floor
dr = r_pool./dr_spacing;
r_pos = linspace(dr,r_pool,dr_spacing);

top_fraction = besseli(0,r_pos./L_d)+(besseli(1,r_pool./L_d)./...
    besselk(1,r_pool./L_d)).*besselk(0,r_pos./L_d);
bottom_fraction = besseli(0,delta./L_d)+(besseli(1,r_pool./L_d).*...
    besselk(0,delta./L_d))./besselk(1,r_pool./L_d);
c_bed_vs_r = Cb.*top_fraction./bottom_fraction; %

end