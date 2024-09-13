function [c_bed_vs_r, c_wall_vs_z, r_pos, h_pos, w_net, ...
    b, u_impact, ws, Cdrag, lambda, r_jet, w_up, ...
    c_capacity_bed_vs_r, ujet, h_total] = ...
    c_from_Qs_fct_tw(Qs, r_pool, Q, h_pool, H, W, D, h_brink, dr_spacing, dz_spacing, Htl)
dbstop if error
%% Notes
% MATLAB code to calculate concentration of sediment around plunge pool
% boundaries following the method of Scheingross and Lamb, in review and
% Scheingross and Lamb (2016).  This is a modified version of the code
% included in the supplement of Scheingross and Lamb (2016, JGR-ES), all
% equation numbers refer to Scheingross and Lamb (2016) unless otherwise
% noted.

% Code written by J. Scheingross

% Note: Use this code at your own risk. This code was written to produce
% the figures presented in Scheingross and Lamb (in review). Different
% applications may require modification of the code. You are responsible
% for inspecting the code and modifying it as neccessary to suit your own
% applications.


% Please report errors in the code to jscheingross . a . t . gmail dot com.

%% Inputs and outputs
%Inputs
% Qs: the imposed sediment supply
% r_pool: plunge pool radius (meters)
% Q: water discharge (m3/s)
% h_pool: plunge pool depth (m)
% H: waterfall drop height (m)
% W: Average channel width upstream of the waterfall (m)
% D: grain size (m)
% dr_spacing, dz_spacing: discretize the plunge pool walls and floors

%Outputs
% c_bed_vs_r: a vector of concentration along the plunge pool floor from the center out
% c_wall_vs_z: a vector of sed. concentration along the pool walls
% r_pos: radial vector across the floor from center to wall
% h_pos: vertical vector from pool floor to pool lip
% w_net: net particle settling velocity (m/s)
% b: waterfall jet half width (m)
% u_impact: waterfall jet velocity at impact with plunge-pool floor (m/s)
% ws: particle settling velocity (m/s)
% Cdrag: drag coefficient for particles
% lambda: length of zone of flow establishment (m)
% r_jet: waterfall jet radius at impact with plunge-pool water surface (m)
% w_up: upwards return flow velocity in plunge pool (m/s)
% c_capacity_bed_vs_r: concentration of sediment along pool floor at transport capacity
% ujet: waterfall jet velocity at impact with plunge-pool water surface (m/s)
% h_total: the distance from the plunge-pool water surface to the pool floor (m)

%% Constants
k1 = 1; %constant in eddy viscosity formula (Eq. 23)
k2 = 0.02; %constant in sediment entrainment function (Eq. 25)
tau_star_c = 0.045; %critical Shields stress
Cf_pool = 0.001; % plunge pool friction factor
Cd = 2.6; %waterfall jet diffusion coefficient, we use the value used in Flores-Cervantes et al (2006)
R = 1.65; %submerged sediment density
g = 9.81; %gravitational acceleration (m/s2)
nu = 1e-6; %kinmatic viscosity of water
rho = 1000; %water density (kg/m3)

%% Calculate flow depth at waterfall brink
q = Q./W; %water discharge per unit width (m2/s)

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

if (lambda/h_total) > 1 ; %find depths within ZOFE
    u_impact = ujet; % Jet impact velocity on the pool floor (Eq. 6a)
else
    u_impact = ujet.*sqrt(lambda/h_total);% Jet impact velocity on the pool floor (Eq. 6b)
end

tau = Cf_pool.*rho.*(u_impact.^2); %plunge pool basal shear stress (Eq. 1)
tau_star = tau./(1650*g.*D); % plunge pool basal Shields stress
ustar = sqrt(tau./1000); % plunge pool basal shear velocity (m/s)

%% Particle settling velocity (Ferguson and Church 2004)
c1 = 20; %constant from Ferguson and Church 2004
c2 = 1.1; %constant from Ferguson and Church 2004
ws= (R.*g.*D.^2)./((c1.*nu)+sqrt(0.75.*c2.*R.*g.*D.^3));
Cdrag = ((2.*c1.*nu./sqrt(3.*R.*g.*D.^3))+sqrt(c2)).^2; %Ferguson eq. 5

%% Wnet
w_up = Q./((pi.*r_pool.^2) - (pi.*r_jet.^2)); %upward advective velocity (m/s) (Eq. 13)
w_net = ws-w_up;
if w_up>0 %
    w_net(w_net<0)=0; % if w_net > 0 that means that the upward current is greater than gravitational
    % settling, so that particles do not gravitationally settle out, set w_net to
    % zero for these cases so that the concentration vertically is always equal to cb
    %(i.e. there is no vertical decay of sediment concentration above the bed)
elseif w_up<0 %this means the jet has a larger diameter than the plunge pool
    w_up = 0;
    w_net = 0; %set w_net to zero for these cases so that the concentration vertically is always equal to cb
end

%% Bedload layer height
z_mixed = 1.44.*D.*sqrt((tau_star./tau_star_c)-1); %Use the Sklar and Dietrich empirical relation to calculate bedload layer height
if z_mixed < D;
    z_mixed = D; %the bedload layer height must be at least as high as the grain size
end

%% Kinematic eddy viscosity and Ld
eddy_visc = k1.*ustar.*lambda; %kinematic eddy viscosity (Eq. 23)
L_d = eddy_visc./w_net; %sediment diffusion length scale (Eq. 16)

%% Sediment source zone size
if h_total < lambda; %find depths in ZOFE
    b = 0.1.*lambda; % In ZOFE set source zone to a constant:
else
    b = 0.1.*h_total; % In ZOEF set source zone to scale with depth:
end
delta = 2.*b;

%% Back calcluate cb with integration over tailwater
h_tw_vec = linspace(h_pool,h_total,100); %set up a vector for water heights
% above z_lip and below the watersurface, we'll integrate over this to
% calculate the averaged sediment concentration at the pool lip

% Calculate the vertical decay of sediment concentration
vertical_decay = exp(-1.*((h_tw_vec-z_mixed)./L_d)); %calculate Eq. 24 for all depths from z_lip to z_water
i2 = h_tw_vec < z_mixed; %check if points in the tailwater are within the bedload layer
% this will only happen if the pool is very shallow.
vertical_decay(i2) = 1; %within the bed load layer, concentration does not vary with z

% Calcluate the radial decay of sediment concentration
k1_radial = 1./(besseli(0,delta./L_d)+(besseli(1,r_pool./L_d)...
    .*besselk(0,delta./L_d))./besselk(1,r_pool./L_d));
k2_radial = k1_radial.*besseli(1,r_pool./L_d)./besselk(1,r_pool./L_d);
radial_decay = k1_radial.*besseli(0,r_pool./L_d)+k2_radial.*besselk(0,r_pool./L_d);

% Take the depth-averaged concentration
dz_tw = mean(diff(h_tw_vec));

if h_total > z_mixed %if the pool deeper than the height of the mixed layer
    chi = sum(radial_decay.*vertical_decay.*dz_tw); %Scheingross and Lamb in review, Eq. 11
    Cb = (Qs./Q).*(h_TW./chi); % Scheingross and Lamb in review, Eq. 12
else
    Cb =  (Qs./Q)./radial_decay;
end

if r_pool < (2*b) %if the pool is narrower than the jet-descending region
    clear Cb
    Cb = (Qs.*h_TW)./(Q.*sum(vertical_decay).*dz_tw);
end

if 1/L_d == 0 %this is the case where w_up = 0 and particles are advected upward
    clear Cb
    Cb = Qs./Q;
end

if Cb > 0.2; %find excessively large (and physically unreasonable) values of near-bed sed. conc.
    Cb = 0.2; %Don't let near-bed sediment concentration exceed 0.2.
end

%% Define the floor and wall vectors
dr = r_pool./dr_spacing;
r_pos = linspace(dr,r_pool,dr_spacing);

dz = h_pool./dz_spacing;
h_pos = linspace(dz,h_pool,dz_spacing);


%% c_0(z)
c_0_z = Cb.*exp(-1.*((h_pos-z_mixed)./L_d)); %Eq. 21 setting z = hpool
i3 = h_pos < z_mixed; %find pool depths less than the bedload layer height
c_0_z(i3) = Cb; %within the bed load layer, keep concentration constant at cb level

%% Sediment concentration along pool floor
top_fraction = besseli(0,r_pos./L_d)+(besseli(1,r_pool./L_d)./...
    besselk(1,r_pool./L_d)).*besselk(0,r_pos./L_d);
bottom_fraction = besseli(0,2.*b./L_d)+(besseli(1,r_pool./L_d).*...
    besselk(0,2.*b./L_d))./besselk(1,r_pool./L_d);
c_bed_vs_r = Cb.*top_fraction./bottom_fraction; %

i4 = r_pos < (2.*b); %find radial distances within the sediment source zone
c_bed_vs_r(i4) = Cb; %keep concentration constant within source zone
i5 = 1./L_d == 0; % these correspond to cases when w_net = 0
c_bed_vs_r(i5) = Cb; %keep concentration constant when w_net = 0 (that is, particles have no gravitaitonal settling)

if 1/L_d == 0 %this is the case where w_up = 0 and particles are advected upward
    c_bed_vs_r(:) = Cb;
end

%% Sediment concentration along pool walls
top_fraction_walls = besseli(0,r_pool./L_d)+...
    (besseli(1,r_pool./L_d)./besselk(1,r_pool./L_d)).*besselk(0,r_pool./L_d);
c_wall_vs_z = c_0_z.*top_fraction_walls./bottom_fraction;

if length(c_wall_vs_z)==1
    jfklsdjfkldsj
end

if r_pool < (2*b); %find radial distances within the sediment source zone
    c_wall_vs_z(:) = Cb; %keep concentration constant within source zone
end
if (1/L_d) == 0; % these correspond to cases when w_net = 0
    c_wall_vs_z(:) = Cb; %keep concentration constant when w_net = 0 (that is, particles have no gravitaitonal settling)
elseif isinf(L_d) == 1
    c_wall_vs_z(:) = Cb; %keep concentration constant when w_net = 0 (that is, particles have no gravitaitonal settling)
end

%% In case Cb is too high from imposed sediment load
if Cb > 0.2; %find excessively large (and physically unreasonable) values of near-bed sed. conc.
    c_bed_vs_r(:) = NaN;
    c_wall_vs_z(:) = NaN;
end
Cb_normal = k2.*((tau_star./tau_star_c)-1).^1.5; % Eq. 25
if tau_star/tau_star_c < 1; %find cases where tau* < tau*c
    Cb_normal = 0; %if shear stress drops below critical, set entrainment to zero
end
if 0.99*Cb > Cb_normal %this means the sediment load is over the transport capacity
    c_bed_vs_r(:) = NaN;
    c_wall_vs_z(:) = NaN;
end

c_rpool_z = Cb_normal.*vertical_decay.*radial_decay; % (Eq. 27)
if r_pool < delta %if we're within the jet-descending region
    c_rpool_z(:) = Cb_normal.*vertical_decay; %make sure there's no radial decay
end
if (1/L_d) == 0; % these correspond to cases when w_net = 0
    c_rpool_z(:) = Cb_normal; %keep concentration constant when w_net = 0 (that is, particles have no gravitaitonal settling)
end
% Take the depth-averaged concentration to computer Qsc_pool
Qsc_pool = Q.*sum(c_rpool_z).*dz_tw./h_TW; %numerical integration of Eq. 29

if Qs > Qsc_pool
    c_bed_vs_r(:) = NaN;
    c_wall_vs_z(:) = NaN;
end

%% capacity for cover term
c_capacity_bed_vs_r = Cb_normal.*top_fraction./bottom_fraction; %

i4 = r_pos < (2.*b); %find radial distances within the sediment source zone
c_capacity_bed_vs_r(i4) = Cb_normal; %keep concentration constant within source zone
i5 = 1./L_d == 0; % these correspond to cases when w_net = 0
c_capacity_bed_vs_r(i5) = Cb_normal; %keep concentration constant when w_net = 0 (that is, particles have no gravitaitonal settling)

if 1/L_d == 0 %this is the case where w_up = 0 and particles are advected upward
    c_capacity_bed_vs_r(:) = Cb_normal;
end

end