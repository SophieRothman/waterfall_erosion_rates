function [E_vs_r, E_vs_z,  u_wall_avg, w_vert_cubed, w_vert, c_wall_avg, Fe_bed, ...
    u_wall_of_z, E_wall_factor, c_impacting_avg, c_settle_avg, w_fall_3, w_settle_3, ...
    w_fall_center, c_impacting_center, w_fall_delta, c_impacting_delta] = ...
    pool_erosion_fct_public(c_wall_vs_z, c_bed_vs_r, w_net, h_pool, H, D, b, r_pos, ...
    u_jet_impact, Cdrag, sigma_t, c_capacity_bed_vs_r, u_jet, r_jet, h_pos, ...
    q, h_brink, Qs, dr_spacing, Htl)
dbstop if error
%% Notes
% This MATLAB code calculates erosion in waterfall plunge pools following
% the theory of Scheingross and Lamb (in review). This code was created to
% produce the figures included in Scheingross and Lamb (in review),
% depending on your purposes, you may need to modify aspects of the code
% (e.g. the dz and dr spacing) for your specific purposes.  You are
% responsible for inspecting the code and making any necessary
% modifications for your own purposes.  Please report errors to Joel at
% jscheingross . a . t . gmail.com
% Code written by J. Scheingross


%% Inputs and outputs
%Inputs
% c_wall_vs_z: sediment concentration on the wall
% c_bed_vs_r: sediment concentration along the bed
% h_pool: plunge pool depth (m)
% H: waterfall drop height (m)
% D: grain size (m)
% b: half-width of waterfall jet (m)
% r_pos: a vector containing positions moving along the pool radius
% u_jet_impact: waterfall jet velocity on the plunge-pool floor (m/s)
% Cdrag: drag coefficient
% sigma_t: rock tensile strength (Pa)
% c_capacity_bed_vs_r: concentration of near-bed sediment as a function of
% radius at transport capacity
% u_jet: jet impact velocity at impact with the pool water surface (m/s)
% r_jet: waterfall jet radius at impact with pool water surface (m)
% h_pos: a vector containing positions moving along the pool walls
% q: water discharge per unit upstream channel width (m2/s)
% h_brink: water depth at the upstream waterfall brink (m)
% Qs: sediment supply (m3/s)
% dr_spacing: number of discritization elements for the plunge-pool floor

%Outputs
% E_vs_r: Vertical erosion on the bed (m/s)
% E_vs_z: Lateral erosion on the walls (m/s)
% all other outputs are quantities saved from the model which can be used
% to examine its performance

%% Constants
g = 9.81; %gravitational acceleration (m/s2)
nu = 1e-6; %kinmatic viscosity of water
rho_s = 2650; %sediment density
rho_w = 1000; %water density
noerosion=0;
%% Calculate tailwater depth (z_water - z_lip)
h_TW = Htl;%(q./sqrt(g)).^(2/3); %tailwater depth (Scheingross and Lamb, 2016, JGR, Eq. 11)
h_total = h_pool + h_TW; % this is the total water depth (z_water - z_sed)

%% Define zmixed zone
Cf_pool = 0.001; % assume a constant pool friction factor (Scheingross and Lamb, 2016, JGR)
tau = Cf_pool.*rho_w.*(u_jet_impact.^2); %plunge pool basal shear stress (Scheingross and Lamb, 2016, JGR, Eq. 1)
tau_star = tau./(1650*g.*D); % plunge pool basal Shields stress
tau_star_c = 0.045;
zmixed = 1.44.*D.*sqrt((tau_star./tau_star_c)-1); %Use the Sklar and Dietrich empirical relation to calculate bedload layer height after Scheingross and Lamb, 2016, JGR
if (tau_star/tau_star_c) < 1
    zmixed = NaN;
    break_the_code % break the code if we're below the threshold for motion
    noerosion=1;
end
if zmixed < D;
    zmixed = D; %the bedload layer height must be at least as high as the grain size
end

%% Lateral wall impacts
u_at_wall = 2.03.*u_jet.*r_jet./max(r_pos); %Eq. (28)
b_radial = 0.09.*max(r_pos); %Eq. (29)
u_wall_of_z = u_at_wall.*1.48.*(h_pos./b_radial).^(1/7).*(1-erf(0.68.*h_pos./b_radial)); %Eq. (27)

%%% Include only impacts by grains in the mixed layer
i_mixed = h_pos < zmixed; % index defining points in the mixed layer
u_wall_avg = 0.6.*(mean(u_wall_of_z(i_mixed).^3)).^(1/3); %Eq. (30)

u_wall_of_z(i_mixed) = u_wall_avg; % Make all impacts in the mixed layer at the same saltation velocity
u_wall_of_z(~i_mixed) = 0; %no impacts above the mixed layer

%%% Include a stokes threshold for viscous dampening
Stokes_thresh = 75; %set the threshold for Stokes dampening following Scheingross et al, 2014, Geology
u_wall_of_z((D.*u_wall_of_z.*rho_s./(9.*nu.*rho_w))<Stokes_thresh) = 0; %Set impact velocity to zero for Stokes # < stokes_threshold

E_wall_factor = (u_wall_of_z).^3.*c_wall_vs_z; %this is cwall and ulat^3 in Eq. (5b)
c_wall_avg = nanmean(c_wall_vs_z(i_mixed)); %average concentration along the wall in the mixed layer

%% Particle veritcal impact velocities 
%% First get impact velocities with an explicit finite-difference scheme
%(this is only for grains falling from the top of the waterfall)

%%%% We'll need a few things to start out with %%%%
r = D/2; %grain radius
particle_A = pi*r^2; %particle XS area
particle_V = (4/3)*pi*r^3; %particle volume
Cd = Cdrag;
rho_w = rho_w.*.7; %multiple by 0.7 to account for aeration of waterfall jets

% For simplicity, we group some constant following Lamb et al (2007)
C1w = (rho_s - rho_w)*(1/rho_s)*g;
C2w = 0.5*Cd*(rho_w/rho_s)*(particle_A/particle_V);

V_brink = q./h_brink; % velocity at waterfall brink assuming conservation of mass
beta = atan(sqrt(2.*g.*H)./V_brink); %angle of jet impact with water surface (Eq. 3 from Scheingross and Lamb 2016)
Cd_jet = 2.6; % Jet diffusion constant following Scheingross and Lamb 2016
lambda = 2.*Cd_jet.^2.*r_jet.*sin(beta); %length scale of ZOFE, Scheingross and Lamb 2016, Eq. 2
%%%% OK, let's keep going %%%%%

% First solve for velocity from the waterfall brink to the plunge-pool
% water surface
dz_drop = 10^-3; %define a length scale over which to discritize the fall through the air
% the theory isn't so sensitive to this, so a coarse length scale is
% generally OK, but be careful...
H_of_z = dz_drop:dz_drop:H; %set up a vector of positions from the waterfall brink to the pool water surface
ujet_of_z = sqrt(V_brink.^2 + 2.*g.*H_of_z);% jet velocity at impact (Eq. 5, Scheingross and Lamb 2016) following Stein and Julien 1993 conservation
wp_air = zeros(size(H_of_z)); %set up a vector to store particle velocity while following through the air
wp_air(1) = V_brink./10; %define an initial velocity at the waterfall brink.  Velocities of 0 had numerical stability issues, so we start with a small, non-zero velocity
for i = 2:length(H_of_z) % run through the explicit finite-difference scheme
    if wp_air(i-1) > ujet_of_z(i-1)
        dwp_air = (dz_drop./wp_air(i-1)).*(C1w-C2w.*(wp_air(i-1)-ujet_of_z(i-1)).^2);
    else
        dwp_air = (dz_drop./wp_air(i-1)).*(C1w+C2w.*(wp_air(i-1)-ujet_of_z(i-1)).^2);
    end
    wp_air(i) = wp_air(i-1)+dwp_air;
end
if isnan(wp_air(end)) == 1 % run some tests to see if I need this
    break_the_code % if there's a nan here, then something went wrong
end

% Second solve for particle velocity falling within the plunge pool
dz = 10^-4; %define the length scale over which to discritize the water column
z = dz:dz:h_total; %set up a vector for the water column
dr = max(r_pos)./dr_spacing; %define the length scale over which to discritize the pool radius

% We'll need to know the jet velocity everywhere in the plunge pool, so the
% next few lines take care of this
lam_hpool = lambda./z;
i1 = lam_hpool > 1 ; %find depths within ZOFE
u_of_z = zeros(size(z));
u_of_z(i1) = u_jet; % Jet impact velocity on the pool floor (Scheingross and Lamb, 2016 Eq. 6a)
u_of_z(~i1) = u_jet.*sqrt(lambda./z(~i1));% Jet impact velocity on the pool floor (Scheingross and Lamb, 2016 Eq. 6b)
b_of_z(~i1) = 0.1.*z(~i1); %jet-half width (i.e., half of the jet-descending region radius, Scheingross and Lamb, 2016, Eq. 13)
b_of_z(i1) = 0.1.*lambda; %jet-half width (i.e., half of the jet-descending region radius, Scheingross and Lamb, 2016, Eq. 13)

w_impact_bed_vector = zeros(size(r_pos)); % set up a vector to store impact velocities within the plunge pool
wp_pool = zeros(size(z));
wp_pool(1) = wp_air(end); %set the initial impact velocity at the water surface to the velocity solved for above for falling from the brink to water surface

% Run through the finite difference scheme
for i=2:length(z)
    r_max = 2.*b_of_z(i);
    i_r = r_pos < r_max;
    ujet_of_r = u_of_z(i-1).*exp(-0.69.*(r_pos(i_r)./b_of_z(i-1)).^2); %this is the quantity inside the integral in Eq. 23
    ujet_r_avg = 2.*sum(ujet_of_r.*r_pos(i_r).*dr)./(max(r_pos(i_r)).^2); % Eq. (23)
    if wp_pool(i-1) > ujet_r_avg
        dwp_zofe = (dz_drop./wp_pool(i-1)).*(C1w-C2w.*(wp_pool(i-1)-ujet_r_avg).^2);
    else
        dwp_zofe = (dz_drop./wp_pool(i-1)).*(C1w+C2w.*(wp_pool(i-1)-ujet_r_avg).^2);
    end
    wp_pool(i) = wp_pool(i-1)+dwp_zofe;
end
w_impact_bed_vector(:) = wp_pool(end);
i_r2 = r_pos > (2.*b); % fine positions outside the jet descending region
w_impact_bed_vector(i_r2) = 0; %there are no impacts of grains falling from the waterfall brink outside the jet-descending region

%% Second, average impact velocities across pool floor and account for grains falling out of suspension
w_settling = zeros(size(w_impact_bed_vector)); %set up a vector to store velocity of grains falling out of suspension within the pool
w_settling(:) = w_net; %set the settling velocity to the net settling velocity following Scheingross and Lamb (2016)

% Because erosion depends on velocity^3, it's useful to have these impact
% velocities in cubed form
bed_impact_cubed_settling = w_settling.^3;
bed_impact_cubed_fall = w_impact_bed_vector.^3;

%%%% Account for ratio of grains falling from the ewaterfall brink vs.
%%%% settling out of suspesnion
ind1 = r_pos<(2*b); %find positions inside the jet-descending region
if (2*b < max(r_pos)) %if the pool radius is larger than the radius of the jet descending region
    c_wf = Qs./(bed_impact_cubed_fall.^(1/3).*pi.*(2.*b).^2); %Eq. 15
else %if the pool radius is smaller than the radius of the jet descending region
    c_wf = Qs./(bed_impact_cubed_fall.^(1/3).*pi.*max(r_pos).^2); %Eq. 15
end
c_wf(isinf(c_wf)) = 0; %replace infinite c_wf values with 0, this comes from having velocities of zero outside of the jet-descending region
c_impacting_avg = nanmean(c_wf(ind1)); % average velocity of impacting grains

w_vert = (bed_impact_cubed_fall.*(c_wf./c_bed_vs_r)+bed_impact_cubed_settling.*(1-(c_wf./c_bed_vs_r))).^(1/3); %Eq. 14a
w_vert((D.*w_vert.*rho_s./(9.*nu.*rho_w))<Stokes_thresh) = 0; %Set impact velocity to zero for Stokes # < stokes_threshold

E_bed_factor = w_vert.^3.*c_bed_vs_r; %this is cbed x wvert^3 in Eq. (5a)
E_bed_factor(r_pos>(2.*b)) = bed_impact_cubed_settling(r_pos>(2.*b)).*c_bed_vs_r(r_pos>(2.*b)); %make sure outside of the jet-descending region all impacts are from particles settling out of suspension only

% Calculate a few extra quantitites just to see how the theory is working
w_vert_cubed = w_vert.^3;
w_fall_3 = nanmean(bed_impact_cubed_fall(ind1));
w_settle_3 = bed_impact_cubed_settling(1);
c_settle = c_bed_vs_r.*(1-c_wf./c_bed_vs_r);
c_settle_avg = nanmean(c_settle(ind1));
w_fall_center = (bed_impact_cubed_fall(1)).^(1/3);
c_impacting_center = c_wf(1);
i2 = find(r_pos < (2*b), 1, 'last' );
w_fall_delta = (bed_impact_cubed_fall(i2)).^(1/3);
c_impacting_delta = c_wf(i2);


%% Calcualte erosion rates
% A few constants following Lamb et al. 2008
A1 = 0.5;
Y = 5e10; %generic young's modulus (Pa)
kv = 1e6; %Single grain abrasion mill kv value
Eros_const = A1.*rho_s.*Y./(kv.*sigma_t.^2); %note that the quantiy Y / kv = ky = 0.05 MPa as used in the text

% Fraction of exposed bedrock
Fe_bed = 1 - c_bed_vs_r./c_capacity_bed_vs_r; %Eq. (24a)
Fe_walls = 1; %no cover on walls 

E_vs_z = Eros_const.*E_wall_factor.*Fe_walls; %lateral erosion along the pool walls [m/s]
E_vs_r = Eros_const.*E_bed_factor.*Fe_bed; %vertical erosion along the pool floor [m/s]
if noerosion
    E_vs_z = 0; %lateral erosion along the pool walls [m/s]
    E_vs_r = 0; %vertical erosion along the pool floor [m/s]
end

end