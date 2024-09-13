function [E_vert_avg_mm_yr, E_lat_avg_mm_yr, r_jet, h_brink, h_TW, lambda, h_sed_hp, lambda_hp, delta_rp, ...
    h_total, c_bed_avg, c_wall_avg,  w_vert_avg, u_wall_avg, ...
    Fe_bed_avg, ws_wup, code, h_sed] = ...
    pool_erosion_for_sophie_tw(Qs, r_pool, Q, h_pool, H, W, D, S, sigmaT, Htl)

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
% Qs: sediment supply from upstream (m3/s)
% r_pool: plunge-pool radius (m)
% Q: water discharge (m3/s)
% h_pool: plunge pool depth (m)
% H: waterfall drop height (m)
% W: width of channel upstream of the waterfall (m)
% D: grain size (m)
% S: average slope of the reach (m/m)
% sigma_t: rock tensile sterngth (Pa)

% Outputs
% E_vert_avg_mm_yr: Area-weighted average vertical erosion rate (mm/yr)
% E_lat_avg_mm_yr: Depth-averaged lateral erosion rate (mm/yr)
% All other parameters are outputs to check on model performance

%% Calculate the flow depth at the brink of the waterfall (don't
%worry about the details here)
g=9.81;
h_normal = ((Q.*(3.*D).^(1/6))./(8.1.*W.*sqrt(g.*S))).^(3/5);
u_n = Q./(h_normal.*W); %normal flow velocity upstream of the waterfall
Fr_n = u_n./sqrt(g.*h_normal); % Fr upstream of the waterfall
if Fr_n > 1
    Fr_lip = ((0.4+Fr_n.^2)./Fr_n.^2).^(3/2).*Fr_n; %Fr at the waterfall lip for supercritical flow
else
    Fr_lip = (1.4./Fr_n.^(2/3)).^(3/2).*Fr_n; %Fr at the waterfall lip for subcritical flow
end
h_brink = ((Q./W)./(sqrt(g).*Fr_lip)).^(2/3); %flow depth at the waterfall brink
%% Set some initial parameters
% Discritization of pool radius and height, for most purposes, using values
% of 100 will work, but finer discritization may be needed in some cases.
dr_spacing = 100;
dz_spacing = 100;

dr = r_pool./dr_spacing;
q = Q./W; %discharge per unit upstream channel width

%% First check to see if, based on the imposed sediment supply, the pool 
% should be filled with sediment or have exposed bedrock 
z = linspace(0.001, h_pool, dz_spacing); %set up a vector of pool depths
for j = 1:length(z)
    [Qsc_pool_vec(j), h_TW, notuse] = Qsc_pool_public_tw(r_pool, Q, z(j), H, W, D, h_brink, ...
        dr_spacing, dz_spacing, Htl); 
    %calculate sediment transport capacity at each pool depth following Scheingross and Lamb (2016)  (r_pool, ...
    %Q, h_pool, H, W, D, h_brink, dr_spacing, dz_spacing)
end
Qsc_Qs = Qsc_pool_vec./Qs; %ratio of sediment transport capacity to sediment supply

if sum(Qsc_Qs>1) == length(Qsc_Qs) % bedrock exposed on pool floor, supply limited conditions
    % the pool should be eroding bedrock
    h_sed = h_pool; % the sediment depth is the bedrock depth
    code = 1; % define a flag to identify this condition
elseif sum(Qsc_Qs<1) == length(Qsc_Qs) % if all depths are transport limited
    % the pool should be completely filled
    h_sed = 0;
    code = 2; % define a flag to identify this condition
else % there's partial alluvial fill
    % Find the level of partial alluvial fill
    [Qsc_Qs_uni, iA] = unique(Qsc_Qs);
    if sum(isnan(Qsc_Qs(iA)))==length(Qsc_Qs(iA))
        break_the_code % break the code everything is NaN
    end
    h_sed = interp1(Qsc_Qs(iA), z(iA), 1);
    code = 3; % define a flag to identify this condition
end

if h_sed < 7*D %if the pool has aggraded within 7 grain diameters of the downstream pool lip
    if h_sed == h_pool %there's no alluvium
        code = 1;
    else
        code = 4; %this is outside of the model domain
    end
end

%% Second, calculate sediment concentrations based on imposed sediment supply
[c_bed_vs_r, c_wall_vs_z, r_pos, h_pos, w_net, b, u_jet_impact, ...
    ws, Cdrag, lambda, r_jet, w_up, c_capacity_bed_vs_r, u_jet, h_total] = ...
    c_from_Qs_fct_tw(Qs, r_pool, Q, h_sed, H, W, D, h_brink, dr_spacing, dz_spacing, Htl);


%% Third, step through all 4 possible cases, and calculate erosion rates in each
if code == 1 % the pool is eroding
    % Calculate erosion rates as function of (r) and (z)
    [E_vs_r, E_vs_z, u_wall_avg, w_vert_cubed, j1, c_wall_avg, cover_bed] = ...
        pool_erosion_fct_public_tw(c_wall_vs_z, c_bed_vs_r, w_net, ...
        h_pool, H, D, b, r_pos, u_jet_impact, Cdrag, ...
        sigmaT, c_capacity_bed_vs_r, u_jet, r_jet, h_pos, ...
        q, h_brink, Qs, dr_spacing, Htl);
    E_lat_avg = nanmean(E_vs_z); %average lateral erosion (Eq. 6b)
    E_vert_avg = sum(E_vs_r.*2.*pi.*r_pos.*dr)./(pi.*r_pool.^2); %average vertical erosion (Eq. 6a)
    Fe_bed_avg = sum(cover_bed.*2.*pi.*r_pos.*dr)./(pi.*r_pool.^2);
elseif code == 2 % the pool is completely filled
    E_vert_avg = 0;
    E_lat_avg = 0;
    w_vert_cubed = NaN;
    u_wall_avg = NaN;
    c_wall_avg = NaN;
    Fe_bed_avg = 0;
elseif code == 3 % the pool is partially alluviated
    % Use Scheingross and Lamb (2016) to calculate sediment concentrations
    % in an alluvial pool
    [blah,  h_TW, c_wall_vs_z, blah2, blah3, c_bed_vs_r] = ...
        Qsc_pool_public_tw(r_pool, Q, h_sed, H, W, D, h_brink, dr_spacing, dz_spacing, Htl);
    c_capacity_bed_vs_r = c_bed_vs_r;
    [E_vs_r, E_vs_z, h_TW, u_wall_avg, w_vert_cubed, j1, c_wall_avg, cover_bed] = ...
        pool_erosion_fct_public_tw(c_wall_vs_z, c_bed_vs_r, w_net, ...
        h_sed, H, D, b, r_pos, u_jet_impact, Cdrag, sigmaT, ...
        c_capacity_bed_vs_r, u_jet, r_jet, h_pos, ...
        q, h_brink, Qs, dr_spacing, Htl);
    E_lat_avg = nanmean(E_vs_z).*h_sed./h_pool; %average lateral erosion (Eq. 6b)
    E_vert_avg = 0; % no vertical erosion since the pool is alluviated
    Fe_bed_avg = 0; % no bedrock exposed
elseif code == 4 % the pool has aggraded within 7D of the downstream lip, this is outside of the model domain
    E_vert_avg = 0;
    E_lat_avg = 0;
    w_vert_cubed = NaN;
    u_wall_avg = NaN;
    c_wall_avg = NaN;
    Fe_bed_avg = 0;
else
    break_the_code % break the code because we should always have code 1, 2, 3, or 4
end

%% Convert erosion rates from (m/s) to (mm/yr)
E_vert_avg_mm_yr = E_vert_avg.*60*60*24*365.*1000;
E_lat_avg_mm_yr = E_lat_avg.*60*60*24*365.*1000;

%% Add some extras to explore model predictions
h_sed_hp = h_sed./h_pool;
lambda_hp = lambda./h_total;%note I'm doing this on the sed layer with tailwater
delta_rp = (2.*b)./r_pool;
w_vert_avg = (sum(w_vert_cubed.*2.*pi.*r_pos.*dr)./(pi.*r_pool.^2))^(1/3);
ws_wup = ws./w_up;
c_bed_avg = sum(c_bed_vs_r.*2.*dr.*r_pos)./(r_pool.^2);
if code == 2
    w_vert_avg = NaN;
    c_bed_avg = NaN;
    u_wall_avg = NaN;
    c_wall_avg = NaN;
end
if code == 3
    w_vert_avg = NaN;
    c_bed_avg = NaN;
end

end