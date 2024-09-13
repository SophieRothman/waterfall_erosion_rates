function [Etl_distd,  Htl, Wtl, Qflood, Qsc, qs]=...
    only_tl_sites(S, Q_gage, If,   A, A_gage, E_d, D_50, D_85, If_s, equ, qsqsc)%E_vertd = sum(E_vert_inst_arrayd.*Drat_d, 'all').*If; %Vertical plunge pool erosion rate
%(Slope, S_poold, Width, W_brink, Qg, If, h_brd, H_dropd, r_poold, A, A_gage, E_b, D);
%E_latd = sum(E_lat_inst_arrayd.*Drat_d, 'all').*If ;%Lateral plunge pool erosion rate
%Etl_brink_d_tt
%% Inputs

% S = 0.046; %slope of channel average
% S_poold = 0.1; % Slope leading into the lip
% W = 16; %channel width upstream of the waterfall (m) [measure in Google Earth]
% W_brink=10; %Width at wf lip
% Q_gage = 108; %Discharge of interest at gage
% h_brd =4.3; %pool depth to bedrock (m) [we don't know how deep the pools actually are, so better
%             %just to keep this set at a constant value]
% H_dropd =9.1; %waterfall drop height (m) [as with above, for a first order appproximation 
%                 %we can just hold waterfall drop height constant between
%                 %reaches]
% r_poold = 4.75; %Pool radius (m)
%A = 7.5e7; %draingae area of waterfall (m2)
%If = 0.01; %intermittancy factor
%% Basin specific constant inputs

%D = 0.01.*[.1, .3, 1, 3.5, 7.5, 15, 35, 60]; %Grain size (m) [maybe try to estimate this from photos? If you're
               %having problems estimatnig you can just leave it constant
               %for now]

%Drat_d=[0.066006601, 0.04620462, 0.290429043, 0.257425743, 0.148514851, 0.089108911, 0.075907591, 0.02640264];

               
Fr_brink=1;
 %Erosion rate (m/My) [Scott says 50 m/My is what we might expect in Dinkey,
            % and erosion rates might get as high as 1000 m/My in the
            % Kawaeh's, but probably ~100-500 m/My is more likely]

%E_d=35;   %Erosion rate basin average (m/Myr)      
%H=.54; %Flow depth (m)- Estimated for Dinkey with the ferguson equation (2007) - Sophie
loopsnum=length(D_50);
detail=200;


%% Gage Specific Constants
% this is tricky.  The best way to do this will be to use drainage area
% scaling for a reference discharge. I suggest you do something like:
% 1) Download gage data for the Kings or Kawaeh (or both if it exists). 
% 2) Make a rating curve and find a characteristic discharge at the gaging
% station (i.e., a 2 year return period flow) - I've included a script in
% this folder (RI_fct.m) that will calculate a two year and 10 year
% recurrence interval discharge if you feed it a vector of the peak
% streamflow values each year (you can get this off the USGS site - let me
% know if you have issues with it). 
% 3) Assume that discharge is proportional to drainage area to calculate
% the discharge at your site.


%A_gage = 1.002e9; %draingae area at the gage (m2) (note USGS reports this on their site)
% 63 ;%linspace(30,200,loopsnum); %(A./A_gage).*Q_gage; %this calculates a characteristic discharge at the site (m3/s)
Qflood =Q_gage*(A/A_gage);


%% Physics CONSTANTS (you probably won't need to change these)
sigmaT = 8e6; %rock tensile strength (Pa) [8 MPa is a reasonable tensile strength for granite]

R = 1.65; %normalized submerged sediment density (non-dimensional)
g = 9.81; %gravitational acceleration

%% INPUTS FOR CALCULATING FLOW HEIGHT IN NON-WATERFALL SCENARIO (FERGUSON 2007)

a1=6.5;   %%empirical constants
a2=2.5;    %%empirical constants

Hrange=linspace(.1, 7.5, detail);  %Array of all the flow heights we're inputting to see what Q correspondsVrange=NaN(1, detail);
Qrange=NaN(1, detail);
Vrange_brink=NaN(1, detail);
Qrange_brink=NaN(1, detail);
%H=zeros(1, loopsnum);   %Array of calculated flow depth for each Qflood -
%FOR MULltiple flood levels
%V=zeros(1, loopsnum);    %Array of calculated velocity for each Qflood

%% Calculations and figures of flow depth


%FOR MULTIPLE FLOOD LEVELS for m=1:detail
%     ushear=sqrt(g*Hrange(m)*S);
%     Vrange(m)=ushear*(a1*a2*(Hrange(m)/D))/(sqrt((a1^2)+((a2^2)*((Hrange(m)/D)^(5/3)))));
%     Qrange(m)=Vrange(m)*Hrange(m)*W;
%     for n=1:loopsnum   %For each height/discharge pair, we cycle through all the Qfloods to see what matches
%     if Qrange(m)<Qflood(n)   %we want to know the H when Q=61.8
%         H(n)=Hrange(m);
%         V(n)=Vrange(m);
%     end
%     end
% end
Davg=mean(D_50);
%using Davg for main channel tl erosion but min for lip to represent that
%lips tend to be sediment free

Widths=equ(1)*(Hrange).^equ(2);
Htl=NaN;
Wtl=NaN;
match=0;
for counter=1:detail
    if match==0
        W=Widths(counter);
        H_match=Hrange(counter);
        for m=1:detail
             ushear=sqrt(g*Hrange(m)*S);
             Vrange(m)=ushear*(a1*a2*(Hrange(m)/D_85))/(sqrt((a1^2)+((a2^2)*((Hrange(m)/D_85)^(5/3)))));
             Qrange(m)=Vrange(m).*Hrange(m).*W;
             for n=1:loopsnum   %For each height/discharge pair, we cycle through all the Qfloods to see what matches
                if Qrange(m)<Qflood   %we want to know the H when Q=61.8
                    H=Hrange(m);
                    V=Vrange(m);
                end
             end

        end
        if abs(H-Hrange(counter))<0.05
            match=1;
            Htl=H;
            Wtl=W;
        end
    end
end

%Finding Hbrink via Froude=1
%Hbrink=(Qflood.^(2/3))/((W_brink.^(2/3))*g.^(1/3))
%Hbrink=(Qflood./(W_brink.*(g^0.5)*Fr_brink)).^(2/3);
%figure

%legend_label1=string(zeros(1,loopsnum));
%legend_label2=string(zeros(1,loopsnum));

% plot(H, Qflood,  '*')
% hold on
% legend_label1=(string(H) + 'm, ' + string(Qflood) + ' m^3/s');
% %plot(Hrange, Qrange, 'k-', 'linewidth', 2)
% xlabel('Flow Depth (m)')
% ylabel('Discharge (m^3/sec)')
% legend(legend_label1, 'Location', 'northwest')
% hold off
%%%
% figure
% for n=1:loopsnum
%     plot(H(n), V(n), '*')
%     hold on
%     legend_label2(n)=(string(H(n)) + 'm, ' + string(V(n)) + ' m/s');
% end
% plot(Hrange, Vrange, 'k-', 'linewidth', 2)
% xlabel('Flow Depth (m)')
% ylabel('Velocity (m/sec)')
% legend( legend_label2, 'Location', 'northwest')
% hold off

%% Convert long term erosion rates into sediment flux AND Sed FLux Arrays
% We're going to do this by multiplying the erosion rate by the drainage 
% area, to get a volumetric sediment flux. 
If_s=1;
Qs_per_yeard = If_s*E_d.*A;  %The factor of 1e6 is convert from m/My to m/yr- canceled out by the conversion from km to m
% OK, now you have a sediment flux in m3/yr. This would be OK if the river
% was constantly transporting sediment, but in reality, the vast majority
% of the sediment transport will occur over a few days, say ~1% of the
% total year.  We're going to use the 1% value as an intermittancy factor
% to calculate a typical sediment flux at a waterfall during high flow.

If=.01;

qs_ed = (Qs_per_yeard*If_s)/If/(365*24*60*60); %Predicted Sediment flux at high flows based on Stock 2005 If is always .01 

%qb_ed=qs_ed./W;


%%REst of needed variables

tau_star=(Htl.*S)./(R.*D_50);

ustar = (tau_star.*R.*g.*D_50).^0.5;  %%Shear velocity

tau_star_c = 0.15*(S^0.25); %critical Shields stress
tstage=tau_star./tau_star_c;


%Calculating the bankful Sediment Load
Qsc=Wtl*sqrt(R.*g.*D_50.^3).*5.7.*((tau_star - tau_star_c).^(3/2)); % m^3/sec
%Using that value to create an array of possible sediment loads

%qt = 5.7.*(R.*g.*D.^3).^0.5.*(tau_star- tau_star_c).^(3/2); 
Qs=linspace(0.01.*Qsc(n), Qsc(n)/2, detail);
qb=Qsc./Wtl;%The TOtal load model requires sediment flux normalized by width
%qsc=qs_ed;
%qb_ed=qs_ed./W; % need to comment out
%qs=qs_ed; %ned to comment out


qb_ed=qb*qsqsc*(.01/If); %critical need to uncomment
qs=qb_ed*Wtl; %need to uncomment


%qb_ed=qsc/W;
%% CALCULATING total load erosion values ASSUMING DETACHMENT LIMITED, NO COVER EFFECTS
Etl_distd=zeros(1, loopsnum);


%for m = 1:length(qb)

% [E_1(m) E_st(m) E_Sklar(m)] = total_load_code_fct_public(D,S,qb(m),H,tau_star,ustar,tau_star_c,tstage);

[Etl_distd] = total_load_code_fct_public(D_50,S,qb_ed,Htl,tau_star,ustar,tau_star_c,tstage);
% you don't need the sklar and dietrich or viscously damped formulations
% for what you're doing, so no need to save those variables
    
%end


 %intermittancy factor of 1%, that is, we expect erosion to occur 
%Etl_distd=Etl_distd.*(Drat_d).*If   ;         % only 1% of the year
%Etl_distk=Etl_distk.*(Drat_k).*If []    ;       % only 1% of the year

%Etl_annual_d = sum(Etl_distd.*(Drat_d), 'all').*If ; % this is the expected erosion rate you'd get 
        % over the course of a year with a 1% intermittancy factor
        
        
%% Figures of Total Load Erosion
% figure
% legend_label3=string(zeros(1, 2));
% 
% loglog(D.*1000, Etl_distd)
% legend_label3(1)='Dinkey Creek Erosion';
% hold on
% loglog(D.*1000, Etl_distk)
% legend_label3(2)='MFK Erosion';    

% plot(qb, E_Sklar, 'b:', 'linewidth', 1)
% plot(qb, E_st, 'r--', 'linewidth', 1) 
%plot(ones(1000).*qb_stock, linspace(0,1,1000).*(E_annual_total_load(loopsnum, detail)), 'r--')
% The above line of code will plot a vertical line of the Stock erosion
% rate so you can compare it to your estimates
% xlabel('Sediment size(mm)')
% ylabel('Erosion Rate (mm/yr)')
% title('Total Load Model')
% legend(legend_label3, 'Location', 'northeastoutside')
% hold off
% %% running waterfall code once to fix pool radius, get initial brink and lip depths
% n=1;
% [E_vert_inst_arrayd(n), E_lat_inst_arrayd(n), h_brinkd(n), h_TWd(n), r_jetd(n)] = pool_erosion_for_sophie(qs_ed, r_poold, ...
%     Qflood, h_brd, H_dropd, W, D(1), S_pool, sigmaT);
% [E_vert_inst_arrayk(n), E_lat_inst_arrayk(n), h_brinkk(n), h_TWk(n), r_jetk(n)] = pool_erosion_for_sophie(qs_ed, r_poolk, ...
%     Qflood, h_brk, H_dropk, W, D(1), S_pool, sigmaT);
% %% Resetting r_pool
% r_poolk=1.6151;
% r_poold=1.6151;
% 


end