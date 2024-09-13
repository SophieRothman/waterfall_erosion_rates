function [E_vert_inst_arrayd, E_lat_inst_arrayd, ...
    Etl_distd,  Htl, h_brink, h_TW, If, Qflood, r_jetd, Wtl, Wbrink, qs_ed, lambda]=...
    all_erosion_sites(S, S_t,  Q_gage, If, h_brd, H_dropd, r_poold, A, A_gage, E_d, D_50, D_85, If_s, equ, qs)%E_vertd = sum(E_vert_inst_arrayd.*Drat_d, 'all').*If; %Vertical plunge pool erosion rate
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
If = 0.01; %intermittancy factor
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
e_basind=035;   %Erosion rate basin average (mm/yr)
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
Qflood =Q_gage*(A/A_gage);% 63 ;%linspace(30,200,loopsnum); %(A./A_gage).*Q_gage; %this calculates a characteristic discharge at the site (m3/s)



%% Physics CONSTANTS (you probably won't need to change these)
sigmaT = 8e6; %rock tensile strength (Pa) [8 MPa is a reasonable tensile strength for granite]

R = 1.65; %normalized submerged sediment density (non-dimensional)
g = 9.81; %gravitational acceleration

%% INPUTS FOR CALCULATING FLOW HEIGHT IN NON-WATERFALL SCENARIO (FERGUSON 2007)

a1=6.5;   %%empirical constants
a2=2.5;    %%empirical constants

Hrange=linspace(.1, 7.5, detail);  %Array of all the flow heights we're inputting to see what Q corresponds
Vrange=NaN(1, detail);
Qrange=NaN(1, detail);
Vrange_brink=NaN(1, detail);
Qrange_brink=NaN(1, detail);
%H=zeros(1, loopsnum);   %Array of calculated flow depth for each Qflood -
%FOR MULltiple flood levels
%V=zeros(1, loopsnum);    %Array of calculated velocity for each Qflood

%% Calculations and figures of flow depth
%THIS IS FOR total load calculations

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
Davg=D_50;
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

Qs_per_yeard = E_d.*1e6.*A./1e6;  %UNITS = m^3/year The factor of 1e6 is convert from m/My to m/yr and also to convert km^2 to m^2 If_s is what year storm it is 
% OK, now you have a sediment flux in m3/yr. This would be OK if the river
% was constantly transporting sediment, but in reality, the vast majority
% of the sediment transport will occur over a few days, say ~1% of the
% total year.  We're going to use the 1% value as an intermittancy factor
% to calculate a typical sediment flux at a waterfall during high flow.


If=.01;
qs_ed = qs;%(Qs_per_yeard*If_s)/If/(365*24*60*60); %Predicted Sediment flux at high flows based on Stock 2005 If is always .01 

qb_ed=qs_ed./Wtl;


%%REst of needed variables

tau_star=(H.*S)./(R.*D_50);
ustar = (tau_star.*R.*g.*D_50).^0.5;  %%Shear velocity
tau_star_c = 0.15*(S^0.25); %critical Shields stress
tstage=tau_star./tau_star_c;

%Calculating the bankful Sediment Load
Qsc=Wtl*sqrt(R.*g.*D_50.^3).*5.7.*(tau_star - tau_star_c).^(3/2); % m^3/sec
%Using that value to create an array of possible sediment loads


Qs=linspace(0.01.*Qsc(n), Qsc(n)/2, detail);
qb=Qsc./(100*Wtl);%The TOtal load model requires sediment flux normalized by width; 50% of Qsc
qstot=qs_ed;
%qb_ed=qsctot/(W);%qb;
%qs_ed=qb_ed*W;
%% CALCULATING total load erosion values ASSUMING DETACHMENT LIMITED, NO COVER EFFECTS
Etl_distd=zeros(1, loopsnum);
Etl_distk=zeros(1,loopsnum);
for n =1
%for m = 1:length(qb)

% [E_1(m) E_st(m) E_Sklar(m)] = total_load_code_fct_public(D,S,qb(m),H,tau_star,ustar,tau_star_c,tstage);

[Etl_distd(n)] = total_load_code_fct_public(D_50,S,qb_ed,H,tau_star(n),ustar(n),tau_star_c,tstage(n));
% you don't need the sklar and dietrich or viscously damped formulations
% for what you're doing, so no need to save those variables

%end
end

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
%% Calculations and figures of flow depth
%THIS IS FOR flow acceleration towards the lip

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
Davg=D_50;
%using Davg for main channel tl erosion but min for lip to represent that
%lips tend to be sediment free
% from ferguson (2000)
%from equation 10 qstar=q/(g^.5 D^(3/2))

%%
Sb=0.1;  %increase the slope on the lip
%Widths=equ(1)*Hrange.^equ(2);
match=0;
D_85=.001; %reduce the roughness on the lip
for counter=1:detail
    if match==0
        W=Widths(counter);
        for m=1:detail
             ushear=sqrt(g*Hrange(m)*Sb);

             Vrange(m)=ushear*(a1*a2*(Hrange(m)/D_85))/(sqrt((a1^2)+((a2^2)*((Hrange(m)/D_85)^(5/3)))));
             Qrange(m)=Vrange(m).*Hrange(m).*W;
             for n=1:loopsnum   %For each height/discharge pair, we cycle through all the Qfloods to see what matches
                if Qrange(m)<Qflood   %we want to know the H when Q=61.8
                    H=Hrange(m);
                    V=Vrange(m);
                end
             end

        end
        if abs(H-Hrange(counter))<0.2
            match=1;
            Wbrink=W;
            Hbrink=H;
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

%% Calculating Waterfall Vertical and Lateral erosion rates in units mm/yr
% Note, even those these values are in mm/yr, they're instanteous rates, 
% that is, they assume that the sediment flux you input is happening over 
% the entire year.
%make a loop to see how it evolves
% tff=1e6;
% dt=10;
% axis=1:dt:tff;
% tf=length(1:dt:tff);
% r_poola_d=zeros(1, tf);
% r_poola_d(1)=r_poold;
% h_bra_d=zeros(1, tf);
% h_bra_d(1)=h_brd;
% H_dropa_d=zeros(1,tf);
% H_dropa_d(1)=H_dropd;
% E_vertk=zeros(1, tf);
% E_latk=zeros(1, tf);
% E_vertd=zeros(1, tf);
% E_latd=zeros(1, tf);
% Sbrink_d=zeros(1, tf);
% Slip_k=zeros(1, tf);
% 
% r_poola_k=zeros(1, tf);
% r_poola_k(1)=r_poolk;
% h_bra_k=zeros(1, tf);
% h_bra_k(1)=h_brk;
% H_dropa_k=zeros(1,tf);
% H_dropa_k(1)=H_dropk;
% aller_wfdink=zeros(tf, loopsnum);
% aller_wfkaw=zeros(tf, loopsnum);
% aller_brink_dink=zeros(tf, loopsnum);
% aller_brink_kaw=zeros(tf, loopsnum);
% aller_lip_dink=zeros(tf, loopsnum);
% aller_lip_kaw=zeros(tf, loopsnum);
% Etl_brink_k=zeros(1, tf);
% Etl_brink_d=zeros(1, tf);
% Etl_lip_k=zeros(1, tf);
% Etl_lip_d=zeros(1, tf);




E_vert_inst_arrayd=zeros(1,loopsnum);
E_lat_inst_arrayd=zeros(1, loopsnum);

    for n=1
    [E_vert_inst_arrayd(n), E_lat_inst_arrayd(n), r_jetd(n), h_brink, h_TW, lambda] = pool_erosion_for_sophie_tw(qs_ed, r_poold, ...(Qs, r_pool, Q, h_pool, H, W, D, S, sigmaT)
        Qflood, h_brd, H_dropd, Wtl, D_50, S, sigmaT, Htl);

    %[Etl_brink_d(n)] = total_load_code_fct_public(D_50,S_poold,qb_ed,Hbrink,tau_star_brink(n),ustar_brink(n),tau_star_c,tstage_brink(n));
%     [Etl_lip_k(n)] = total_load_code_fct_public(D(n),Slip_k(i),qb_ed,H,tau_stars(i),ustar(i),tau_star_c(i),tstage(i));
%     [Etl_lip_d(n)] = total_load_code_fct_public(D(n),Slip_d(i),qb_ek,H,tau_star(i),ustar(i),tau_star_c(i),tstage(i));
    
    end
r_jetd=r_jetd(1);
    % You can use the Intermittency factor a calculate what the average erosion
    % rate is over the entire year (taking into account periods of sediment
    % transport and erosion and periods of no erosion)

    %Summing accross grain sizes and making vertical erosion zero if any
    %grain sizes aren't moving

%E_vertd = sum(E_vert_inst_arrayd.*Drat_d, 'all').*If; %Vertical plunge pool erosion rate

%E_latd = sum(E_lat_inst_arrayd.*Drat_d, 'all').*If ;%Lateral plunge pool erosion rate
%Etl_brink_d_ttl=sum(Etl_brink_d.*Drat_d, 'all').*If ; %total erosion at brink

    
    %Recording erosion values
% aller_wfdink(i)=E_vert_inst_arrayk;
% aller_lat_wfdink(i)=E_vert_inst_arrayk;
% aller_wfkaw(i)=E_lat_inst_arrayk;
% aller_lat_wfkaw(i)=E_lat_inst_arrayk;
% aller_brink_dink(i)=zeros(tf, loopsnum);
% aller_brink_kaw(i)=zeros(tf, loopsnum);
% aller_lip_dink(i)=zeros(tf, loopsnum);
% aller_lip_kaw(i)=zeros(tf, loopsnum);

%Updating variables
% Slip_d(i+1)=
% Slip_k(i+1)=
% 
% Sbrink_d(i+1)=
% Sbring_k(i+1)=
% 
% h_bra_k(i+1)=h_bra_k(i)+ (E_vert_mm_yrk)*.5*dt/1000;
% h_bra_d(i+1)=h_bra_d(i)+ (E_vert_mm_yrd)*.5*dt/1000;
% 
% H_dropa_d(i+1)=H_dropa_d(i)+((E_vert_mm_yrd)*.5*dt/1000)+(E_lat_mm_yrd*dt/1000)*S_pool;%((E_vert_mm_yrd-Etl_annual_d)*.5*dt/1000)+(E_lat_mm_yrd*dt/1000)*S_pool;
% H_dropa_k(i+1)=H_dropa_k(i)+((E_vert_mm_yrk)*.5*dt/1000)+(E_lat_mm_yrk*dt/1000)*S_pool;%((E_vert_mm_yrk-Etl_annual_k)*.5*dt/1000)+(E_lat_mm_yrk*dt/1000)*S_pool;
% 
% r_poola_d(i+1)=r_poola_d(i)+E_lat_mm_yrd*dt/1000;
% r_poola_k(i+1)=r_poola_k(i)+ E_lat_mm_yrk*dt/1000;



%(note, I included vertical and lateral erosion here, but all you should
%really care about for this exercise is vertical erosion)
end