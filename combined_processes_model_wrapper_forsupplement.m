%clc, close all, clear all
%This page written by Sophie Rothman
%However, it draws on three codes (which call on numerous other codes)
%which are all written by Joel Scheingross
%
%% load parameters
data = readmatrix('D:\Users\srothman\Documents\field_materials\persite_pred_vlidar.xlsx');
width_equ  = readmatrix('D:\Users\srothman\Documents\field_materials\width_power3.xlsx');
%data_control=readmatrix()
%load('E:\sitewide_parameters.mat')
%using the 2 and 10 year discharges from dinkey and mfk

% flood parameters

%Qg_k=1xn array of floods of different magnitudes at the MFK at the gage, for example [182.6, 619.34, 949.7, 1142.96, 1386.4];
%Qg_d=1xn array of floods of different magnitudes at Dinkey Creek at the gage, for example [ 44.6, 123.74, 183.6, 218.622, 262.7];
%Ifyearly=1xn array if intermittancy factor indicating what fraction of the
    %year the flood event happens over for examPLE [0.01, 0.01, 0.01, 0.01,
    %0.01]
%If_sed=1xn array indicating the return period of each flood, for example [2, 5, 10, 15, 25];
%A_gaged = 131.3125; %draingae area at the Dinkey Creek gage (km2) (note USGS reports this on their site)
%A_gagek=1082.62; drainage area at the MFK gage (km2)

%per site parameters 

% Names=mx1 array of site names, for example
% {'MFK21_CRN_06'
% 'MFK21_CRN_04'
% 'MFK21_CRN_11'
% 'MFK22_CRN_02'
% 'MFK21_CRN_09'
% 'DNK21_CRN_06'
% 'DNK22_CRN_01'
% 'DNK22_CRN_02'
% 'DNK21_CRN_0304'
% 'DNK22_CRN_03'};



%width_equ=mx2 array to indicate the power law width equation for each
    %site, in which columns house a and b values for each site, where
    %width=a*depth^b
% qsqsctl_10=mx1 array of Qs/Qsc (sediment transported over transport capacity)
%for all sites equal to that of the site in a 10 year flood -for example we
%use
% [0.130552909865198
% 0.331485570313698
% 0.230813975763316
% 0.146033901286548
% 0.382395606570716
% 0.0289392864621729
% 0.0110829207907272
% 5.41802238019842e-05
% 0.00113872735573026
% 0.0139815834453041];

%Q_gage =mxn array of all the flood sizes for each site for exmaple
%   [Qg_k
%     Qg_k
%     Qg_k
%     Qg_k
%     Qg_k
%     Qg_d
%     Qg_d
%     Qg_d
%     Qg_d
%     Qg_d];

% A_gaget=nx1 array of of gage drainage area for example
%   [A_gagek
%     A_gagek
%     A_gagek
%     A_gagek
%     A_gagek
%     A_gaged
%     A_gaged
%     A_gaged
%     A_gaged
%     A_gaged];
%Needed data


%S_tot=mx1 array of total reach slope
%S_nowf=mx1 array of S_nonwf values per site (i.e. slope of planar reach
%between waterfalls)
% D85=mx1 array of average D84 size per river, in m
% D50=mx1 array of average D50 size per river, in m
% da= mx1 array of drainage area per site, in km^2
% dh= mx1 array of reach-averaged drop height per site, in m
% depth= mx1 array of reach-averaged plunge pool depth per site, in m
% radius= mx1 array of reach-averaged plunge pool radius per site, in m
% wf_freq=mx1 array of waterfall frequency 
% e_ba=mx1 array of basin-averaged erosion rate in the given river, for
% sediment supply, in m/million years
% length_np=mx1 array of channel length without plunge-pools for each reach (m)
% length_ch=mx1 array of channel length for each reach (m)
% S_brink=S_tot%

%% making



E_vert_inst_arrayd=NaN( length(Names), length(Qg_k));
E_lat_inst_arrayd=NaN( length(Names), length(Qg_k));
E_lat_inst_arrayd=NaN( length(Names), length(Qg_k));
E_vert_inst_array=NaN( length(Names), length(Qg_k));
Etl_brink_d=NaN( length(Names), length(Qg_k));

Etl_distd=NaN(length(Names), length(Qg_k));
Etl_nowf=NaN(length(Names), length(Qg_k));
Qflood=NaN( length(Names), length(Qg_k));
Qfloodtl=NaN( length(Names), length(Qg_k));
Htl=NaN( length(Names), length(Qg_k));
Htl_nowf=NaN( length(Names), length(Qg_k));
Wtl_nowf=NaN( length(Names), length(Qg_k));
Qflood_nowftl=NaN( length(Names), length(Qg_k));
Qsc_pool=NaN( length(Names), length(Qg_k));
Names=categorical(Names);

%% Setting up function parameters - Here we want to calculate erosion at each waterfall and then 
%NOW USING tl width as Tailwater depth

% Using the qs/qsc ratio from total load 1 yr sediment in a 10 yr
% flood then applying that qsqsc to the total load for the other floods and
% then using that sediment discharge

for ii=1:length(Names) %which site
    
    Width = W(ii); %channel width upstream of the waterfall (m) [measure in Google Earth]
    
    A = da(ii); %draingae area of waterfall (m2
    Slope = S_tot(ii); %slope of channel average
    S_poold = S_brink(ii); % Slope leading into the lip
    E_d=e_ba(ii);
    A_gage=A_gaget(ii);

    D_50=D50(ii);
    D_85=D85(ii);
    Slope_np=S_nowf(ii);

    
    qsqsc=qsqsctl_10(ii);
    %Discharges for 2, 5, 10, 50 and 100 yr floods
    h_brd =depth(ii); %pool depth to bedrock (m) [we don't know how deep the pools actually are, so better
                %just to keep this set at a constant value]
    H_dropd =dh(ii); %waterfall drop height (m) [as with above, for a first order appproximation 
                    %we can just hold waterfall drop height constant between
                    %reaches]
    r_poold = radius(ii); %Pool radius (m)
    
    W_brink=W(ii); %Width at wf lip
    equ=width_equ(ii, 2:3);

    for n=1:length(Qg_k) %which discharge
        %qsctot=qsc_all(ii, n);
        Qg=Q_gage(ii, n);
        If_s=If_sed(n);

        [Etl_nowf( ii, n),   Htl_nowf( ii, n), Wtl_nowf(ii, n), Qflood_nowftl(ii, n), qsc_tl(ii, n), qs_tl(ii,n)]=...
                only_tl_sites(Slope,  Width,  Qg, If,  A, A_gage, E_d, D_50, D_85, If_s, equ, qsqsc);

        [E_vert_inst_arrayd(ii, n),  E_lat_inst_arrayd( ii, n), ...
        Etl_distd( ii, n), Htl( ii, n), Hbrink(  ii, n),  h_TW(  ii, n), If, Qflood(ii, n),  r_jet( ii, n),  Width_tl( ii, n),  Widthbrink( ii, n),  qs_all(ii, n)]=...
        all_erosion_sites(Slope_np, S_poold, Width, W_brink, Qg, If, h_brd, H_dropd, r_poold, A, A_gage, E_d, D_50, D_85, If_s, equ, qs_tl(ii,n));

        %Qsc_pool(ii, n) = Qsc_pool_public(r_poold, Qflood(ii,n), h_brd, H_dropd, Widthbrink(ii,n), D_50, Hbrink(ii,n), ...
        %100, 100);
                             
        %E_vert_inst_arrayd(iii,  ii, n)
        %Qflood(iii, ii, n)
    end
%         E_tl_avg_perdischarge=sum(Drat_d.*Etl_distd, 2 );
%         E_vert_perdischarge=sum(Drat_d.*E_vert_inst_arrayd, 2 );
%         E_tl_brink_perdischarge=sum(Drat_d.*Etl_brink_d, 2);
%         E_lat_perdischarge=sum(Drat_d.*E_lat_inst_arrayd, 2);
%         E_all_perdischarge=[transpose(E_vert_perdischarge)
%         transpose(E_lat_perdischarge)
%         transpose(E_tl_brink_perdischarge)
%         transpose(E_tl_avg_perdischarge)];
end

%     if isnan(Wlip(1,1))
% 
%         for n=1:length(Qg_k) %which discharge
%             Width = W(ii); %channel width upstream of the waterfall (m) [measure in Google Earth]
%             
%             A = da(ii); %draingae area of waterfall (m2
%             Slope = S(ii); %slope of channel average
%             E_d=E_dtot(ii);
%             A_gage=A_gaget=(ii);
%             Qg=Q_gage(ii, n);
%             If_s=If_sed(n);
%             [Etl_distd(1, ii, n),   H(1,  ii, n), Qflood(1, ii, n)]=...
%                 only_tl_formodeling(Slope,  Width,  Qg, If,  A, A_gage, E_d, D, If_s);
%             %Qflood(1, 1, n)
%         end
% 
%     end

%end
Utl=Qflood./(Width_tl.*Htl);
Fr_tl=Utl./((9.8.*Htl).^0.5);

%% integrating across discharges
E_vert_inst_array_if=E_vert_inst_arrayd;
E_lat_inst_array_if=E_lat_inst_arrayd;
Etl_if=Etl_distd;
Etl_brink_if=Etl_brink_d;
Etl_ifnowf=Etl_nowf;

for ijk=1:length(E_vert_inst_arrayd(1, :))
    E_vert_inst_array_if( :, ijk)=E_vert_inst_arrayd( :, ijk).*Ifyearly(ijk)./If_sed(ijk);
    E_lat_inst_array_if( :, ijk)=E_lat_inst_arrayd( :, ijk).*Ifyearly(ijk)./If_sed(ijk);
    Etl_if( :, ijk)=Etl_distd( :, ijk).*Ifyearly(ijk)./If_sed(ijk);
    Etl_ifnowf(:, ijk)=Etl_nowf( :, ijk).*Ifyearly(ijk)./If_sed(ijk);
end

%E_vert_tot_perwf=sum(E_vert_inst_array_if, 2, 'omitnan');
%E_lat_tot_perwf=sum(E_lat_inst_array_if, 2, 'omitnan');
%Etl_brink_tot_perwf=sum(Etl_brink_if, 2, 'omitnan');
%Etl_tot_perwf=sum(Etl_if, 2, 'omitnan');
%Etl_nowf_perwf=sum(Etl_ifnowf, 2, 'omitnan');
Qcat = categorical(Qflood);
Dcat=categorical(D_50);
%Qcat = reordercats(Qcat,Q_gage);
Ecat={'WF Vert Erosion', 'WF Lateral Erosion',  'Total Load Avg Channel'};
Ecat=categorical(Ecat);
%% FIgures for each type of erosion (length(da), length(Q_gage), length(Names));
% Qcat = categorical(Qflood);
% Dcat=categorical(D);
% %Qcat = reordercats(Qcat,Q_gage);
% Ecat={'WF Vert Erosion', 'WF Lateral Erosion', 'Total Load @Brink', 'Total Load Avg Channel'};
% Ecat=categorical(Ecat);
% 
% % WF Vertical
% figure
% bar( E_vert_inst_arrayd(:, :, 1), 'grouped')
% title('WF Vert Erosion by and flood size')
% xlabel(' Waterfall #')
% ylabel('Erosion without intermittancy')
% legend('Yearly', '2 Year', '5 Year', '10 Year', '50 Year', '100 Year', 'Location','northwest')
% 
% figure
% bar( E_lat_inst_arrayd(:, :, 1), 'grouped')
% title('WF Lat Erosion by grain size and flood size')
% xlabel('Waterfall #')
% ylabel('Erosion without intermittancy')
% legend('Yearly', '2 Year', '5 Year', '10 Year', '50 Year', '100 Year', 'Location','northwest')
% 
% figure
% bar(Dcat, Etl_brink_d, 'grouped')
% title('Total Load Brink Erosion by grain size and flood size')
% xlabel('Grain Diameter (m)')
% ylabel('Erosion without intermittancy')
% legend('Yearly', '2 Year', '5 Year', '10 Year', '50 Year', '100 Year', 'Location','northwest')
% 
% 
% 
% 
% %Total Load Stream Average
% figure
% bar( Etl_distd(:, :, 1), 'grouped')
% title('Total Load Erosion by flood size')
% xlabel(' Grain Diameter (m)')
% ylabel('Erosion without intermittancy')
% legend('Yearly', '2 Year', '5 Year', '10 Year', '50 Year', '100 Year', 'Location','northwest')
% % %% Figures per Discharge across grain sizes
% % 
% % E_tl_avg_perdischarge=sum(Drat_d.*Etl_distd, 2 );
% % E_vert_perdischarge=sum(Drat_d.*E_vert_inst_arrayd, 2 );
% % E_tl_brink_perdischarge=sum(Drat_d.*Etl_brink_d, 2);
% % E_lat_perdischarge=sum(Drat_d.*E_lat_inst_arrayd, 2);
% % E_all_perdischarge=[transpose(E_vert_perdischarge)
% %     transpose(E_lat_perdischarge)
% %     transpose(E_tl_brink_perdischarge)
% %     transpose(E_tl_avg_perdischarge)];
% % figure
% % bar(Qcat, E_all_perdischarge)
% % legend('WF Vert Erosion', 'WF Lateral Erosion', 'Total Load @Brink', 'Total Load Avg Channel', 'location', 'northoutside')
% % xlabel('Flood DIscharge (m^3/sec')
% % ylabel('Erosion (no intermitancy)')
% % title('Erosion per flood type, across grain sizes')
% % 
% % %% Applying Intermittancy to flood levels
% % E_vert_tot=sum(E_vert_perdischarge.*Ifyearly, 'all');
% % E_lat_tot=sum(E_lat_perdischarge.*Ifyearly, 'all');
% % E_tl_avg_tot=sum(E_tl_avg_perdischarge.*Ifyearly, 'all');
% % E_tl_brink_tot=sum(E_tl_brink_perdischarge.*Ifyearly, 'all');
% % E_all=[E_vert_tot ...
% %     E_lat_tot ...
% %     E_tl_brink_tot ...
% %     E_tl_avg_tot];
% % figure
% % bar(Ecat, E_all)
% % ylabel('Erosion (mm/yr)')
% % title('MFK 11')
%% Making Figures per wf
% 
% for i=1:length(Names)
%     barstuff=[E_vert_tot_perwf(i), E_lat_tot_perwf( i),  Etl_tot_perwf(i)];
%     figure
%     
%     bar(barstuff, 'grouped')
%     %xticklabels(1:7)
%     legend(Ecat)
%     title(Names(i))
%     hold off
% end
%hold on
%bar(log(D), E_lat_inst_arrayd, 'grouped')
%hold off
%% Combining erosion to be per site

Etot=NaN(length(Names), 3, length(E_vert_inst_array_if(1,:)) );%Column 1 is waterfall vert, Column 2 is watervall Lat, column 3 is total load
Enowf=NaN(length(Names), 1);
wfvert_weighted=NaN(length(W(:,1)), length(Names) ); %these arrays contain erosion rates for each site multiplied by the distance along which they apply
wflat_weighted=NaN(length(W(:,1)), length(Names) );
tl_weighted=NaN(length(W(:,1)), length(Names) );
length_tl=NaN(length(Names), 1);
for i=1:length(Names)
    for j=1:length(E_vert_inst_array_if(1,:))
        Etot(i, 1, j)=(E_vert_inst_array_if(i, j)*wf_freq(i)*radius(i)*2)./length_ch(i);
        Etot(i, 2, j)=(E_lat_inst_array_if(i, j)*wf_freq(i)*depth(i)*2)./length_ch(i);
        Etot(i, 3, j)=(Etl_if(i, j)*length_np(i))./length_ch(i);
        Enowf(i, j)=Etl_ifnowf(i, j)*length_ch(i);

%     length_tl(i)=Reach_len(i)-(sum(poolrad(:,i), 'omitnan')*2);
%     for j=1:length(Wlip(:,1))
%         wfvert_weighted(j,i)=E_vert_tot_perwf(j,i)*poolrad(j,i)*2;
%         wflat_weighted(j,i)=E_lat_tot_perwf(j,i)*depth(j,i)*2;
%         tl_weighted(j,i)=Etl_tot_perwf(j,i)*(length_tl(i));
%     end
%     Etot(i, 1)=sum(wfvert_weighted(:,i), 'all', 'omitnan');
%     Etot(i, 2)=sum(wflat_weighted(:,i), 'all', 'omitnan');
%     Etot(i,3)=tl_weighted(1,i);
    end
end

Ecat1=['WF Vert', 'WF Lat', "Total Load"];
Ecat1=categorical(Ecat1);
types={'2 yr flood'
    '5 yr flood'
    '10 yr flood'
    '15 yr flood'
    '25 yr flood'};
% for j=1:length(E_vert_inst_array_if(1,:))
%     figure
%     bar(Etot(:,:,j)./length_ch(j), 'grouped')
%     xticklabels(Names)
%     
%     legend(Ecat1)
%     title(types(j))
%     hold off
% end



%% combining for sitewide predictions
Esite=NaN(length(Names), 1);
for i=1:length(Names)
    for j=1:length(E_vert_inst_array_if(1,:))
        Esite(i, j)=sum(Etot(i, :, j), 'all', 'omitnan');
        Esite_wfonly(i, j)=sum(Etot(i, 1:2, j), 'all', 'omitnan');
    end
end

% for j=1:length(E_vert_inst_array_if(1,:))    
%     figure
%     bar(Esite(:,j))
%     xticklabels(Names)
%     ylabel('total erosion (mm*m/yr)')
%     title(types(j))
% end

%% esite at a point
Esite_point=NaN(length(Names), length(E_vert_inst_array_if(1,:)));
for i=1:length(Names)
    for j=1:length(E_vert_inst_array_if(1,:))
        Esite_point(i, j)=Esite(i, j);
        Esite_point_wfonly(i, j)=Esite_wfonly(i, j);
    end
end

% for j=1:length(E_vert_inst_array_if(1,:))
%     figure
%     bar(Esite_point(:,j))
%     xticklabels(Names)
%     ylabel('total erosion at a point')
%     title(types(j))
% end

%% current erosion estimates
% e_data=[195
% 80.2
% 295
% 451
% 181
% 51.9
% 44.8
% 92.1
% 172
% 99.9
% ];
Names={'MFK21-06'
'MFK21-04'
'MFK21-11'
'MFK22-02'
'MFK21-09'
'DNK21-06'
'DNK22-01'
'DNK22-02'
'DNK21-0304'
'DNK22-03'};

e_data=[172
    72.6
    258
    351
    159
    49.6
    39.2
    85
    154
    90.8]; %updated rates, taking into account topographic shielding, water shielding, and slope shielding

error=[10.200142
4.30003
15.200389
36.9
14.1
2.95000513
3.198
6.93
9.080124
8.96];


e_data=e_data./1000;
for j=1:length(E_vert_inst_array_if(1,:))
    figure
    scatter(e_data, Esite_point(:,j))
    %hold on
    %scatter(e_data, Esite_point_wfonly(:,j))
    xlabel(' measured erosion rate (mm/yr)')
    ylabel('predicted erosion rate (mm/yr')
    legend('full model', 'wf only')
    title(types(j))
end

%% pure waterfall and total load comparison
% 
% for j=1:length(E_vert_inst_array_if(1,:))
%     figure
%     scatter(e_data, E_vert_inst_array_if(:, j))
%     xlabel(' measured erosion rate (mm/yr)')
%     ylabel('predicted erosion rate (mm/yr')
%     %text(.25,2e4,'waterfall vertical erosion')
%     title(types(j), 'waterfall vertical erosion at a point')
% 
%     figure
%     scatter(e_data, (E_vert_inst_array_if(:, j)+E_lat_inst_array_if(:, j)))
%     xlabel(' measured erosion rate (mm/yr)')
%     ylabel('predicted erosion rate (mm/yr')
%     title(types(j), 'waterfall vertical and lateral erosion')
% 
%     figure
%     scatter(e_data, Etl_ifnowf(:,j))
%     xlabel(' measured erosion rate (mm/yr)')
%     ylabel('predicted erosion rate (mm/yr')
%     title(types(j), 'pure total load')
% end
%% looking at qs/qsc
for i=1:length(e_data)
    for j=1:length(Q_gage(1,:))
        %qsqsc_pool(i,j)=qs_all(i,j)/Qsc_pool(i,j);
        qsqsc_tl(i,j)=qs_tl(i,j)/qsc_tl(i,j); %qsc_tl(ii, n), qs_tl(ii,n)
    end
end


%% normalized wf erosion rate v waterfall frequency
e_data_norm=NaN(length(e_data));
e_data_norm(1:5)=e_data(1:5)./92.24;
e_data_norm(6:10)=e_data(6:10)./32.17;

figure
scatter(wf_freq, e_data_norm)
ylabel('Normalized erosion rates')
xlabel('Wf occurence')


%%
figure

h = tiledlayout(2,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
%box on
%fitobj=fitlm(e_data(1:10), Esite_point(1:10, 3), "Exclude",8)
%plot([0, 2], [fitobj.Coefficients{1, 'Estimate'}, 2*fitobj.Coefficients{2, 'Estimate'}+fitobj.Coefficients{1, 'Estimate'}])

%errorbar(Esite_point(1:5, 3), e_data(1:5), error(1:5)./1000, 'k-','LineStyle','none')
%hold on

errorbar( e_data(1:5), Esite_point(1:5, 3),error(1:5)./1000,'horizontal','ks', 'markerFaceColor', [0.4118    0.7804    0.9412],'LineStyle','none', 'MarkerSize', 6, 'Markeredgecolor', 'k')
hold on
errorbar(e_data(6:10),Esite_point(6:10, 3), error(6:10)./1000,'horizontal','ko', 'markerFaceColor', [0.8118    0.4118    0.7059], 'LineStyle','none', 'MarkerSize', 6, 'Markeredgecolor', 'k')
plot([0, 0.9], [0,0.9],  'k--')
plot([0, 0.9], [0,3.2],  '--')
hold on

xlabel('Channel averaged erosion rate (mm/yr)')
ylabel('Combined model prediction (mm/yr)')
legend('MFK', 'Dinkey Creek', '1:1')
%ylabel('Combined model prediction (mm/yr)')
xlim([0,0.4])
ylim([0, .4])
%box off

nexttile

errorbar(e_data(1:10),E_vert_inst_array_if(:, 2),error(1:10)./1000, 'horizontal','k^', 'linestyle','none', 'Markeredgecolor', 'k', 'MarkerFaceColor',[0.4941    0.1843    0.5569], 'MarkerSize', 6)
hold on
errorbar(e_data(1:10),E_lat_inst_array_if(:,  2),error(1:10)./1000,'horizontal',  'kd', 'linestyle','none',  'Markeredgecolor', 'k', 'MarkerFaceColor', [1.0000   ,0.4118   ,0.1608],'MarkerSize', 6)
errorbar( e_data(1:10),Etl_if(:, 2),error(1:10)./1000,'horizontal', 'kp', 'linestyle','none','MarkerSize', 6)
plot([0, 3], [0,3], 'k--')
legend('E_{wfvert}', 'E_{wflat}', 'E_{tl}', '1:1')
%sum(Etot(i, :, j), 'all', 'omitnan');
xlabel('Channel averaged erosion rate (mm/yr)')
ylabel('E_{tl}, E_{wfvert}, E_{wflat} (mm/yr)')
xlim([0, .4])
ylim([0,2.2])
%% good fig but just combined model
figure

%h = tiledlayout(2,1, 'TileSpacing', 'tight', 'Padding', 'compact');
%nexttile
%box on
%fitobj=fitlm(e_data(1:10), Esite_point(1:10, 3), "Exclude",8)
%plot([0, 2], [fitobj.Coefficients{1, 'Estimate'}, 2*fitobj.Coefficients{2, 'Estimate'}+fitobj.Coefficients{1, 'Estimate'}])

%errorbar(Esite_point(1:5, 3), e_data(1:5), error(1:5)./1000, 'k-','LineStyle','none')
mask=[true, true, true, true, true, true, true, false, true, true]
[rho, p]=corr(e_data(mask), Esite_point(mask, 3))
%hold on

errorbar( e_data(1:5), Esite_point_wfonly(1:5, 3),error(1:5)./1000,'horizontal','ks', 'markerFaceColor', [0.4118    0.7804    0.9412],'LineStyle','none', 'MarkerSize', 6, 'Markeredgecolor', 'k')
hold on
errorbar(e_data(6:10),Esite_point(6:10, 3), error(6:10)./1000,'horizontal','ko', 'markerFaceColor', [0.8118    0.4118    0.7059], 'LineStyle','none', 'MarkerSize', 6, 'Markeredgecolor', 'k')
plot([0, 1], [0,1],  'k--')
plot([0, 1], [0,5],  ':', Color=[0.4 0.4, 0.4])
plot([0, 1], [0,(1/5)],  '-.', Color=[0.4, 0.4, 0.4])
hold on

xlabel('Channel averaged erosion rate (mm/yr)')
ylabel('Combined model prediction (mm/yr)')
legend('MFK', 'Dinkey Creek', '1:1', '1:5', '5:1')
%ylabel('Combined model prediction (mm/yr)')
xlim([0,0.8])
ylim([0, .8])
%box off

% nexttile
% 
% errorbar(e_data(1:10),Etot(:, 1, 3),error(1:10)./1000, 'horizontal','k^', 'linestyle','none', 'Markeredgecolor', 'k', 'MarkerFaceColor',[0.4941    0.1843    0.5569], 'MarkerSize', 6)
% hold on
% errorbar(e_data(1:10),Etot(:, 2, 3),error(1:10)./1000,'horizontal',  'kd', 'linestyle','none',  'Markeredgecolor', 'k', 'MarkerFaceColor', [1.0000   ,0.4118   ,0.1608],'MarkerSize', 6)
% errorbar( e_data(1:10),Etot(:, 3, 3),error(1:10)./1000,'horizontal', 'kp', 'linestyle','none','MarkerSize', 6)
% %plot([0, 3], [0,3], 'k--')
% legend('E_{wfvert}', 'E_{wflat}', 'E_{tl}', '1:1')
% %sum(Etot(i, :, j), 'all', 'omitnan');
% xlabel('Channel averaged erosion rate (mm/yr)')
% ylabel('E_{tl}, E_{wfvert}, E_{wflat} (mm/yr)')
% xlim([0, .4])
% ylim([0,0.5])
%% fig with combined, tl, and wf
Esitegood=Esite_point(1:10, 3)<5;
[rho, pval]=corr(e_data(Esitegood), Esite_point(Esitegood, 3))
[rho, pval]=corr(e_data(Esitegood), E_vert_inst_array_if(Esitegood, 3)+E_lat_inst_array_if(Esitegood, 3))
[rho, pval]=corr(e_data(Esitegood), Etl_ifnowf(Esitegood, 3))
edata_fit=e_data(Esitegood);
Esite_fit=Esite_point(Esitegood, 3);
figure
fit1=fitlm(edata_fit,  Etl_ifnowf(Esitegood, 3))
h = tiledlayout(3, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
%nexttile
%box on
%fitobj=fitlm(e_data(1:10), Esite_point(1:10, 3), "Exclude",8)
%plot([0, 2], [fitobj.Coefficients{1, 'Estimate'}, 2*fitobj.Coefficients{2, 'Estimate'}+fitobj.Coefficients{1, 'Estimate'}])

%errorbar(Esite_point(1:5, 3), e_data(1:5), error(1:5)./1000, 'k-','LineStyle','none')
%hold on
nexttile
errorbar( e_data(1:5), Esite_point(1:5, 3),error(1:5)./1000,'horizontal','ks', 'markerFaceColor', [0.4118    0.7804    0.9412],'LineStyle','none', 'MarkerSize', 6, 'Markeredgecolor', 'k')
hold on
errorbar(e_data(6:10),Esite_point(6:10, 3), error(6:10)./1000,'horizontal','ko', 'markerFaceColor', [0.8118    0.4118    0.7059], 'LineStyle','none', 'MarkerSize', 6, 'Markeredgecolor', 'k')
plot([0, 1], [0,1],  'k--')
plot([0, 1], [0,5],  ':', Color=[0.4 0.4, 0.4])
plot([0, 1], [0,(1/5)],  '-.', Color=[0.4, 0.4, 0.4])
hold on

xlabel('Channel averaged erosion rate (mm/yr)')
ylabel('Combined model prediction (mm/yr)')
%legend('MFK', 'Dinkey Creek', '1:1', '1:5', '5:1')
%ylabel('Combined model prediction (mm/yr)')
xlim([0,0.4])
ylim([0, 7])


nexttile
errorbar( e_data(1:5), (E_vert_inst_array_if(1:5, 3)+E_lat_inst_array_if(1:5, 3)),error(1:5)./1000,'horizontal','ks', 'markerFaceColor', [0.4118    0.7804    0.9412],'LineStyle','none', 'MarkerSize', 6, 'Markeredgecolor', 'k')
hold on
errorbar(e_data(6:10),(E_vert_inst_array_if(6:10, 3)+E_lat_inst_array_if(6:10, 3)), error(6:10)./1000,'horizontal','ko', 'markerFaceColor', [0.8118    0.4118    0.7059], 'LineStyle','none', 'MarkerSize', 6, 'Markeredgecolor', 'k')
plot([0, 1], [0,1],  'k--')
plot([0, 1], [0,5],  ':', Color=[0.4 0.4, 0.4])
plot([0, 1], [0,(1/5)],  '-.', Color=[0.4, 0.4, 0.4])
hold on

xlabel('Channel averaged erosion rate (mm/yr)')
ylabel('Plunge-pool erosion perdiction (mm/yr)')
legend('MFK', 'Dinkey Creek', '1:1', '1:5', '5:1')
%ylabel('Combined model prediction (mm/yr)')
xlim([0,0.4])
%ylim([0, 10])

nexttile
errorbar( e_data(1:5), Etl_ifnowf(1:5,3),error(1:5)./1000,'horizontal','ks', 'markerFaceColor', [0.4118    0.7804    0.9412],'LineStyle','none', 'MarkerSize', 6, 'Markeredgecolor', 'k')
hold on
errorbar(e_data(6:10),Etl_ifnowf(6:10, 3), error(6:10)./1000,'horizontal','ko', 'markerFaceColor', [0.8118    0.4118    0.7059], 'LineStyle','none', 'MarkerSize', 6, 'Markeredgecolor', 'k')
plot([0, 1], [0,1],  'k--')
plot([0, 1], [0,5],  ':', Color=[0.4 0.4, 0.4])
plot([0, 1], [0,(1/5)],  '-.', Color=[0.4, 0.4, 0.4])
hold on

xlabel('Channel averaged erosion rate (mm/yr)')
ylabel('Total load prediction (mm/yr)')
%legend('MFK', 'Dinkey Creek', '1:1', '1:5', '5:1')
%ylabel('Combined model prediction (mm/yr)')
xlim([0,0.4])
%ylim([0, 0.4])
%box off

% nexttile
% 
% errorbar(e_data(1:10),Etot(:, 1, 3),error(1:10)./1000, 'horizontal','k^', 'linestyle','none', 'Markeredgecolor', 'k', 'MarkerFaceColor',[0.4941    0.1843    0.5569], 'MarkerSize', 6)
% hold on
% errorbar(e_data(1:10),Etot(:, 2, 3),error(1:10)./1000,'horizontal',  'kd', 'linestyle','none',  'Markeredgecolor', 'k', 'MarkerFaceColor', [1.0000   ,0.4118   ,0.1608],'MarkerSize', 6)
% errorbar( e_data(1:10),Etot(:, 3, 3),error(1:10)./1000,'horizontal', 'kp', 'linestyle','none','MarkerSize', 6)
% %plot([0, 3], [0,3], 'k--')
% legend('E_{wfvert}', 'E_{wflat}', 'E_{tl}', '1:1')
% %sum(Etot(i, :, j), 'all', 'omitnan');
% xlabel('Channel averaged erosion rate (mm/yr)')
% ylabel('E_{tl}, E_{wfvert}, E_{wflat} (mm/yr)')
% xlim([0, .4])
% ylim([0,0.5])