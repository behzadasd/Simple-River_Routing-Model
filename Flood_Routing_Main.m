%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Flood Routing Model     %%%
%%%  ISI-MIP - WBM RunOff data  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Behzad Asadieh, Ph.D.           %%%
%%% University of Pennsylvania      %%%
%%% basadieh@sas.upenn.edu          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
clear;
clc;

load HydroSTN_Grid

CellArea_0=CellArea;
CellArea_0(isnan(CellArea_0))=0;

Models_Names = {'GFDL-ESM2M', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC-ESM-CHEM','NorESM1-MR'};
dir=[pwd '\']; % Current Directory Path
%dir_data_in1=[dir 'NetCDF Raw Data - RUNOFF\']; % Directory to raed raw data from
dir_data_in1='C:\Users\Behzad\Desktop\ISI-MIP RunOff and Discharge Data\Data Initialization - RunOff (Daily)\NetCDF Raw Data - RUNOFF\';
dir_mat_out=[pwd '\Routed Data - Global\'];
dir_mat_out_model=[pwd '\Matlab Manipulated Data - Global\'];

Input_File_Names = ls( fullfile(dir_data_in1,'*.nc') ); % Lists the .csv files names in the specified directory

[~,file_name,~]=fileparts(Input_File_Names(1,:)); % Name of NetCDF file to be loaded (without the extention)
file_name_dir_in=[dir_data_in1 file_name '.nc']; % Directory and Name of .nc file to be loaded
%ncdisp(file_name_dir_in)
GCM_Lat=ncread(file_name_dir_in , 'lat');
GCM_Lon=ncread(file_name_dir_in , 'lon');
GCM_Lat_n=size(GCM_Lat,1);
GCM_Lon_n=size(GCM_Lon,1);

%%% Latitude Bounds %%%
GCM_Lat_bound=zeros(GCM_Lat_n,2); GCM_Lat_bound(1,1)=90; GCM_Lat_bound(1,2)=GCM_Lat(1,1)+((GCM_Lat(2,1)-GCM_Lat(1,1))/2);
for i=2:GCM_Lat_n-1
    GCM_Lat_bound(i,1)=GCM_Lat_bound(i-1,2);
    GCM_Lat_bound(i,2)=GCM_Lat(i,1)+((GCM_Lat(i+1,1)-GCM_Lat(i,1))/2);
end
GCM_Lat_bound(GCM_Lat_n,1)=GCM_Lat_bound(GCM_Lat_n-1,2); GCM_Lat_bound(GCM_Lat_n,2)=-1*GCM_Lat_bound(1,1);
%%% Longitude Bounds %%%
GCM_Lon_bound=zeros(GCM_Lon_n,2); GCM_Lon_bound(1,1)=-180; GCM_Lon_bound(1,2)=GCM_Lon(1,1)+((GCM_Lon(2,1)-GCM_Lon(1,1))/2);
for i=2:GCM_Lon_n-1
    GCM_Lon_bound(i,1)=GCM_Lon_bound(i-1,2);
    GCM_Lon_bound(i,2)=GCM_Lon(i,1)+((GCM_Lon(i+1,1)-GCM_Lon(i,1))/2);
end
GCM_Lon_bound(GCM_Lon_n,1)=GCM_Lon_bound(GCM_Lon_n-1,2); GCM_Lon_bound(GCM_Lon_n,2)=180;

%%% Flood Routing %%%
multiplier_i=ceil((1:GCM_Lat_n*5)/5); % will be used to multiply the run-off matrix size by 5
multiplier_j=ceil((1:GCM_Lon_n*5)/5);

yrs_n=35; % Number of the Years that data is available for
start_date='1971';
end_date='2005';

year_first=1971; % Initial year
leap_yrs=[1972; 1976; 1980; 1984; 1988; 1992; 1996; 2000; 2004; 2008];

for M_i=2:2 % Number of ISI-MIP GCM Models involved
    
    Model_Name=Models_Names{M_i};
    disch_0=zeros(Lat_n, Lon_n); disch_0(isnan(CellID))=NaN; % disch_0 is the initial water at the gridcells before the routing starts
    
    ISIMIP_WBM_HydroSTN_Discharge_1971_2005_Mean=NaN(Lat_n, Lon_n, yrs_n);
    ISIMIP_WBM_HydroSTN_Discharge_1971_2005_Max=NaN(Lat_n, Lon_n, yrs_n);
    ISIMIP_WBM_HydroSTN_Discharge_1971_2005_Min=NaN(Lat_n, Lon_n, yrs_n);
    
    for F_i=1:4 % Number of .nc files per model to be processed
        
        F = (M_i-1)*4+F_i; % count of the file to be loaded in the folder
        [~,file_name,~]=fileparts(Input_File_Names(F,:)); % Name of NetCDF file to be loaded (without the extention)
        file_name_dir_in=[dir_data_in1 file_name '.nc']; % Directory and Name of .nc file to be loaded
        
        if F_i==4
            if M_i==2
                n2=4; % HadGEM2-ES model (M_i=2) has 4 years of data at the last file - yrs_n=34
            else
                n2=5; % other models have 5 years of data at the last file - yrs_n=35
            end
        else
            n2=10;
        end
        
        j_2=0;
        for yr_i=1:n2
            
            k = (F_i-1)*10 + yr_i ; % Number of the corresponding year in total (from 1 to 35)
            year_no=year_first-1+k; % The year (form 1971 to 2005)
            
            if ismember(year_no, leap_yrs) % If the current year is a leap year...
                days_n=366;
            else
                days_n=365;
            end
            
            j_1 = j_2 + 1; % starting day
            j_2 = j_1 + days_n -1;
            
            Disch_routed=zeros(Lat_n, Lon_n, 367);
            
            RUNOFF_daily=ncread(file_name_dir_in , 'mrro', [1 1 j_1], [Inf Inf days_n], [1 1 1] );
            RUNOFF_daily=permute(RUNOFF_daily,[2,1,3]); % flipping file -90degrees [ from (lon,lat,time) to (lat,lon,time) ]
            RUNOFF_daily=RUNOFF_daily * 24 * 3600; % Converting the runoff from kg/m2/s to kg/m2/day == mm/day   (kg/m2 = mm)
            
            Disch_routed(:,:,1)=disch_0; % disch_0 is the initial water at the gridcells before the routing starts - The first layer containing the disch_0 will be eliminated at the end of routing for this time period
            for t=1:days_n
                
                RUNOFF_daily_6min=RUNOFF_daily(:,:,t);
                RUNOFF_daily_6min=RUNOFF_daily_6min(multiplier_i,:); % Rows are replicated 5 times at this point
                RUNOFF_daily_6min=RUNOFF_daily_6min(:,multiplier_j); % Columns are also replicated 5 times at this point - each 30min runoff grid is converted to 5*5 6min grids
                RUNOFF_daily_6min(isnan(RUNOFF_daily_6min))=0;
                
                Disch_routed(:,:,t+1)=RUNOFF_daily_6min .* CellArea_0 * 1000; % Unit is m3/day % (Runoff * 0.001) to convert from mm/day to m/day   % (CellArea * 1000 * 1000 ) to convert from km2 to m2
                Disch_routed(1,1,t+1)=0; % the cell Disch_routed(1,1) stores the flow into ocean at each step
                
                for lt=1:Lat_n
                    for ln=1:Lon_n
                        
                        if ~isnan(ToCell_i(lt,ln))
                            
                            tci=ToCell_i(lt,ln);
                            tcj=ToCell_j(lt,ln);
                            Disch_routed(tci,  tcj, t+1)=Disch_routed(tci,  tcj, t+1) + Disch_routed(lt, ln, t);
                            
                        end
                        
                    end
                end
                
                disp(t)
                
            end
            
% % %             Disch_routed(:,:,1)=[]; % Disch_routed(:,:,1) is actually the discharge at the gridcells at the begenning (belonging to the previous time period/year) and should be eliminated from the current time period/year
            disch_0=Disch_routed(:,:,days_n+1);
            
            ISIMIP_WBM_HydroSTN_Discharge_1971_2005_Mean(:,:,k)=nanmean(Disch_routed,3);
            ISIMIP_WBM_HydroSTN_Discharge_1971_2005_Max(:,:,k)=nanmax(Disch_routed, [],3);
            ISIMIP_WBM_HydroSTN_Discharge_1971_2005_Min(:,:,k)=nanmin(Disch_routed, [],3);
            
% % %             file_mat_dir_out=[dir_mat_out 'ISIMIP_WBM_HydroSTN_FloodRouting_Discharge_' Model_Name '_' num2str(year_no) '.mat']; % Directory and Name of .mat file to be saved after analysis
% % %             save(file_mat_dir_out, 'Disch_routed', 'Model_Name', 'year_no', 'CellID', 'XCoord', 'YCoord', 'BasinID', 'ToCell', 'CellArea', 'CellLength', ...
% % %                 'CellXCoord', 'CellYCoord', 'Lat', 'Lat_n', 'Lon', 'Lon_n', 'Lat_bound', 'Lon_bound');
            
            clear RUNOFF_daily % clearing  RUNOFF_daily from workspace for more RAM memory
            clear Disch_routed
            
            disp(['Flood Routing - Model : ' Model_Name '   year # ' num2str(k) '   Processed' ])
            
        end
        
       
    end
    
    file_mat_dir_out_model=[dir_mat_out_model 'ISIMIP_WBM_HydroSTN_FloodRouting_Discharge_' Model_Name '_1971-2005' '.mat']; % Directory and Name of .mat file to be saved after analysis
    save(file_mat_dir_out_model, 'ISIMIP_WBM_HydroSTN_Discharge_1971_2005_Mean', 'ISIMIP_WBM_HydroSTN_Discharge_1971_2005_Max', 'ISIMIP_WBM_HydroSTN_Discharge_1971_2005_Min', 'Lat', 'Lat_n', ...
        'Lon', 'Lon_n', 'Lat_bound', 'Lon_bound', 'yrs_n', 'start_date', 'end_date', 'Model_Name')
    disp(['*** *** *** Saved   ' Model_Name ' *** *** ***'])
    
end


Variable=disch_ave /(24*3600) ;
%Variable= nanmean( ISIMIP_WBM_HydroSTN_Discharge_1971_2005_Mean(:,:,2:end), 3)  /(24*3600) ;
Variable(isnan(CellID))=NaN;
Variable(Variable<1)=1;
Variable=log10( Variable );
load coast
geoshow(flipud(lat), flipud(long),'DisplayType', 'polygon', 'FaceColor', [0.94 0.97 1])
xlim([-180 180]); ylim([-90 90]);
hold on
h=imagesc(Lon', Lat, Variable);
set(h,'alphadata',~isnan(Variable)) % Sets NaN values no color (colors them white)
set(gca,'YDir','normal') % Prevents fliping latitude directions because of the upper arrays being greater than lower arrays
title({'Log of Discharge Average [log(m^3/s)] - Year 1972 ' ; 'Flood Routing of WBM/ISI-MIP RunOff (HadGEM2-ES model)'});
set(gca,'FontSize',34, 'FontName', 'MS Sens Serif') % Axis Numbers and rages Font
set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
colormap jet
colorbar('location','Eastoutside')
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window


toc;
