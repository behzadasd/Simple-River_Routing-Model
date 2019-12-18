%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Flood Routing Grid-cells Initiation  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Behzad Asadieh, Ph.D.           %%%
%%% University of Pennsylvania      %%%
%%% basadieh@sas.upenn.edu          %%%
%%% github.com/behzadasd            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
clear;
clc;

%[HydroSTN, HydroSTN_txt, ~]=xlsread('HydroSTN-excel.xlsx', 'HydroSTN');
load HydroSTN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Creating 6-minute resolution grids (0.1 Degree) %%%
Lat=(89.95:-0.1:-89.95)'; Lat_n=size(Lat,1);
Lon=(-179.95:0.1:179.95)'; Lon_n=size(Lon,1);

%%% Latitude Bounds %%%
Lat_bound=zeros(Lat_n,2); Lat_bound(1,1)=90; Lat_bound(1,2)=Lat(1,1)+((Lat(2,1)-Lat(1,1))/2);
for i=2:Lat_n-1
    Lat_bound(i,1)=Lat_bound(i-1,2);
    Lat_bound(i,2)=Lat(i,1)+((Lat(i+1,1)-Lat(i,1))/2);
end
Lat_bound(Lat_n,1)=Lat_bound(Lat_n-1,2); Lat_bound(Lat_n,2)=-1*Lat_bound(1,1);
%%% Longitude Bounds %%%
Lon_bound=zeros(Lon_n,2); Lon_bound(1,1)=-180; Lon_bound(1,2)=Lon(1,1)+((Lon(2,1)-Lon(1,1))/2);
for i=2:Lon_n-1
    Lon_bound(i,1)=Lon_bound(i-1,2);
    Lon_bound(i,2)=Lon(i,1)+((Lon(i+1,1)-Lon(i,1))/2);
end
Lon_bound(Lon_n,1)=Lon_bound(Lon_n-1,2); Lon_bound(Lon_n,2)=180;

%%% Hydro STN %%%
CellID=NaN(Lat_n, Lon_n);
XCoord=NaN(Lat_n, Lon_n);
YCoord=NaN(Lat_n, Lon_n);
BasinID=NaN(Lat_n, Lon_n);
ToCell=NaN(Lat_n, Lon_n);
CellArea=NaN(Lat_n, Lon_n);
CellLength=NaN(Lat_n, Lon_n);
CellXCoord=NaN(Lat_n, Lon_n);
CellYCoord=NaN(Lat_n, Lon_n);

ToCell_i=NaN(Lat_n, Lon_n); % Stores the row number of ToCell - ( Flow of the cells flowing into ocean is stored at the cell (1,1) )
ToCell_j=NaN(Lat_n, Lon_n); % Stores the column number of ToCell - ( Flow of the cells flowing into ocean is stored at the cell (1,1) )


for i=1:size(HydroSTN,1)
    
    ct=0; % The ct counter is 1 if the corresponding coordinates of HydroSTN(i,:) is found, otherwise its 0
    for ln=1:Lon_n
        if (HydroSTN(i,2) > Lon_bound(ln,1)) && (HydroSTN(i,2) < Lon_bound(ln,2))
            
            for lt=1:Lat_n
                if (HydroSTN(i,3) < Lat_bound(lt,1)) && (HydroSTN(i,3) > Lat_bound(lt,2))
                    
                    CellID(lt,ln)=HydroSTN(i,1);
                    XCoord(lt,ln)=HydroSTN(i,2);
                    YCoord(lt,ln)=HydroSTN(i,3);
                    BasinID(lt,ln)=HydroSTN(i,4);
                    ToCell(lt,ln)=HydroSTN(i,5);
                    CellArea(lt,ln)=HydroSTN(i,6);
                    CellLength(lt,ln)=HydroSTN(i,7);
                    CellXCoord(lt,ln)=HydroSTN(i,8);
                    CellYCoord(lt,ln)=HydroSTN(i,9);
                    
                    ct=1; % The corresponding coordinates of HydroSTN(i,:) is found, set ct=1 to exit the search loops
                end
                
                if ct==1
                    break
                end
            end
            
        end
        
        if ct==1
            break
        end
        
    end
    
    if ~(mod(i,10000)) % mod(i,10000) returns remainder of i/10000 - this will work every 10000 cell
        disp(['CellID ', num2str(i), ' processed']) 
    end

end

%%% Calculating row and column number of ToCell
for lt=1:Lat_n
    for ln=1:Lon_n
        
        if ~isnan(ToCell(lt,ln))
            
            for ci=lt-1:lt+1
                for cj=ln-1:ln+1
                    if CellID(ci,cj)==ToCell(lt,ln)
                        ToCell_i(lt,ln)=ci;
                        ToCell_j(lt,ln)=cj;
                        
                    end
                end
            end
            
            if ToCell(lt,ln)==-1 % That means the cell flows into ocean
                ToCell_i(lt,ln)=1; % Flow of the cells flowing into ocean is stored at the cell (1,1)
                ToCell_j(lt,ln)=1;
            end
            
        end
        
    end
end





Variable=CellID;

figure
load coast
geoshow(flipud(lat), flipud(long),'DisplayType', 'polygon', 'FaceColor', [0.94 0.97 1])
%geoshow('landareas.shp', 'FaceColor', [0.97 0.97 0.97]);  %%% Enable this to Change the LandArea color
%states = shaperead('usastatehi', 'UseGeoCoords', true);    %%% Enable this to add USA state boundaries
%geoshow(states,'DefaultFaceColor', 'white', 'DefaultEdgeColor', 'blue', 'FaceColor', 'white');
xlim([-180 180]); ylim([-90 90]);
hold on
h=imagesc(Lon', Lat, Variable);
set(h,'alphadata',~isnan(Variable)) % Sets NaN values no color (colors them white)
%set(gca, 'CLim', var_limit);
%colormap jet
set(gca,'YDir','normal') % Prevents fliping latitude directions because of the upper arrays being greater than lower arrays
%colorbar('location','Eastoutside')
set(gca,'FontSize',34, 'FontName', 'MS Sens Serif') % Axis Numbers and rages Font




toc;
