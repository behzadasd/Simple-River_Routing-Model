%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Flood Routing Grid-cells Initiation North America  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Behzad Asadieh, Ph.D.           %%%
%%% University of Pennsylvania      %%%
%%% basadieh@sas.upenn.edu          %%%
%%% github.com/behzadasd            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
clear;
clc;

load HydroSTN_Grid
CellID_NA=CellID (191:750 , 151:1220);
XCoord_NA=XCoord (191:750 , 151:1220);
YCoord_NA=YCoord (191:750 , 151:1220);
BasinID_NA=BasinID (191:750 , 151:1220);
ToCell_NA=ToCell (191:750 , 151:1220);
CellArea_NA=CellArea (191:750 , 151:1220);
CellLength_NA=CellLength (191:750 , 151:1220);
CellXCoord_NA=CellXCoord (191:750 , 151:1220);
CellYCoord_NA=CellYCoord (191:750 , 151:1220);
Lat_NA=Lat (191:750, :); Lat_n_NA=size(Lat_NA,1);
Lon_NA=Lon (151:1220, :); Lon_n_NA=size(Lon_NA,1);
Lat_bound_NA=Lat_bound (191:750, :);
Lon_bound_NA=Lon_bound (151:1220, :);

%%% Calculating row and column number of ToCell
ToCell_i_NA=NaN(Lat_n_NA, Lon_n_NA); % Stores the row number of ToCell - ( Flow of the cells flowing into ocean is stored at the cell (1,1) )
ToCell_j_NA=NaN(Lat_n_NA, Lon_n_NA); % Stores the column number of ToCell - ( Flow of the cells flowing into ocean is stored at the cell (1,1) )

for lt=1:Lat_n_NA
    for ln=1:Lon_n_NA
        
        if ~isnan(ToCell_NA(lt,ln))
            
            for ci=lt-1:lt+1
                for cj=ln-1:ln+1
                    if CellID_NA(ci,cj)==ToCell_NA(lt,ln)
                        ToCell_i_NA(lt,ln)=ci;
                        ToCell_j_NA(lt,ln)=cj;
                        
                    end
                end
            end
            
            if ToCell_NA(lt,ln)==-1 % That means the cell flows into ocean
                ToCell_i_NA(lt,ln)=1; % Flow of the cells flowing into ocean is stored at the cell (1,1)
                ToCell_j_NA(lt,ln)=1;
            end
            
        end
        
    end
end



