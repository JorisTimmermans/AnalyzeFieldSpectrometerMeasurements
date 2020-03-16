function [  option                      ,...
            wl,nwl                      ,...
            optipar,soilspec,soilpar    ,...
            leafgreen,leafbrown,vegpar  ,...
            ang_s, ang_o                ,...
            Satellite]     	=   input_parameters(Parameter_file)

global path_home path_version path_code 
global path_input path_output
%% Files
%paths
path_home                       =   '../LUT retrieval/';
path_version                    =   'SLC_create_LUT/';

% path_code                       =   [path_home,path_version,'slc_code'];               %path of all inputs   
path_input                      =   [path_home,'data/input/'];                         %path of all inputs   
path_output                     =   [path_home,path_version,'output/'];                %path of all inputs   


% cd(path_code)
% path(path,path_code) 

% files
leaf_file                       =   'optipar.txt';                                      % Leaf file
soil_file                       =   'soilparam.txt';                                    % Soil file

%% Options in Calculations
% Digit1                       	=   0;                                                  %If <> 0 THEN no Hapke BRDF model for the soil is applied (soilspec array contains Lambertian reflectances)
% Digit2                        =   0;                                                  %If <> 0 THEN no soil moisture effect modelling is applied
% Digit3                        =   0;                                                  %If <> 0 THEN skip leaf and canopy model 
                                                                                        % (equivalent to LAI = 0)
option                          =   int32(000);                                         %Hapke_Calculation parameter (000 =not)

%% Angles
tts                             =   30;                                                  %Sun Zenith angle (degrees)
tto                             =   0;                                                  %Observing zenith angle (degrees)
psi                             =   0;                                                  %relative azimuth angle (degrees)

%% Spectrum
resolution                      =   1;    
wl_vis                          =   350:resolution:850;
wl_nir                          =   850:resolution:2500;

wl                              =   single(unique([wl_vis wl_nir]));
nwl                             =   int32(size(wl,2));

%% Reading Parameters
path2file                       =   [Parameter_file];
[lai, lidfa, lidfb, hot, fb, diss, cv, zeta, ...
 cab_g, cw_g, cdm_g, cs_g, n_g,...
 cab_b, cw_b, cdm_b, cs_b, n_b,...
 H_SM] =   ReadParameters(path2file);

%% Reading and interpolate coefficient files for different wavelengths
% Leaf parameters => Prospect parameters
soiltype                        =	1;

optipar_raw                     =   load([path_input,'prospect_parameters/',leaf_file]);

wavl                            =   optipar_raw(:,1)*1000;
Leaf_refr                       =   interp1(wavl,optipar_raw(:,2),wl);
Kdm                             =   interp1(wavl,optipar_raw(:,3),wl);                  % dry matter
Kab                             =   interp1(wavl,optipar_raw(:,4),wl);                  % Chlorophyll
Kw                              =   interp1(wavl,optipar_raw(:,5),wl);                  % Water
Ks                              =   interp1(wavl,optipar_raw(:,6),wl);                  % Senescent
H20_refr                        =   interp1(wavl,optipar_raw(:,7),wl);                  % refraction
H20_abs                         =   interp1(wavl,optipar_raw(:,8),wl);                  % absorption
clear optipar_raw

% Soil parameters => Hapke parameters
soilpar_raw                     =   load([path_input,'soil_spectrum/',soil_file]);        
H_b                             =   soilpar_raw(2,soiltype);                            %Hapke B parameter
H_c                             =   soilpar_raw(3,soiltype);                            %Hapke c parameter
H_B0                            =   soilpar_raw(4,soiltype);                            %Hapke B0 parameter
H_h                             =   soilpar_raw(5,soiltype);                            %Hotspot parameter
soil_refr                       =   interp1(wavl,soilpar_raw(6:end,soiltype),wl);       %Soil nr
clear soilpar_raw

%% Sensor Specifics
global sensor
global bandwidth
global minwl maxwl 
Satellite.sensor                        =   sensor;
switch lower(sensor)
    case 'modis'
        path2files                      =   [path_input,'SensorsSpecifics/'];
        files                           =   dir([path2files,'*.dat']);
        Nr_bands                        =   length(files);
        for i=1:Nr_bands
            [wl_raw,S_raw]              =   textread([path2files,files(i).name],'%f%*f%f%*f','headerlines',3);

            Satellite.Raw(i).wl         =   wl_raw;
            Satellite.Raw(i).S          =   S_raw;
            S(:,i)                      =   interp1(wl_raw,S_raw,wl);
        end
        S(isnan(S))                     =   0;

        Satellite.Nr_bands              =   Nr_bands;
        Satellite.Wavelength            =   wl;
        Satellite.Bandwidth             =   resolution;
        Satellite.Sensitivity           =   S;        
    case 'rapideye'
        rapideye_min                    =   [ 440, 520, 630, 690, 760];
        rapideye_max                    =   [ 510, 590, 685, 730, 850];
        Satellite.Wavelength_c          =   mean([rapideye_min;rapideye_max]);
        Satellite.Nr_bands              =   length(rapideye_min);
        Satellite.Wavelength            =   wl;
        Satellite.Bandwidth             =   resolution;
        
        Satellite.Sensitivity           =   zeros(length(wl),Satellite.Nr_bands);
        for iband = 1:Satellite.Nr_bands
            band_wl                     =   rapideye_min(iband):rapideye_max(iband);
            sensitivity_wl              =   ones(size(band_wl));
            Satellite.Sensitivity(:,iband)=   interp1(band_wl, sensitivity_wl,Satellite.Wavelength,'nearest');
        end
        Satellite.Sensitivity(isnan(Satellite.Sensitivity))=0;
    
    case 'spot'
        spot_min                        =   [ 438, 615, 772, 1564];
        spot_max                        =   [ 486, 696, 914, 1635];
        Satellite.Wavelength_c          =   mean([spot_min;spot_max]);
        Satellite.Nr_bands              =   length(spot_min);
        Satellite.Wavelength            =   wl;
        Satellite.Bandwidth             =   resolution;
        
        Satellite.Sensitivity           =   zeros(length(wl),Satellite.Nr_bands);
        for iband = 1:Satellite.Nr_bands
            band_wl                     =   spot_min(iband):spot_max(iband);
            sensitivity_wl              =   ones(size(band_wl));
            Satellite.Sensitivity(:,iband)=   interp1(band_wl, sensitivity_wl,Satellite.Wavelength,'nearest');
        end
        Satellite.Sensitivity(isnan(Satellite.Sensitivity))=0;
        
    case {'hyperspectral', 'rs3500'}
        Satellite.Bandwidth             =   bandwidth;
        hypspec_min                     =   (minwl:Satellite.Bandwidth:maxwl) - Satellite.Bandwidth/2;
        hypspec_max                     =   (minwl:Satellite.Bandwidth:maxwl) + Satellite.Bandwidth/2;
        Satellite.Wavelength_c          =   mean([hypspec_min;hypspec_max]);
        Satellite.Nr_bands              =   length(hypspec_min);
        Satellite.Wavelength            =   wl;        
        
        Satellite.Sensitivity           =   zeros(length(wl),Satellite.Nr_bands);
        for iband = 1:Satellite.Nr_bands
            band_wl                     =   hypspec_min(iband):hypspec_max(iband);
            sensitivity_wl              =   ones(size(band_wl));
            Satellite.Sensitivity(:,iband)=   interp1(band_wl, sensitivity_wl,Satellite.Wavelength,'nearest');
        end
        Satellite.Sensitivity(isnan(Satellite.Sensitivity))=0;

end



%% Creating Variable Space for creation of Lookup Table
[Lai, Lidfa, Lidfb, Hot, Fb, Diss, Cv, Zeta,...
 Cab_g, Cw_g, Cdm_g, Cs_g, N_g,...
 Cab_b, Cw_b, Cdm_b, Cs_b, N_b] =   ndgrid(lai, lidfa, lidfb, hot, fb, diss, cv, zeta,...
                                           cab_g, cw_g, cdm_g, cs_g, n_g,...
                                           cab_b, cw_b, cdm_b, cs_b, n_b);

vegpar                          =   single([Lai(:), Lidfa(:), Lidfb(:), Hot(:), Fb(:), Diss(:), Cv(:), Zeta(:)])';   %ALL
leafgreen                       =   single([Cab_g(:), Cw_g(:), Cdm_g(:), Cs_g(:), N_g(:)])';                         %ALL Green leaf parameter
leafbrown                       =   single([Cab_b(:), Cw_b(:), Cdm_b(:), Cs_b(:), N_b(:)])';                         %ALL Brown leaf parameter

ang_s                           =   single(tts);                                        
ang_o                           =   single([tto, psi])';
soilspec                        =   single(soil_refr);
optipar                         =   single([Leaf_refr; Kdm; Kab; Kw; Ks; H20_refr; H20_abs]);                       %All PROSPECTs parameters
soilpar                         =   single([H_b, H_c, H_B0, H_h, H_SM])';                                            %All Hapke parameters


%% END
return

function [lai, lidfa, lidfb, hot, fb, diss, cv, zeta, ...
         cab_g, cw_g, cdm_g, cs_g, n_g,...
         cab_b, cw_b, cdm_b, cs_b, n_b,...
         H_SM]       =   ReadParameters(path2file)
      
fid                             =   fopen(path2file);

% Read Vegetation parameters
% keyboard
text.lai                        =   fgetl(fid);
isep                            =   find(text.lai=='=')+1;
lai			                    =   str2num(text.lai(isep:end));
text.lidfa                      =   fgetl(fid);
isep                            =   find(text.lidfa=='=')+1;
lidfa		                    =   str2num(text.lidfa(isep:end));
text.lidfb                      =   fgetl(fid);
isep                            =   find(text.lidfb=='=')+1;
lidfb		                    =   str2num(text.lidfb(isep:end));
text.hot                        =   fgetl(fid);
isep                            =   find(text.hot=='=')+1;
hot			                    =   str2num(text.hot(isep:end));
text.fb                         =   fgetl(fid);
isep                            =   find(text.fb=='=')+1;
fb			                    =   str2num(text.fb(isep:end));
text.diss                       =   fgetl(fid);
isep                            =   find(text.diss=='=')+1;
diss		                    =   str2num(text.diss(isep:end));
text.cv                         =   fgetl(fid);
isep                            =   find(text.cv=='=')+1;
cv			                    =   str2num(text.cv(isep:end));
text.zeta                    	=   fgetl(fid);
isep                            =   find(text.zeta=='=')+1;
zeta		                    =   str2num(text.zeta(isep:end));
fgetl(fid);

% Read Green leaves
text.cab_g                     	=   fgetl(fid);
isep                            =   find(text.cab_g=='=')+1;
cab_g		                    =   str2num(text.cab_g(isep:end));
text.cw_g                       =   fgetl(fid);
isep                            =   find(text.cw_g=='=')+1;
cw_g		                    =   str2num(text.cw_g(isep:end));
text.cdm_g                      =   fgetl(fid);
isep                            =   find(text.cdm_g=='=')+1;
cdm_g		                    =   str2num(text.cdm_g(isep:end));
text.cs_g                       =   fgetl(fid);
isep                            =   find(text.cs_g=='=')+1;
cs_g		                    =   str2num(text.cs_g(isep:end));
text.n_g                        =   fgetl(fid);
isep                            =   find(text.n_g=='=')+1;
n_g			                    =   str2num(text.n_g(isep:end));
fgetl(fid);

% Read Brown leaves
text.cab_b                      =   fgetl(fid);
isep                            =   find(text.cab_b=='=')+1;
cab_b		                    =   str2num(text.cab_b(isep:end));
text.cw_b                       =   fgetl(fid);
isep                            =   find(text.cw_b=='=')+1;
cw_b		                    =   str2num(text.cw_b(isep:end));
text.cdm_b                      =   fgetl(fid);
isep                            =   find(text.cdm_b=='=')+1;
cdm_b		                    =   str2num(text.cdm_b(isep:end));
text.cs_b                       =   fgetl(fid);
isep                            =   find(text.cs_b=='=')+1;
cs_b		                    =   str2num(text.cs_b(isep:end));
text.n_b                        =   fgetl(fid);
isep                            =   find(text.n_b=='=')+1;
n_b                             =   str2num(text.n_b(isep:end));
fgetl(fid);

%Soil Moisture
text.H_SM                        =   fgetl(fid);
isep                            =   find(text.H_SM=='=')+1;
H_SM                             =   str2num(text.H_SM(isep:end));


fclose(fid);
return