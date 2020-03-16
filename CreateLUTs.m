clear all
close all
fclose all;
clc


%% Parameters
tto                         =   30;

%% define datadir
dirname         = SelectDirectory;

%% Read Soil spectrum    
files           =   dir([dirname,'/soil/*.sed']);
for j=1:length(files)
    filename    =   files(j).name;
    fid         =   fopen([dirname,'/soil/',filename],'r');
    data        =   textscan(fid,'%f %f','headerlines',27);
    fclose(fid);
    
    wl          =   data{1};
    Refl(:,j)   =   data{2};
end
%filter out atmospheric effects
iatm            =   (wl>1800 & wl<2000) | (wl>2300);
Refl(iatm,:)    =   NaN;

% average all soil reflectance measurements
soil_           =   nanmean(Refl,2)/100;

% smooth out any spectral variations
soilspec           =   single(conv(soil_,1./(5*ones(5,1)),'same')');
% plot(wl,soilspec)

%% check if 'parameter exists'
lut_dir= [dirname,'\LUT\'];
if ~exist(lut_dir,'dir')
    warndlg('No LUT-parameters defined of requested days')
    return
else
    Parameter_file  =    [lut_dir,'Parameters.txt'];
end

%% cloudiness
cloudiness                      = questdlg('How cloudy was it','cloudiness','Blue sky','Some clouds','Overcast','Blue sky');
switch lower(cloudiness)
    case 'blue sky'
        % full BRDF reflectance
        ibrdf                   =   4;
    case 'some clouds'
        % full BRDF reflectance
        ibrdf                   =   [1,3];
    case 'overcast'
        % full HDRF reflectance
        ibrdf                   =   3;
    otherwise
end

%% Angles
angles_dir  =   [dirname,'/Angles/'];
angles      =   load([angles_dir,'Angles.mat'],'-ascii');
zenith      =   angles(1);
azimuth     =   angles(2);
%% Setup Parameters for LUT creation
LUT_dir = '../LUT retrieval/SLC_create_LUT/slc_code/';
path(path,LUT_dir)


global sensor bandwidth minwl maxwl
minwl           =   min(wl);
maxwl           =   max(wl);
bandwidth       =   unique(diff(wl));
sensor          =   'hyperspectral';       %'Rapideye', %'Spot', %'modis'

[  option                      ,...
   wl_,nwl                     ,...
   optipar_,soilspec_,soilpar_,...
   Leafgreen,Leafbrown,Vegpar  ,...
   ~, ~,Satellite]     =   input_parameters_meas(Parameter_file);

N                           =   size(Vegpar,2);

% ierror                      = isnan(soilspec_)
% keyboard
%% Allocating memory

%% Allocating memory
fprintf(1,'Allocating Memory\n')
Rsoil                           =   single(zeros(N,4,nwl));
RTleaf                          =   single(zeros(N,4,nwl));
Rcan                            =   single(zeros(N,6,nwl));
fvc                             =   single(zeros(N,1));

%% Running SLC
fprintf(1,'\nRunning Program\n')
for i=1:N
    tts                         =   zenith;
    psi                         =   azimuth;
 
    leafgreen_                  =   Leafgreen(:,i);
    leafbrown_                  =   Leafbrown(:,i);
    vegpar_                     =   Vegpar(:,i);
    
%     wv                          =   400:1:2500;
    [Rsoil(i,:,:)   ,...
     RTleaf(i,:,:)  ,...
     Rcan(i,:,:)]               =   SLC2(option,nwl,...
                                        optipar_,soilspec,soilpar_    ,...
                                        leafgreen_,leafbrown_,vegpar_  ,...
                                        [tts; tto; psi]);	

	%Post processing
    Rsot                        =   permute(Rcan(i,ibrdf,:),[3 1 2]);
    Rsot(isnan(Rsot))           =   0;
    
    %sensor simulator
    for iband=1:Satellite.Nr_bands
        BandsReflectance(iband,i)=   sum(Rsot.*Satellite.Sensitivity(:,iband).*Satellite.Bandwidth)./sum(      Satellite.Sensitivity(:,iband).*Satellite.Bandwidth +1e-9);
    end    
end

%%
% % Show difference in soil/leaf/canopy reflectance
% rsoil     = interp1(wl_,permute(Rsoil(1,1,:),[3 1 2]), wv);
% rleaf     = interp1(wl_,permute(RTleaf(1,1,:),[3 1 2]), wv);
% rcan      = interp1(wl_,permute(Rcan(:,1,:),[3 1 2]), wv);
% plot(wv,rsoil,'.r'5, ...
%      wv,rleaf,'.g', ...
%      wv,rcan,'.b')

%% save output
save([lut_dir,'LUT_SLC.mat'],'Leafgreen','Leafbrown','Vegpar','BandsReflectance','wl') 

%% plot
plot(wl,BandsReflectance)