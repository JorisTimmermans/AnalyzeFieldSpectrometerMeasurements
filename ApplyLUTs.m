clear all
close all
fclose all;
clc

%% define datadir
dirname         = SelectDirectory;

%% Read Canopy spectrum    
files           =   dir([dirname,'/*.sed']);
for j=1:length(files)
    filename    =   files(j).name;
    fid         =   fopen([dirname,'/',filename],'r');
    data        =   textscan(fid,'%f %f','headerlines',27);
    fclose(fid);
    
    wl          =   data{1};
    Refl(:,j)   =   data{2}/100;
end

%filter out atmospheric effects
iatm            =   (wl>1800 & wl<2000) | (wl>2300);
Refl(iatm,:)    =   NaN;

%% load LUT data
lut_dir         =   [dirname,'\LUT\'];
if ~exist([lut_dir,'LUT_SLC.mat'],'file')
    warndlg('No LUT file found, please run "CreateLUTs.m" ')
    return
end
load([lut_dir,'LUT_SLC.mat'],'BandsReflectance','Vegpar','Leafbrown','Leafgreen')

%% plot data
subplot(2,1,1)
plot(wl,Refl)
title('Measurements')
subplot(2,1,2)
plot(wl,BandsReflectance)
title('LUT')

nLUT = size(BandsReflectance,2);
nwl             =   length(wl);
filter          =   ones(nwl,1);
for j=1:size(Refl,2)
    can_refl    =   Refl(:,j);
    Can_Refl    =   can_refl(:, ones(1,nLUT))*pi;
    
    
%     plot(wl,BandsReflectance,'b',wl,can_refl,'r')
%     keyboard
    E           =   BandsReflectance - Can_Refl;
    
    RMSE        =   sqrt(nanmean(E(:,:).^2));
    [~,imin]    =   min(RMSE);
    
%     rMSE        =   sqrt(nanmean(E./Can_Refl.^2));
%     [~,imin]    =   min(RMSE);

    for ilut=1:size(BandsReflectance,2)
        c           =   corrcoef(BandsReflectance(50:1440,ilut), can_refl(50:1440,1));
        C(ilut)     =   c(1,2);
    end
    [~,imin]    =   max(C);
    
%     keyboard
    leafbrown   =   Leafbrown(:,imin);
    leafgreen   =   Leafgreen(:,imin);
    vegpar      =   Vegpar(:,imin);
    
    lai(j)         =   vegpar(1);
%     lidfa       =   vegpar(2);
%     lidfb       =   vegpar(3);
%     hot         =   vegpar(4);
%     fbrown      =   vegpar(5);
%     diss        =   vegpar(6);
%     clumping    =   vegpar(7);
%     zeta        =   vegpar(8);
    
    cab(j)         =   leafgreen(1);
    cw(j)          =   leafgreen(2);
    cdm(j)         =   leafgreen(3);
    cs(j)          =   leafgreen(4);
    N(j)           =   leafgreen(5);
end

%%
figure
subplot(3,1,1)
plot(cab)
ylabel('cab')
subplot(3,1,2)
plot(lai)
ylabel('lai')
subplot(3,1,3)
plot(cw)
ylabel('cw')
% save([lut_dir,'LUT_SLC.mat'],'Leafgreen','Leafbrown','Vegpar','BandsReflectance','wl') 

% save(