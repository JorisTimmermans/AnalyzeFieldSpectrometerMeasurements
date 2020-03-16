function [] = InvestigateMeasurements()
    %% parameters
    wl                      =   350:10:2500;
    longitude               =   5.96944;
    latitude                =   52.21;
    rotation                =   0;
    dst                     =   1;
    time_zone               =   +1;

    %% Rename directories; 
    measurement_days        =   dir('*_*');
    measurement_days        =   measurement_days([measurement_days.isdir]);
    Ndays                   =   length(measurement_days);

    for i=1:Ndays
        string              =   measurement_days(i).name;
        string_new          =   strsplit(string,'_');

        month_str           =   {'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'};
        month               =   find( strcmp(month_str,lower(string_new{2})) );

        if ~isempty(month)
            string_new{2}       =   sprintf('%02.0f',month);
            string_new      =   sprintf('%s_%s_%s',string_new{1:3});
            movefile(string,string_new,'f');
        end
    end


    %% sort folders on date
    measurement_days        =   dir('*_*');
    measurement_days        =   measurement_days([measurement_days.isdir]);
    Ndays                   =   length(measurement_days);

    dates                   =   1:Ndays;
    for i=1:Ndays
        string              =   measurement_days(i).name;
        string(string=='_') =   '-';
        dates(i)            =   datenum(string);
    end
    [~,is]                  =   sort(dates);
    measurement_days        =   measurement_days(is);

    %% Plot average values for each folder
    figure('Position',[20 20 1424 800])
    colors                  =   lines(Ndays);
    for i = 1:Ndays
        measurement_day     = measurement_days(i).name;

        acquisitions        =   dir([measurement_day,'/*.sed']);
        Nacq                =   length(acquisitions);

        if Nacq>1
            Refl            =   zeros(Nwl, Nacq)*NaN;
            Zenith          =   zeros(1, Nacq)*NaN;
            Azimuth         =   zeros(1, Nacq)*NaN;
            for j=1:Nacq
                acquisition =   acquisitions(j).name;
                path2file   =   [measurement_day,'/',acquisition];

                fid         =   fopen(path2file,'r');
                comment     =   fgetl(fid);
                version     =   fgetl(fid);
                filename    =   fgetl(fid);
                instrument  =   fgetl(fid);
                detectors   =   fgetl(fid);
                meas_type   =   fgetl(fid);
                aq_date     =   fgetl(fid);
                aq_time     =   fgetl(fid);
                frewind(fid)
                data        =   textscan(fid,'%f %f','headerlines',27);
                fclose(fid);

    %             date_s      =   datenum([aq_date(07: 16),' ',aq_time(07:17 )]);
    %             date_n      =   datenum([aq_date(18:end),' ',aq_time(19:end)]);
    %             
    %             datestr     =   datevec(date_s);
    %             datestr     =   [aq_date(07: 16),' ',aq_time(07:17 )]
                datestr     =   [aq_date(18:end),' ',aq_time(19:end)];
                angles      =   solarPosition(datestr,latitude,longitude, time_zone,rotation,dst);
                zenith      =   angles(1);
                azimuth     =   angles(2);

                wv          =   data{1};
                refl        =   data{2};

                [wv,is]     =   unique(wv);
                refl        =   refl(is);

                Refl(:,j)   =   interp1(wv,refl,wl);
                Zenith(j)   =   zenith;
                Azimuth(j)  =   azimuth;
            end

            % Filter white reflectances
            iwhite          =   mean((Refl>50))>0.8;
            Refl(:,iwhite)  =   NaN;

            % Filter too large reflectances
            ierror          =   Refl>100;
            Refl(ierror)    =   NaN;

            Refl_q          =   zeros(length(wl),3)*NaN;
            for iwl=1:length(wl)
                r           =   Refl(iwl,:);
                ierror      =   isnan(r);
                if any(~ierror)
                    Refl_q(iwl,1) =   min(r(~ierror));
                    Refl_q(iwl,2) =   median(r(~ierror));
                    Refl_q(iwl,3) =   min(r(~ierror));
                end
             end
    %         Refl_q          =   quantile(Refl,[0.1, 0.5, 0.90],2);

            subplot(2,2,1:2,'nextplot','add')
            isubset             =   1:5:Nwl;
            plot(wl(isubset),Refl_q(isubset,2), 'Color',colors(i,:))
            xlim([350 1300])
            title('Intercomparison of measurements in All folders')

        end
    end

    % Visualize the Spectrum of last folder
    plot(wl(isubset),Refl_q(isubset,2),'o', 'Color',colors(i,:))
    legend({measurement_days.name},'location','best')

    %% Plot all measurements
    subplot(2,2,3:4,'nextplot','add')
    plot(wl,Refl,'linewidth',2)
    xlim([350 1300])
    title(['Intercomparison of ', measurement_day])
    % legend({acquisitions.name},'location','best','fontsize',6)
    
    refl_dir  =   [measurement_day,'/Refl/'];
    if ~exist(refl_dir,'dir')
        mkdir(refl_dir)
    end
    save([refl_dir,'Reflectances.mat'],'wl','Refl')

    %% Investigate Angles
    zenith  = nanmean(Zenith);
    azimuth = nanmean(Azimuth);
    
    try
        figure('Position',[20 20 1424 800])
        subplot(2,1,1)
        plotyy(1:length(Zenith),Zenith, 1:length(Zenith), Azimuth)
        legend(sprintf('Zenith_{avg}=%3.2f',zenith), sprintf('Azimuth_{avg}=%3.2f',azimuth))
        title('Solar Angles')

        subplot(2,1,2)
        Az  =   (90+Azimuth)/180*pi; %Matlab orients north 90s off. 
        El  =   (90-Zenith)/180*pi;
        R   =   ones(size(Azimuth));
        [x,y,z] = sph2cart(Az,El,R); 
        plot(x,y,'yo','Markersize',10,'MarkerFacecolor','R'), 
        axis equal tight
        axis([-1 1 -1 1])
    end
    angles_dir  =   [measurement_day,'/Angles/'];
    if ~exist(angles_dir,'dir')
        mkdir(angles_dir)
    end
    angles = [zenith,azimuth];
    save([angles_dir,'Angles.mat'],'angles','-ascii')
    
    %% step by step comparison
    error_dir   =   [measurement_day,'/','error'];
    soil_dir    =   [measurement_day,'/','soil'];
    if ~exist(error_dir,'dir')
        mkdir(error_dir)
        mkdir(soil_dir)

        figure('Position',[20 20 1424 800])
        for j=1:size(Refl,2)
            hold off
            plot(wl,Refl,'k'), hold on
            plot(wl,Refl(:,j),'r','linewidth',2)
            title(acquisitions(j).name)
            ylim([0,100])

            choice = questdlg('Was this a proper measurements?', ...
                                 'Acquisition quality', ...
                                 'Yes-Canopy','Yes-Soil', 'No', 'Yes-Canopy');
            switch lower(choice)
                case 'stop'
                    break
                case 'yes-soil'    
                    try
                        %
                        filename        =   acquisitions(j).name;
                        filename_meta   =   [filename(1:end-3),'raw'];
                        movefile([measurement_day,'/',filename     ],[soil_dir,'/',filename])
                        movefile([measurement_day,'/',filename_meta],[soil_dir,'/',filename_meta])
                    end

                case 'no'
                    try
                        %
                        filename        =   acquisitions(j).name;
                        filename_meta   =   [filename(1:end-3),'raw'];
                        movefile([measurement_day,'/',filename     ],[error_dir,'/',filename])
                        movefile([measurement_day,'/',filename_meta],[error_dir,'/',filename_meta])
                    end
            end
        end
    end

end

function [angles,projection] = solarPosition(datetime,latitude,longitude, ...
                                             time_zone,rotation,dst)
    %SOLARPOSITION Calculate solar position using most basic algorithm
    %   This is the most basic algorithm. It is documented in Seinfeld &
    %   Pandis, Duffie & Beckman and Wikipedia.
    %
    % [ANGLES,PROJECTION] = SOLARPOSITION(DATE,TIME,LATITUDE,LONGITUDE,TIME_ZONE)
    % returns ZENITH & AZIMUTH for all DATE & TIME pairs at LATITUDE, LONGITUDE.
    % ANGLES = [ZENITH,AZIMUTH] and PROJECTION = [PHI_X, PHI_Y]
    % PHI_X is projection on x-z plane & PHI_Y is projection on y-z plane.
    % DATETIME can be string, vector [YEAR, MONTH, DAY, HOURS, MINUTES, SECONDS],
    %   cellstring or matrix N x [YEAR, MONTH, DAY, HOURS, MINUTES, SECONDS] for N
    %   times.
    % LATITUDE [degrees] and LONGITUDE [degrees] are the coordinates of the site.
    % TIME_ZONE [hours] of the site.
    % ROTATION [degrees] clockwise rotation of system relative to north.
    % DST [logical] flag for daylight savings time, typ. from March to November
    %   in the northern hemisphere.
    %
    % References:
    % http://en.wikipedia.org/wiki/Solar_azimuth_angle
    % http://en.wikipedia.org/wiki/Solar_elevation_angle
    %
    % Mark A. Mikofski
    % Copyright (c) 2013
    %
    %% datetime
    datetime = datenum(datetime); % [days] dates & times

    date = floor(datetime); % [days]
    [year,~,~] = datevec(date);
    time = datetime - date; % [days]

    %% constants
    toRadians = @(x)x*pi/180; % convert degrees to radians
    toDegrees = @(x)x*180/pi; % convert radians to degrees
    %% Equation of time
    d_n = mod(date-datenum(year,1,1)+1,365); % day number
    B = 2*pi*(d_n-81)/365; % ET parameter
    ET = 9.87*sin(2*B)-7.53*cos(B)-1.5*sin(B); % [minutes] equation of time
    % approximate solar time
    solarTime = ((time*24-double(dst))*60+4*(longitude-time_zone*15)+ET)/60/24;
    latitude_rad = toRadians(latitude); % [radians] latitude
    rotation_rad = toRadians(rotation); % [radians] field rotation
    t_h = (solarTime*24-12)*15; % [degrees] hour angle
    t_h_rad = toRadians(t_h); % [radians]
    delta = -23.45 * cos(2*pi*(d_n+10)/365); % [degrees] declination
    delta_rad = toRadians(delta); % [radians]
    theta_rad = acos(sin(latitude_rad)*sin(delta_rad)+ ...
        cos(latitude_rad)*cos(delta_rad).*cos(t_h_rad)); % [radians] zenith
    theta = toDegrees(theta_rad); % [degrees] zenith
    elevation = 90 - theta; % elevation
    day = elevation>0; % day or night?
    cos_phi = (cos(theta_rad)*sin(latitude_rad)- ...
        sin(delta_rad))./(sin(theta_rad)*cos(latitude_rad)); % cosine(azimuth)
    % azimuth [0, 180], absolute value measured from due south, so east = west = 90,
    % south = 0, north = 180
    phi_south = acos(min(1,max(-1,cos_phi)));
    % azimuth [0, 360], measured clockwise from due north, so east = 90,
    % south = 180, and west = 270 degrees
    phi_rad = NaN(size(phi_south)); % night azimuth is NaN
    % shift from ATAN to ATAN2, IE: use domain from 0 to 360 degrees instead of
    % from -180 to 180
    phi_rad(day) = pi + sign(t_h(day)).*phi_south(day); % Shift domain to 0-360 deg
    % projection of sun angle on x-z plane, measured from z-direction (up)
    phi_x = toDegrees(atan2(sin(phi_rad-rotation_rad).*sin(theta_rad), ...
        cos(theta_rad))); % [degrees]
    % projection of sun angle on y-z plane, measured from z-direction (up)
    phi_y = toDegrees(atan2(cos(phi_rad-rotation_rad).*sin(theta_rad), ...
        cos(theta_rad))); % [degrees]
    phi = toDegrees(phi_rad); % [degrees] azimuth
    angles = [theta, phi]; % [degrees] zenith, azimuth
    projection = [phi_x,phi_y]; % [degrees] x-z plane, y-z plane
end
