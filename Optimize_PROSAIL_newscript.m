function [Output]=Optimize_PROSAIL()   
    clc
    close all
    global wl    
    
    global limits
    limits.Cab = [20  120];
    limits.Car = [1e-9 30];
    limits.Ant = [1e-9 10];
    limits.Cm = [1e-9 0.3];
    limits.Cs = [1e-9 10];
    limits.LAI = [1e-9 6];
    limits.Cw = [1e-9 0.3];
    
    
    %% Parameters
    global Cab Car Ant Cs Cw Cm LAI N bsoil
    global selection
    global mu
    
    N                               =   1.6;
    Cab                             =   120;
    Car                             =   40;
    Ant                             =   15;
    Cs                              =   0.0;
    Cw                              =   0.1;
    Cm                              =   0.1;
    LAI                             =   1.6;
    bsoil                           =   0.5;
    
%   uuuuu selection                   =   {'Cab','Car','Ant','Cm','Cw','LAI','Cs'};
%   statevector_init            =   [Cab,Car,Ant,Cm,Cw,LAI,Cs];
    
    selection1                  =   {'Cab','Car','Ant'};
    selection2                  =   {'Cm','Cs','LAI','Cw'};
    
    test                        =   [60.123, 15.4321, 0.01123, 0, 0.2, 0.2, 2.3243];
    
    mu                          =   1e-9;               %dampening factor
    %%
    global tts tto psi 
    tto                         =   30;
    longitude                   =   5.96944;
    latitude                    =   52.21;
    rotation                    =   0;
    dst                         =   1;
    time_zone                   =   +1;

    wl                          =   400:2500;
    global soiltype         
    soiltype                    =   10;
    
    %% define datadir
    path(path,'../LUT retrieval/PROSAIL_D_Matlab/code/')
    dirname                         =   SelectDirectory;

    %% Read Soil spectrum    
    global rsoil
    files                           =   dir([dirname,'/soil/*.sed']);
    for j=1:length(files)
        filename                    =   files(j).name;
        fid                         =   fopen([dirname,'/soil/',filename],'r');
        data                        =   textscan(fid,'%f %f','headerlines',27);
        fclose(fid);

        wv                          =   data{1};
        Refl_soil(:,j)              =   data{2};
    end


    %% Read Canopy spectrum    
    files                       =   dir([dirname,'/*.sed']);
    for j=1:length(files)
        filename                =   files(j).name;
        fid                     =   fopen([dirname,'/',filename],'r');
        comment                 =   fgetl(fid);
        version                 =   fgetl(fid);
        filename                =   fgetl(fid);
        instrument              =   fgetl(fid);
        detectors               =   fgetl(fid);
        meas_type               =   fgetl(fid);
        aq_date                 =   fgetl(fid);
        aq_time                 =   fgetl(fid);
        frewind(fid)

        data                    =   textscan(fid,'%f %f','headerlines',27);
        fclose(fid);

        wv                      =   data{1};
        Refl_can(:,j)           =   data{2};
        date                    =   [aq_date(18:end),' ',aq_time(19:end)];
        Angles(:,j)             =   solarPosition(date,latitude,longitude, time_zone,rotation,dst)';
        
    end

    %% Post processing
    %filter out atmospheric effects
    iatm                        =   wv<399 | (wv>1800 & wv<2000) | (wv>2300);
    Refl_soil(iatm,:)           =   NaN;
    Refl_can(iatm,:)            =   NaN;

    %filter out strange behaviour (possibly also atmospheric, this needs to be investigated...)
    iatm                        =   (wv>1350 & wv<1415) | wv>1340;
    Refl_soil(iatm,:)           =   NaN;
    Refl_can(iatm,:)            =   NaN;

    % average all soil reflectance measurements
    iatm                        =   (wv>1800 & wv<2000) | (wv>2300);
    soil_                       =   nanmean(Refl_soil,2);

    % smooth out any spectral variations
    soilspec                    =   single(conv(soil_,1./(5*ones(5,1)),'same')')/100;
    Canspec                     =   Refl_can/100;

    %% 
%     h11= subplot(2,1,1);
%     plot(wv,Canspec(:,1))
    for imeas=1:size(Canspec,2)
        canspec                     =   Canspec(:,imeas);
        spectraljump                =   canspec(wv==995)-canspec(wv==979);
        if abs(spectraljump)>0.02
            canspec(wv>=995)            =   canspec(wv>=995)-spectraljump;
            canspec(wv>=979 & wv<=995)  =   NaN;
        end
        Canspec(:,imeas)            =   canspec;
        
    end
%     h12= subplot(2,1,2);
%     plot(wv,Canspec(:,1))
%     linkaxes([h11,h12],'xy')
%     keyboard
    %%
    legendstr                       =   selection;
    
%     figure
%     axes('nextplot','add')
    for imeas=1:size(Canspec,2)
        imeas
        tts                         =   Angles(1,imeas);
        psi                         =   Angles(2,imeas);
        canspec                     =   Canspec(:,imeas);
                   
        rcan                        =   interp1(wv,canspec,wl,'nearest')';
        rsoil                       =   interp1(wv,soilspec,wl,'nearest')';
        
        %% only use for testing purposes
%         rdot_test                   =   RTmodel(test, 30, tto, psi, rsoil); 
%         rcan                        =   rdot_test + 1e-2;
        
        for pass=1:2
            %% first retrieve first set of parameters (sensitive mostly to VIS)
            selection                       =   selection2;
            statevector_init                =   zeros(size(selection));
            for j=1:length(selection), 
                statevector_init(j)         =   eval(selection{j}); 
            end
            [statevector,error1]            =   InvertPROSAIL(statevector_init,rcan);

    %         write retrieved parameters to default
            for j=1:length(selection)
                eval([selection{j},sprintf('=%f;',statevector(j))])
            end

            %% afterwards retrieve first set of parameters (sensitive mostly to IR)
            selection                       =   selection1;
            statevector_init                =   zeros(size(selection));        
            for j=1:length(selection), 
                statevector_init(j)         =   eval(selection{j}); 
            end
            [statevector,error2]            =   InvertPROSAIL(statevector_init,rcan);
            
        end
        
        %% save state
        Output.N_(imeas)                           =   N;
        Output.Cab_(imeas)                         =   Cab;
        Output.Car_(imeas)                         =   Car;
        Output.Ant_(imeas)                         =   Ant;
        Output.Cs_(imeas)                          =   Cs;
        Output.Cw_(imeas)                          =   Cw;
        Output.Cm_(imeas)                          =   Cm;
        Output.LAI_(imeas)                         =   LAI;
        Output.bsoil_(imeas)                       =   bsoil;
        Output.error_(imeas)                       =   error1 | error2;
    
%         selection                           =   [selection1, selection2];
%         
%         for j=1:length(selection)
%             eval([selection{j},'_(imeas)=',selection{j},';'])
%         end
%     
    end
    
    %% save output
    selection                                       =   [selection1,selection2];
    outputdir                                       =   [dirname,'/Traits/'];
    if ~exist(outputdir,'dir')
        mkdir(outputdir)        
    end
    save([outputdir,'RetrievedTraitsNewfull4.mat'],'Output','selection','selection1','selection2') 
    
    %% Plot
%     load('2019_05_29\Traits\RetrievedTraits2.mat')
    selection               =   [selection1,selection2];
    Ncol                    =   3;
    Nrow                    =   ceil(length(selection)/Ncol);
    
    h1                      =   figure('Position',[20 20 1024 800]);
    ierror                  =   Output.error_;
    for j=1:length(selection)
        varname             =   [selection{j},'_'];
        V                   =   Output.(varname);
        [Nhist,xhist]       =   hist(V(~ierror),20);
        h11                 =   subplot(Ncol,Nrow,j);
        h111                =   bar(xhist,Nhist);
        
        h114                =   title([varname(1:end-1), sprintf(' (Q_{50}=%4.3f)',median(V(~ierror)))]);
%         h115                =   legend();        
    end
    print(h1,[outputdir,'RetrievedTraitsNewfull4.png'],'-dpng','-r300')
%     keyboard
end
    
%%
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

function [rdot]=RTmodel(statevector, tts, tto, psi, rsoil)
    global Cab Car Ant Cs Cw Cm LAI N bsoil
    global selection
    for j=1:length(selection)
%         assignin('base', selection{j}, statevector(j))
        eval([selection{j},sprintf('=%f;',statevector(j))])
    end
    % execute prosail
    [rdot,~,~,~]       =   PRO4SAIL(N,Cab,Car,Ant,Cs,Cw,Cm,...
                                                -0.35,-0.15,1,LAI,0.05,tts,tto,psi,rsoil*bsoil);
    %%
%     close all
%     global wl
%     figure
%     for ant=00:5:20
%         for cab=40:5:80
%             for car=0:5:20
%                 [rdot,rsot,rddt,rsdt]       =   PRO4SAIL(N,cab,car,ant,Cs,Cw,Cm,...
%                                                     -0.35,-0.15,1,LAI,0.05,tts,tto,psi,rsoil*0.2);
%                 hold on
%                 plot(wl,rdot)
%             end
%         end
%     end
% 
%     figure
%     for cw=linspace(0.01,0.1,5)
%         for cm=linspace(0.01,0.1,5)
%             for cs=0:5:10
%                 [rdot,rsot,rddt,rsdt]       =   PRO4SAIL(N,Cab,Car,Ant,Cs,cw,cm,...
%                                                     -0.35,-0.15,1,LAI,0.05,tts,tto,psi,rsoil*0.2);
%                 hold on
%                 plot(wl,rdot)
%             end
%         end
%     end

end

function [statevector,error]=InvertPROSAIL(statevector_init,rcan)
    lastwarn('')
    global tts tto psi rsoil
    global selection limits
    global wl
    global converge
    global mu
    %% Parameters
    maxiter                     =   50;
    converge                    =   1e-4;
    
    %% Initialization
    statevector                 =   statevector_init;
    dx                          =   statevector*0 + Inf;
    iter                        =   0;        
    cont                        =   1;
    error                       =   0;
    %% Retrieval
    plotoption                  =   1;
    if plotoption==1;
        figure
    end
	
    while cont==1;
        dx_old                  =   dx;
        iter                    =   iter + 1;

        [rdot]                  =   RTmodel(statevector, tts, tto, psi, rsoil);

        for j=1:length(statevector)
            statevector_prime   =   statevector;
            statevector_prime(j)=   statevector_prime(j)*(1+1e-1) +1e-6;
            ds                  =   statevector_prime(j)-statevector(j);
            
            rdot_prime          =   RTmodel(statevector_prime, tts, tto, psi, rsoil);
            dr_prime            =   rdot_prime-rdot;
            dr_dx               =   dr_prime/ds;
            J(:,j)              =   dr_dx';
        end

        %% correct for bias using the difference in blue band
        iwlbias             =   450;
        bias                =   rcan(wl==iwlbias) - rdot(wl==iwlbias);
        rcan_bias           =   rcan-bias;
        dr                  =   rdot-rcan_bias;

% %             perhaps use continuum removal?             
%             [~, ipks] = findpeaks(double(rdot));
%             iredge      =   find(wl>730 & wl<750);
%             ipks        =   unique([ipks(:); iredge(:)]);                        
%             cc          =   rdot - interp1(wl(ipks), rdot(ipks), wl,'linear')';
%             plot(wl(ipks), rdot(ipks),'o', wl, rdot)
%             ccr                     =   (cc-refl_avg(:,j,jj)')/100;


        % now apply jacobian on difference Obs with Sim. 
        dr(isnan(dr))=0;
        J(isnan(J)) = 0;

%           According to Verhoef 2000, Simultaneous retrieval of soil, leaf, canopy and 
%           atmospheric parameters from hyperspectral information in the red edge through 
%           model inversion

%             h11 = subplot(2,1,1);
%             plot(wl,rcan, wl, rdot)
%             h12 = subplot(2,1,2);
%             plot(wl,J)
%             linkaxes([h11,h12],'x')

        M   =   (J'*J + mu* eye(length(statevector)));
        dx  =    M^-1*J'*dr;
        
        [q1,q2] = lastwarn;
        switch q2
            case {'MATLAB:nearlySingularMatrix','MATLAB:illConditionedMatrix'}
                if max(abs(dx)>100)
                    error   =   1;
                  %  keyboard
                end
%                 plot(wl,rdot, wl, rcan-bias)
%                 keyboard
                
        end
        % update state vector
        statevector                 =   statevector - dx';
        statevector                 =   max(statevector,1e-9);
        statevector                 =   min(statevector,1.2e2);
        
        if plotoption==1;
            subplot(2,1,1,'nextplot','add')
            plot(wl,rcan-bias,'r','linewidth',2)
            plot(wl,rdot,'b')
            drawnow
            pause(0.5)
            plot(wl,rdot,'g')
            subplot(2,1,2,'nextplot','add')
            plot(iter,statevector,'.')
            legend(selection)
            drawnow
            
            statevector
%             keyboard
        end
        
        
%         keyboard
        %% set values to be within limits
        for iselec = 1:length(selection)
        limit = limits.(selection{iselec});
        statevector(iselec) = max(statevector(iselec),limit(1));
        statevector(iselec) = min(statevector(iselec),limit(2));
        end
        %% Continue
        cont1   =   max(abs(dx))>converge;
        cont2   =   iter<maxiter;
        cont3   =   max(abs(dx_old(:)-dx(:)))>0.001;

        cont    =   cont1 * cont2 * cont3;        
    end

%     keyboard
    if nansum(abs(dr))>10
        error =1;
    end
    
    if nanmean((abs(dr)./rcan_bias)*100)>30
        error = 1;
    end
%     rcan_bias

    %%
    if plotoption
        subplot(2,1,1)
        plot(wl,rdot,'b','linewidth',2)
    end
    
end