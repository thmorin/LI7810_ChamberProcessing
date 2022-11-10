function [HM_CO2, NRMSE_CO2, HM_CH4,NRMSE_CH4,BubbleFlux]=ProcessLI7810Chamber(meas_date,start_time,end_time,L_data_sheet,Start_offset,End_offset,ChamberVolume,ChamberArea, plotYN,Site,Start_temp,End_temp,Pressure,~,~)
% ProcessLI7810Chamber - Written by Tim Morin, 7/5/2022
% Added ebullition and HM fit, assuming Picarro laser use - Erin Hassett and Tim Morin 7/26/2022
% Adapted for LI-7810 laser use again. Cleaned up code - Tim Morin 11/9/2022
% 
% meas_date - provide as a datetime, because that is how Picarro will read the time it must compare it to
% nb_start_time - datetime for start of chamber measurement, recorded in notebook with stopwatch
% nb_end_time - datetime for end of chamber measurement, recorded in notebook with stopwatch
% Start_temp - degrees C - Temp at start of measurement period
% End_temp - degrees C - Temp at end of measurent period
% start_time- epoch time recorded by picarro (Seconds since 1/1/1970)
% stop_time- epoch time recorded by picarro
% Start_offset - How many points to ignore at the start of the measurement period
% End_offset - How many points to ignore at the end of the measurement period
% Pressure - Pa - Ambient pressure during measurement
% Offset - Time measurement offset to GMT= unused currently

%% Import raw Picarro data
laserdata=ReadLI7810(L_data_sheet);

%% filter to relevant day
dayFilt=laserdata.DATE==meas_date;
laserdata.TIME(~dayFilt)=NaT;

%% Remove bad diagnostic data. NOTE: Picarro alarm doesn't necessarily seem to indicate bad data
% bad=laserdata.DIAG~=0;
% laserdata.CO2_dry(bad | ~dayFilt)=nan;
% laserdata.CH4_dry(bad | ~dayFilt)=nan;

%% Filter to relevant time of measurment
tfilt = laserdata.TIME>=start_time & laserdata.TIME<=end_time;%in between start time and end time
TIME=laserdata.SECONDS(tfilt);
TIME=TIME-TIME(1);

CO2=laserdata.CO2(tfilt);
CH4=laserdata.CH4(tfilt);
H2O=laserdata.H2O(tfilt);
%% Correct for chamber volume and area
ChamberVolume=ChamberVolume/100/100/100;%originally in cm3
ChamberArea=ChamberArea/100/100; %conversion from cm2

%% Ideal gas law correction
TMP=(linspace(Start_temp,End_temp,length(TIME))+273.15)';   %Conversion to K, puts into vector
R=8.314;                  %m3 Pa K-1 mol-1
CO2=CO2*Pressure/R./TMP;  %umol/m^3
CH4=CH4*Pressure/R./TMP;  %umol/m^3

%% Trim out any data at beginning or end of measurement period that were affected by chamber settling
[~,IndSt]=min(abs(TIME-(Start_offset+1)));
[~,IndEn]=min(abs(TIME-(End_offset+1))); 

T_TIME=TIME(IndSt:(length(TIME)-IndEn));
% T_H2O =H2O (Start_offset+1:(length(TIME)-End_offset));
T_CO2 =CO2 (IndSt:(length(TIME)-IndEn));
T_CH4 =CH4 (IndSt:(length(TIME)-IndEn)); 

if length(T_CH4) < 50 %if you only had 50 points for a sample, don't count it
    %% Too few points to run regression on
    HM_CH4 = nan;
    HM_CO2=nan;
    BubbleFlux=nan;
    NRMSE_CH4= nan;
    NRMSE_CO2=nan;
else
    %% Detect, sum up, and remove bubbles
    BubThr = Inf; %August 2022 - TIM: NEED TO DIAL THIS IN;

    CHDif = diff(T_CH4);  %creates matrix of changes in methane concentration as the sample progresses; 
    % fprintf('Chamber %s had max diff of %0.03f\n',Site,max(CHDif));
    BubbleAdj = (diff(T_CH4)./ diff(T_TIME) > BubThr) | (diff(T_CH4)./diff(T_TIME) < -BubThr); %makes 1 or 0 if above or below threshold (| = or);creates a matrix of when bubbles cross the threshold in either (+) or (-) direction

    CH4Add = zeros(size(CHDif));  %matrix with zeros + same size as CHDif
    CH4Add(BubbleAdj) = CHDif(BubbleAdj);%pairs bubbleadjust with CHDif; when they have same paired event, 
                                           %it puts it into CH4add; (i.e. when
                                           %Bubble Adj =1; pulls value from
                                           %CHDiff)= all methane above thresh
    TotBubble = cumsum(CH4Add);%adds up the total amount of methane attributed to bubbling events for the sample
    T_CH4_NoBubbles = [T_CH4(1); T_CH4(2:end) - TotBubble]; %new matrix; T_CH4 start and stop; starts at 1 and goes from 2:end.
                                                            %creates a matrix of methane exluding CH4 from bubbling events= total diffusive flux
    BubbleFlux = (TotBubble(end)/(T_TIME(IndEn) - T_TIME(IndSt)))*ChamberVolume/ChamberArea; % umol m-2 s-1
    %%  METHANE - Nonlinear fit: Hutchinson-Mosier fit one-dimensional fit for diffusive fluxes
    % Assumes that concentrations will approach an asymptote as they draw near the limit where gas will go back into surface as fast as emission
    
    % Remove CH4 outliers
    T_CH4_rm=T_CH4;
    [~,TF] = rmoutliers(detrend(T_CH4_rm), "mean");
    T_CH4_rm(TF==1)=nan;
    
    % Initial conditions for HM parameters
    Start_co_CH4 = mean(T_CH4_rm(1:5),'omitnan'); % Initial concentration value
    Start_cs_CH4 = max(T_CH4_rm);                 % Saturation concentration value
    Start_k_CH4 = 0.0035;                         % First order-rate constant
    x0=[Start_cs_CH4,Start_co_CH4,Start_k_CH4];
    
    x=nlinfit(T_TIME(~isnan(T_CH4_rm)),T_CH4_rm(~isnan(T_CH4_rm)),@myfun2,x0);
    %x = lsqcurvefit(@myfun2,x0,T_TIME(~isnan(T_CH4_rm)),T_CH4_rm(~isnan(T_CH4_rm))); %Apply curve fit
    
    % Post fit values
    cs_CH4=x(1);
    co_CH4=x(2);
    k_CH4=x(3);
    
    % Determine flux - Assumes first derivative at time zero is diffusive flux in the absence of chamber presence
    HM_CH4 = -k_CH4.*(co_CH4 - cs_CH4)*ChamberVolume/ChamberArea;

    %%  CARBON DIOXIDE - Nonlinear fit: Hutchinson-Mosier fit one-dimensional fit for diffusive fluxes
    % Similar assumption, but for negative fluxes roughly approximates that plants will have a hard time taking up past the asymptotic value
    
    % Remove CO2 outliers
    T_CO2_rm=T_CO2;
    [~,TF] = rmoutliers(detrend(T_CO2), "mean");
    T_CO2_rm(TF==1)=nan;
    
    % Initial conditions for HM parameters
    Start_co_CO2 = mean(T_CO2_rm(1:5),'omitnan'); % Initial concentration value
    Start_cs_CO2 = 0;                             % Saturation concentration value
    Start_k_CO2 = 0;                              % First order-rate constant
    z0=[Start_cs_CO2,Start_co_CO2,Start_k_CO2];
    
    
    z=nlinfit(T_TIME(~isnan(T_CO2_rm)),T_CO2_rm(~isnan(T_CO2_rm)),@myfun3,z0);
%     z = lsqcurvefit(@myfun3,z0,T_TIME(~isnan(T_CO2_rm)),T_CO2_rm(~isnan(T_CO2_rm)));
    
    % Post fit values
    cs_CO2=z(1);
    co_CO2=z(2);
    k_CO2=z(3);

    % Determine flux - Assumes first derivative at time zero is diffusive flux in the absence of chamber presence
    HM_CO2 =-k_CO2.*(co_CO2 -cs_CO2)*ChamberVolume/ChamberArea; 

    % see Pendersen et al.2000 (https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2389.2010.01291.x) and Kutzbach et al. 2007(https://www.biogeosciences.net/4/1005/2007/)

    %% Optionally plot the chamber raw data and the fit
    if plotYN==1
        f1=figure();
        f1.Units='Inches';
        f1.Position=[2 2 4 4];

        %% Useful lines to plot what you're doing with
        subplot(2,1,1);
        y=myfun2(x,T_TIME);
        y2=myfun3(z,T_TIME);

        %%RMSE without outliers
        %CH4 RMSE
        err= y - T_CH4_rm;
        squareError= err.^2;
        sumSqaureError= sum(squareError, 'omitnan'); %added ~isnan= doesn't work
        n=sum(~isnan(T_CH4_rm));
        RMSE_CH4=sqrt(sumSqaureError/n);
        NRMSE_CH4=RMSE_CH4/mean(T_CH4_rm, 'omitnan')*100; 

        %CO2 RMSE
        err_CO2= y2 - T_CO2_rm;
        squareError_CO2= err_CO2.^2;
        sumSqaureError_CO2= sum(squareError_CO2,'omitnan'); % keeps returning NANs
        n_CO2=sum(~isnan(T_CO2_rm));
        RMSE_CO2=sqrt(sumSqaureError_CO2/n_CO2);
        NRMSE_CO2=RMSE_CO2/mean(T_CO2_rm, 'omitnan')*100; %added omitnan

        %CH4 plots without outliers
        hold on;
        plot(T_TIME,y,'-r','MarkerSize',15); %this gives the curve fit with HM line (should pair with no bubbles plot)
        plot(T_TIME,T_CH4_rm,'.c','MarkerSize',10); %This is the plot with no bubbles 
        errorbar(T_TIME,T_CH4_rm,RMSE_CH4*ones(length(T_CH4_rm),1));
        title(Site);
        ylabel('CH_4 (ppm)');
        legend ({'H-M Fit','Diffusive Flux'},... %
             'FontSize',18,'TextColor','black' );

        %CO2 without outliers
        subplot(2,1,2);
        hold on;
        plot(T_TIME,y2,'-r','MarkerSize',15); %CO2 HM line
        plot(T_TIME,T_CO2_rm,'.c','MarkerSize',10);
        ylabel('CO_2 (\mumol m^{-3})');
        xlabel('Time (s)');
        errorbar(T_TIME,T_CO2_rm,RMSE_CO2*ones(length(T_TIME),1));
        legend ({'H-M Fit','No Outliers',},... 
            'FontSize',18,'TextColor','black' );
    end
end
end
function F =myfun2(x,xdata)
F=x(1)+(x(2)-x(1)).*exp(-x(3).*xdata);
end
function F2 =myfun3(z,xdata)
F2=z(1)+(z(2)-z(1)).*exp(-z(3).*xdata);
end