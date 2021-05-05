clc; clear all; 
close all;

%% 1. Specify Experimental Properties

nExpID      = 'ZNH1';
nSaveID     = '21-04-27.mat';  % Contains T_save
nSample     = '21-04-27';
nRef        = 'H2O';

dt          = 10;                        % Sampling time steps in [s]
texp        = 4;                       % Time for a single cooling/heating experiment in [h]

Tamb_max = 45;
Tamb_min = 24;

msample     = 10e-3;                  % in [kg] 100g topp
mref        = 17.2e-3;                  % in [kg]
mtsample    = 168.2e-3;                  % in [kg]
mtref       = 139.7e-3;                  % in [kg]
nPlot       = 'k-';

idx_sensor = {'c';...     % Top Center Bottom
                 9;...    % Reference readings
                3};       % Sample readings

NS = 1;     % No. of sensors per sample
NC = 1;     % No. of cycles

       
%% 2. Plot
load(nSaveID);          % load raw_readings (import from Excel)
nType_cool = 'c';  % Part of filename 
nType_heat = 'h';  % Part of filename 
Untitled=Untitled(1:1:46931,:);
figure
%plot(Untitled(:,2)), hold on
plot(Untitled(:,3)), hold on
plot(Untitled(:,4)), hold on
legend('50g toppen','50g botten')
figure
plot(Untitled(:,5)), hold on
plot(Untitled(:,6)), hold on
legend('100g toppen','100g botten')
figure
plot(Untitled(:,7)), hold on
plot(Untitled(:,8)), hold on
legend('150g toppen','150g botten')
%plot(Untitled(:,9)), hold on

%legend('omgivining','50g toppen','50g botten','100g toppen','150g toppen','150g botten','referens')

xlabel('Mätningar')
ylabel('T (°C)')
            
%% 3. Set Evaluation Parameters (inspect from plot)

eval_T      = [27 40 34];    % [Tmin Tmax Tnorm]
Tmin        = eval_T(1);    
Tmax        = eval_T(2);
Tnorm       = eval_T(3);
dT          = 1e-3;

Tmax_sc     = 36;         % Region to inspect for SC
Tmin_sc     = 33;

% - Smoothing: Necessary if data is not monotoneous increasing/decreasing
isSmSample  = 0;                % Turn on/off smoothing
isSmRef     = 1;            
Tsm_sample_max_h  = 35;           % Temperature limits for smoothing (heating)
Tsm_sample_min_h  = 24;           
Tsm_sample_max_c  = [35];           % Temperature limits for smoothing (cooling)
Tsm_sample_min_c  = [24];        

T_p = 36.4;                       % Evaluation point for h_p: (T_norm - T_p)
T_solid = [25 25.5];                % Ranges for average cp
T_liquid = [38 38.5];            


%% 4. Define indices for individual cooling/heating experiment

offset = texp*3600/dt;    % Number of indices for texp

idx_1 = [14200:16600;... % Cooling
            11000:13400];    % Heating         
idx_Cycle = {idx_1};

%% 5. Create "THistory" object and store in object array

for i=1:NS
    for j=1:NC
        for k=1:2
%          for k=1 % only cooling
%          for k=2 % only heating
            
            % Cycle no. and cool/heat case
            idx_Cycle_temp = idx_Cycle{j};          % Cycle j in array "idx_Cycle" 
            idx_Cycle_temp = idx_Cycle_temp(k,:);   % Cycle j, type k in array "idx_Cycle"
            
            % Extract from data array
            T_sample   = Untitled(idx_Cycle_temp,idx_sensor{2,i});
            T_ref      = Untitled(idx_Cycle_temp,idx_sensor{3,i});
            
            % Set instance properties
            nSensor    = idx_sensor{1,i};
            nCycle     = num2str(j);
            if k==1; nType=nType_cool; else nType=nType_heat;end 
            
            % Create "ExpID_THistory" object and store in object array
            objArray(i,j,k) = THistoryExperiment(nExpID,nSample,nRef,msample,mref,mtsample,mtref,dt,... % Superclass properties
                nSensor,nCycle,nType,... % Subclass properties
                T_sample,T_ref,Tmax,Tmin,Tnorm,isSmSample,isSmRef,Tsm_sample_max_h,Tsm_sample_min_h,...
                Tsm_sample_max_c,Tsm_sample_min_c,Tmax_sc,Tmin_sc,Tamb_max,Tamb_min,dT) % Subclass properties
            % Evaluation
            objArray(i,j,k) = eval_h_diff(objArray(i,j,k));

        end
    end
end

a = ExpID_THistory_array(objArray);


%% Generate Plots
a.plot_h_array([Tmin Tmax])
