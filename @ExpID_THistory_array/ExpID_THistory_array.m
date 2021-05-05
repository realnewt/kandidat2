classdef ExpID_THistory_array  
    
% This class evaluates an objArray of ExpID_THistory objects and
% calculates the mean values according to the dimensions of the objArray:
% objArray(Sensor,Cycle,1 or 2: Cooling or Heating)
    
% Summary of possible operations on object array
%       - generate vector or obj.Properties: 
%               ans = [objArray(:).Prop1] 
%          <==> ans = [objArray(1).Prop1, objArray(2).Prop1, ...]
%
%       - experimental standard deviation: 
%               std = sqrt(1/n-1 * sum[(q_i - q_mean)^2]) 
%                   Source: JCGM 100:2008 // B.2.17 

     properties (SetAccess = private, GetAccess = public)
        
         
         % Object array parameters
         NSensor;    % Number of sensors
         NCycle;     % Total number of cycles
         objArray;
         
         nObjID_array;     % 'Name of object': 'nSample_ExpID' (e.g. RT28HC_B2I)

         
         % Evaluation parameters
         
         %%%-Value
            T_p_array;
            Tnorm_array;
            T_s_range_array;
            T_l_range_array;
            Tmax_array;
            Tmin_array;
            Tamb_max_array;
            Tamb_min_array;
            
            % Standard deviations
            std_dTdt_ref_array_c=[];   % Average std for dTdt(ref) for the evaluation temperature intervall
            std_dTdt_ref_array_h=[];   % Average std for dTdt(ref) for the evaluation temperature intervall
            
            std_dTdt_sample_I_array_c=[];   % Average std for dTdt(sample, part I:liquid) for the evaluation temperature intervall
            std_dTdt_sample_II_array_c=[];   % Average std for dTdt(sample, part II:phase change) for the evaluation temperature intervall
            std_dTdt_sample_III_array_c=[];   % Average std for dTdt(sample, part III:solid) for the evaluation temperature intervall

            
            % Average for all sensors and cycles
            
                % All
                    h_p_avg_all;          % Average (cooling & heating) enthalpy value at T_p for all sensors and cycles 
                    cp_s_avg_all;         % Average (cooling & heating) solid spec. heat capacity for all sensors and cycles
                    cp_l_avg_all;         % Average (cooling & heating) liquid spec. heat capacity for all sensors and cycles
                    
                    h_p_std_all;          % Standard deviation enthalpy value at T_p for all sensors and cycles
                    cp_s_std_all;         % Standard deviation solid spec. heat capacity for all sensors and cycles            
                    cp_l_std_all;         % Standard deviation liquid spec. heat capacity for all sensors and cycles  
                    
            
                % Cooling
                    h_p_avg_all_c;          % Average (cooling) enthalpy value at T_p for all sensors and cycles                    
                    cp_s_avg_all_c;         % Average (cooling) solid spec. heat capacity for all sensors and cycles            
                    cp_l_avg_all_c;         % Average (cooling) liquid spec. heat capacity for all sensors and cycles
                    
                    h_p_std_all_c;          % Standard deviation (cooling) enthalpy value at T_p for all sensors and cycles
                    cp_s_std_all_c;         % Standard deviation (cooling) solid spec. heat capacity for all sensors and cycles            
                    cp_l_std_all_c;         % Standard deviation (cooling) liquid spec. heat capacity for all sensors and cycles                    
                    

                % Heating
                    h_p_avg_all_h;          % Average (heating) enthalpy value at T_p for all sensors and cycles                    
                    cp_s_avg_all_h;         % Average (heating) solid spec. heat capacity for all sensors and cycles            
                    cp_l_avg_all_h;         % Average (heating) liquid spec. heat capacity for all sensors and cycles
                    
                    h_p_std_all_h;          % Standard deviation (heating) enthalpy value at T_p for all sensors and cycles
                    cp_s_std_all_h;         % Standard deviation (heating) solid spec. heat capacity for all sensors and cycles            
                    cp_l_std_all_h;         % Standard deviation (heating) liquid spec. heat capacity for all sensors and cycles
            
         %%%-Vector
            % per Sensor: Average for all cycles
                % Cooling
                    h_p_avg_sensor_c;       % Average (cooling) enthalpy value at T_p for all cycles per sensor (Vector)
                    cp_s_avg_sensor_c;      % Average (cooling) solid spec. heat capacity for all cycles per sensor (Vector)
                    cp_l_avg_sensor_c;      % Average (cooling) liquid spec. heat capacity for all cycles per sensor (Vector)
                    
                    h_p_std_sensor_c;       % Standard deviation (cooling) enthalpy value at T_p for all cycles per sensor (Vector)
                    cp_s_std_sensor_c;      % Standard deviation (cooling) solid spec. heat capacity for all cycles per sensor (Vector)
                    cp_l_std_sensor_c;      % Standard deviation (cooling) liquid spec. heat capacity for all cycles per sensor (Vector)


                % Heating
                    h_p_avg_sensor_h;       % Average (heating) enthalpy value at T_p for all cycles per sensor (Vector)
                    cp_s_avg_sensor_h;      % Average (heating) solid spec. heat capacity for all cycles per sensor (Vector)
                    cp_l_avg_sensor_h;      % Average (heating) liquid spec. heat capacity for all cycles per sensor (Vector)
                    
                    h_p_std_sensor_h;       % Standard deviation (heating) enthalpy value at T_p for all cycles per sensor (Vector)
                    cp_s_std_sensor_h;      % Standard deviation (heating) solid spec. heat capacity for all cycles per sensor (Vector)
                    cp_l_std_sensor_h;      % Standard deviation (heating) liquid spec. heat capacity for all cycles per sensor (Vector)
                    
     end
     
     methods
         
         %% Class constructor for 'obj'
         function obj = ExpID_THistory_array(objArray)

             % ------- Error Checks
                 T_p_array_ = [objArray(:,:,:).T_p];                % (1 x NS*NC*2) vector
                 Tnorm_array_ = [objArray(:,:,:).Tnorm];            % (1 x NS*NC*2) vector
                 T_s_range_array_ = [objArray(:,:,:).T_s_range];    % (1 x NS*NC*2) vector
                 T_l_range_array_ = [objArray(:,:,:).T_l_range];    % (1 x NS*NC*2) vector
                 Tmax_array_ = [objArray(:,:,:).Tmax];            % (1 x NS*NC*2) vector
                 Tmin_array_ = [objArray(:,:,:).Tmin];            % (1 x NS*NC*2) vector
                 Tamb_max_array_ = [objArray(:,:,:).Tamb_max];            % (1 x NS*NC*2) vector
                 Tamb_min_array_ = [objArray(:,:,:).Tamb_min];            % (1 x NS*NC*2) vector
                 
                 % Check whether equal
                 T_p_array_ = unique(T_p_array_);                   % All values have to be equal
                 Tnorm_array_ = unique(Tnorm_array_);               % All values have to be equal
                 T_s_range_array_ = unique(T_s_range_array_);       % All values have to be equal
                 T_l_range_array_ = unique(T_l_range_array_);       % All values have to be equal
                 Tmax_array_ = unique(Tmax_array_);               % All values have to be equal
                 Tmin_array_ = unique(Tmin_array_);               % All values have to be equal
                 Tamb_max_array_ = unique(Tamb_max_array_);               % All values have to be equal
                 Tamb_min_array_ = unique(Tamb_min_array_);               % All values have to be equal
                 
                 if length(T_p_array_)>1 || length(Tnorm_array_)>1
                     error('T_p or Tnorm values in obj.Array are not all equal'); 
                 else
                     obj.T_p_array = T_p_array_;
                     obj.Tnorm_array = Tnorm_array_;
                     obj.T_s_range_array = T_s_range_array_;
                     obj.T_l_range_array = T_l_range_array_;
                     obj.Tmax_array = Tmax_array_;
                     obj.Tmin_array = Tmin_array_;
                     obj.Tamb_max_array = Tamb_max_array_;
                     obj.Tamb_min_array = Tamb_min_array_;

                    obj.nObjID_array = [objArray(1,1,1).nSample '\_' objArray(1,1,1).nExpID];
                 end
             
             % ------- 
                 [NS_, NC_, NHC_]= size(objArray);

                 obj.objArray = objArray;
                 obj.NSensor = NS_;
                 obj.NCycle = NC_;
             % -------
             
             % -------
             % Average for all sensors and cycles
                % All
                    obj.h_p_avg_all     = mean([objArray(:,:,:).h_p]);      % single value
                    obj.cp_s_avg_all	= mean([objArray(:,:,:).cp_s]);     % single value                   
                    obj.cp_l_avg_all	= mean([objArray(:,:,:).cp_l]);     % single value
                    
                    obj.h_p_std_all     = std([objArray(:,:,:).h_p]);      % single value
                    obj.cp_s_std_all	= std([objArray(:,:,:).cp_s]);     % single value                
                    obj.cp_l_std_all	= std([objArray(:,:,:).cp_l]);     % single value
                    
                % Cooling
                    obj.h_p_avg_all_c   = mean([objArray(:,:,1).h_p]);      % single value           
                    obj.cp_s_avg_all_c	= mean([objArray(:,:,1).cp_s]);     % single value                   
                    obj.cp_l_avg_all_c	= mean([objArray(:,:,1).cp_l]);     % single value
                    
                    obj.h_p_std_all_c	= std([objArray(:,:,1).h_p]);      % single value
                    obj.cp_s_std_all_c	= std([objArray(:,:,1).cp_s]);     % single value                
                    obj.cp_l_std_all_c	= std([objArray(:,:,1).cp_l]);     % single value                    
                    
                % Heating
                    obj.h_p_avg_all_h   = mean([objArray(:,:,2).h_p]);      % single value           
                    obj.cp_s_avg_all_h	= mean([objArray(:,:,2).cp_s]);     % single value                   
                    obj.cp_l_avg_all_h	= mean([objArray(:,:,2).cp_l]);     % single value
                    
                    obj.h_p_std_all_h	= std([objArray(:,:,2).h_p]);      % single value  
                    obj.cp_s_std_all_h	= std([objArray(:,:,2).cp_s]);     % single value                
                    obj.cp_l_std_all_h	= std([objArray(:,:,2).cp_l]);     % single value 
             % -------
             
            %%%-Vector
                % per Sensor: Average for all cycles
                         for i=1:NS_ % Loop through all sensors := Average over all cycles
                             % Cooling
                                 h_p_avg_sensor_c_(i)  = mean([objArray(i,:,1).h_p]);        % (1 x NS_) vector
                                 cp_s_avg_sensor_c_(i) = mean([objArray(i,:,1).cp_s]);       % (1 x NS_) vector
                                 cp_l_avg_sensor_c_(i) = mean([objArray(i,:,1).cp_l]);       % (1 x NS_) vector
                                 
                                 h_p_std_sensor_c_(i)  = std([objArray(i,:,1).h_p]);        % (1 x NS_) vector
                                 cp_s_std_sensor_c_(i) = std([objArray(i,:,1).cp_s]);       % (1 x NS_) vector
                                 cp_l_std_sensor_c_(i) = std([objArray(i,:,1).cp_l]);       % (1 x NS_) vector
                                 
                             % Heating   
                                 h_p_avg_sensor_h_(i)  = mean([objArray(i,:,2).h_p]);        % (1 x NS_) vector
                                 cp_s_avg_sensor_h_(i) = mean([objArray(i,:,2).cp_s]);       % (1 x NS_) vector
                                 cp_l_avg_sensor_h_(i) = mean([objArray(i,:,2).cp_l]);       % (1 x NS_) vector 
                                 
                                 h_p_std_sensor_h_(i)  = std([objArray(i,:,2).h_p]);        % (1 x NS_) vector
                                 cp_s_std_sensor_h_(i) = std([objArray(i,:,2).cp_s]);       % (1 x NS_) vector
                                 cp_l_std_sensor_h_(i) = std([objArray(i,:,2).cp_l]);       % (1 x NS_) vector
                         end
             
             obj.h_p_avg_sensor_c = h_p_avg_sensor_c_;
             obj.cp_s_avg_sensor_c = cp_s_avg_sensor_c_;
             obj.cp_l_avg_sensor_c = cp_l_avg_sensor_c_;
             
             obj.h_p_std_sensor_c = h_p_std_sensor_c_;
             obj.cp_s_std_sensor_c = cp_s_std_sensor_c_;
             obj.cp_l_std_sensor_c = cp_l_std_sensor_c_;

             obj.h_p_avg_sensor_h = h_p_avg_sensor_h_;
             obj.cp_s_avg_sensor_h = cp_s_avg_sensor_h_;
             obj.cp_l_avg_sensor_h = cp_l_avg_sensor_h_;
             
             obj.h_p_std_sensor_h = h_p_std_sensor_h_;
             obj.cp_s_std_sensor_h = cp_s_std_sensor_h_;
             obj.cp_l_std_sensor_h = cp_l_std_sensor_h_;
             
         end
         
         %% Evaluations
         function obj = dTdt_sample_calculate_std_c(obj,T_range)
            
             if nargin < 2
                 Tmax_=obj.Tmax_array;
                 Tmin_=obj.Tmin_array;
             else
                 Tmax_ = max(T_range);
                 Tmin_ = min(T_range);
             end
             
             objArray_ = obj.objArray;
             [NS_, NC_, NHC_]= size(objArray_);
             
             vSensor = 1:NS_;
             vCycle = 1:NC_;
             % Only cooling cases
             vCase = 1;
             
             n = 1;
             for k = vCase    % cooling or heating
                 for i = vSensor          % Sensor
                     for j = vCycle       % Cycle
                         
                         % extract index
                         i_eval_c_part_I_ = objArray_(i,j,k).i_eval_c_part_I;
                         i_eval_c_part_II_ = objArray_(i,j,k).i_eval_c_part_II;
                         T_eval_ = objArray_(i,j,k).T_eval;
                                                  
                         T_eval_c_part_I_end_(n) = min(T_eval_(i_eval_c_part_I_)); 
                         T_eval_c_part_II_start_(n) = max(T_eval_(i_eval_c_part_II_));
                         
                         n=n+1;
                     end
                 end
             end
             
             T_eval_c_part_I_end_ = max(T_eval_c_part_I_end_);
             T_eval_c_part_II_start_ = min(T_eval_c_part_II_start_);
             
             % Manual setting parts
                 % Part I: Liquid -> (earliest) Start of Recalescence
                 T_eval_c_part_I_start_ = Tmax_;
                 T_eval_c_part_I_end_ = 27;

                 % Part II: (latest) end of Recalescence -> End of Latent Heat
%                  T_eval_c_part_II_start_ = T_eval_c_part_II_start_
                 T_eval_c_part_II_start_ = 27.6;
                 T_eval_c_part_II_end_ = 24;

                 % Part III: End of Latent Heat -> Sensible Heat
                 T_eval_c_part_III_start_ = 24;
                 T_eval_c_part_III_end_ = Tmin_;
             
             % Calculating for each Part
             dT = 0.01;
             
             T_eval_I_ = T_eval_c_part_I_start_:-dT:T_eval_c_part_I_end_;
             T_eval_II_ = T_eval_c_part_II_start_:-dT:T_eval_c_part_II_end_;
             T_eval_III_ = T_eval_c_part_III_start_:-dT:T_eval_c_part_III_end_;
             
             n=1;
             for k = vCase    % cooling or heating
                 for i = vSensor          % Sensor
                     for j = vCycle       % Cycle
                         
                         % For each element in T_eval_:                             
                             % for each i,j,k case, extract raw data:
                                dTdt_s_     = objArray_(i,j,k).dTdt_sample_eval;        % [1 N] vector
                                T_s_        = (objArray_(i,j,k).T_eval)';                  % [1 N] vector
                               
                             % Find indices for I, II, III
                                % Part I
                                    idx_high_I = find(T_s_<=max(T_eval_I_),1);
                                    idx_low_I  = find(T_s_<=min(T_eval_I_),1);
                                    i_eval_I = (idx_high_I:idx_low_I);
                                    
                                    T_s_I_ = T_s_(i_eval_I);
                                    dTdt_s_I_ = dTdt_s_(i_eval_I);
                                % Part II
                                
                                    idx_high_II = find(T_s_>=max(T_eval_II_),1,'last'); % needs attention
                                    idx_low_II  = find(T_s_<=min(T_eval_II_),1);
                                    i_eval_II = (idx_high_II:idx_low_II);
                                    
                                    T_s_II_ = T_s_(i_eval_II);
                                    dTdt_s_II_ = dTdt_s_(i_eval_II);
                                    
                                    idx_del = find(dTdt_s_II_>=0);

                                    if isempty(find(diff(T_s_II_)>=0,1))==0
                                        n
                                        continue;
                                    end
%                                     figure; plot(T_s_II_,dTdt_s_II_,'.:');
                                    
                                 % Part III
                                    idx_high_III = find(T_s_<=max(T_eval_III_),1);
                                    idx_low_III  = find(T_s_<=min(T_eval_III_),1);
                                    i_eval_III = (idx_high_III:idx_low_III);
                                    
                                    T_s_III_ = T_s_(i_eval_III);
                                    dTdt_s_III_ = dTdt_s_(i_eval_III);   
                                    
                             % - Interpolate dTdt -> dTdt_eval_
                                F_interp_I = griddedInterpolant(fliplr(T_s_I_),fliplr(dTdt_s_I_));
                                dTdt_s_I_eval_ = F_interp_I(T_eval_I_)'; % [N 1] column vector
                                
                                F_interp_II = griddedInterpolant(fliplr(T_s_II_),fliplr(dTdt_s_II_));
                                dTdt_s_II_eval_=F_interp_II(T_eval_II_)'; % [N 1] column vector
                                 
                                F_interp_III = griddedInterpolant(fliplr(T_s_III_),fliplr(dTdt_s_III_));
                                dTdt_s_III_eval_=F_interp_III(T_eval_III_)'; % [N 1] column vector
                                
                             % - Store results
                                dTdt_eval_I(:,n) = dTdt_s_I_eval_;
                                dTdt_eval_II(:,n) = dTdt_s_II_eval_;
                                dTdt_eval_III(:,n) = dTdt_s_III_eval_;

                                figure
                                hold on
                                plot(T_s_I_,dTdt_s_I_,'.:')
                                plot(T_eval_I_,dTdt_s_I_eval_,'o')
                                plot(T_s_II_,dTdt_s_II_,'.:')
                                plot(T_eval_II_,dTdt_s_II_eval_,'o')
                                plot(T_s_III_,dTdt_s_III_,'.:')
                                plot(T_eval_III_,dTdt_s_III_eval_,'o')
                                hold off

                         n=n+1;
                     end
                 end
             end

                
             %      - Calculate mean dTdt
             dTdt_eval_I_mean = mean(dTdt_eval_I,2);    % mean(A,2): column vector with mean of each row
             dTdt_eval_II_mean = mean(dTdt_eval_II,2);    % mean(A,2): column vector with mean of each row
             dTdt_eval_III_mean = mean(dTdt_eval_III,2);    % mean(A,2): column vector with mean of each row

             %      - Calculate std dTdt
             dTdt_eval_I_std = std(dTdt_eval_I,0,2);    % std(A,0,2): column vector with mean of each row
             dTdt_eval_II_std = std(dTdt_eval_II,0,2);    % std(A,0,2): column vector with mean of each row
             dTdt_eval_III_std = std(dTdt_eval_III,0,2);    % std(A,0,2): column vector with mean of each row

             std_I_mean = mean(dTdt_eval_I_std)
             std_II_mean = mean(dTdt_eval_II_std)
             std_III_mean = mean(dTdt_eval_III_std)

             figure;hold on; 
             plot(T_eval_I_,dTdt_eval_I_mean);plot(T_eval_I_,dTdt_eval_I_mean+std_I_mean);plot(T_eval_I_,dTdt_eval_I_mean-std_I_mean);
             plot(T_eval_II_,dTdt_eval_II_mean);plot(T_eval_II_,dTdt_eval_II_mean+std_I_mean);plot(T_eval_II_,dTdt_eval_II_mean-std_II_mean);
             plot(T_eval_III_,dTdt_eval_III_mean);plot(T_eval_III_,dTdt_eval_III_mean+std_III_mean);plot(T_eval_III_,dTdt_eval_III_mean-std_III_mean);
             hold off

             figure;hold on; 
                plot(T_eval_I_,dTdt_eval_I_std);
                c1=refline([0 std_I_mean]);
                c1.Color = 'b';
                plot(T_eval_II_,dTdt_eval_II_std);
                c2=refline([0 std_II_mean]);
                c2.Color = 'r';
                plot(T_eval_III_,dTdt_eval_III_std);
                c3=refline([0 std_III_mean]);
                c3.Color = 'r';
                hold off
                 
            obj.std_dTdt_sample_I_array_c=std_I_mean;   % Average std for dTdt(sample, part I:liquid) for the evaluation temperature intervall
            obj.std_dTdt_sample_II_array_c=std_II_mean;    % Average std for dTdt(sample, part II:phase change) for the evaluation temperature intervall
            obj.std_dTdt_sample_III_array_c=std_III_mean;   % Average std for dTdt(sample, part III:solid) for the evaluation temperature intervall    
                
         end
         
         function obj = dTdt_ref_calculate_std(obj,T_range)
             
             % This function estimates the mean standard deviation for the
             % reference first derivative dTdt.
             %  Methodology:    - T_eval is defined and each data set dTdt_ref(i,j,k) is interpolated to each T_eval element.
             %                  - Then a std from (i,j,k) can be calculated for each T_eval element.
             %                  - The average std is calculated from all std values
             
             if nargin < 2
                 Tmax_=obj.Tmax_array;
                 Tmin_=obj.Tmin_array;
             else
                 Tmax_ = max(T_range);
                 Tmin_ = min(T_range);
             end
                 
             dT = 0.01;
             
             objArray_ = obj.objArray;
             [NS_, NC_, NHC_]= size(objArray_);
             
             vSensor = 1:NS_;
             vCycle = 1:NC_;
             vCase = 1:NHC_;

             
             % Preallocate dTdt_eval_ vector for heating and cooling case
             T_r_eval_c_ = Tmax_:-dT:Tmin_;
             T_r_eval_h_ = Tmin_:+dT:Tmax_;
             dTdt_eval_c = zeros(length(T_r_eval_c_),NS_*NC_);
             dTdt_eval_h = zeros(length(T_r_eval_h_),NS_*NC_);
             n = 1; % index for NS_*NC_ cases

             for k = vCase    % cooling or heating
                     for i = vSensor          % Sensor
                         for j = vCycle       % Cycle
                             
                             % For each element in T_eval_:
                             %      - Interpolate dTdt -> dTdt_eval_
                             
%                              i_eval_r_ = objArray_(i,j,k).i_eval_r
                             dTdt_r_ = objArray_(i,j,k).dTdt_ref(objArray_(i,j,k).i_eval_r);
                             T_r_ = objArray_(i,j,k).T_ref(objArray_(i,j,k).i_eval_r);
                             
%                              figure;plot(T_r_,dTdt_r_);
                             
                             if k==2 % Heating
                                 F_interp_r = griddedInterpolant(T_r_',dTdt_r_');
                                 dTdt_r_eval_=F_interp_r(T_r_eval_h_)'; % [N 1] column vector
                                 
                                 dTdt_eval_h(:,n) = dTdt_r_eval_;
                                 
%                                  figure;hold on; 
%                                  plot(T_r_,dTdt_r_);
%                                  plot(T_r_eval_h_,dTdt_r_eval_);hold off;
                                 
                             elseif k==1
                                 F_interp_r = griddedInterpolant(fliplr(T_r_'),fliplr(dTdt_r_'));
                                 dTdt_r_eval_=F_interp_r(T_r_eval_c_)'; % [N 1] column vector

                                 dTdt_eval_c(:,n) = dTdt_r_eval_;
                                 
                                 dTdt_eval_h(:,n) = dTdt_r_eval_;
%                                  figure;hold on; 
%                                  plot(T_r_,dTdt_r_);
%                                  plot(T_r_eval_c_,dTdt_r_eval_);hold off;

                             end

                            n = n+1;
                         end
                     end
                     n = 1;
             end
             
             %      - Calculate mean dTdt
             dTdt_eval_c_mean = mean(dTdt_eval_c,2);    % mean(A,2): column vector with mean of each row
             dTdt_eval_h_mean = mean(dTdt_eval_h,2);    % mean(A,2): column vector with mean of each row
             %      - Calculate std dTdt
             dTdt_eval_c_std = std(dTdt_eval_c,0,2);    % std(A,0,2): column vector with mean of each row
             dTdt_eval_h_std = std(dTdt_eval_h,0,2);    % std(A,0,2): column vector with mean of each row
             
             std_c_mean = mean(dTdt_eval_c_std)
             std_h_mean = mean(dTdt_eval_h_std)
             
             figure;hold on; plot(T_r_eval_c_,dTdt_eval_c_mean);plot(T_r_eval_c_,dTdt_eval_c_mean+std_c_mean);plot(T_r_eval_c_,dTdt_eval_c_mean-std_c_mean);hold off
             figure;hold on; plot(T_r_eval_h_,dTdt_eval_h_mean);plot(T_r_eval_h_,dTdt_eval_h_mean+std_h_mean);plot(T_r_eval_h_,dTdt_eval_h_mean-std_h_mean);hold off
             figure;hold on; 
                plot(T_r_eval_c_,dTdt_eval_c_std);
                c=refline([0 std_c_mean]);
                c.Color = 'b';
                plot(T_r_eval_h_,dTdt_eval_h_std);
                h=refline([0 std_h_mean]);
                h.Color = 'r';
                hold off
             
             obj.std_dTdt_ref_array_c = std_c_mean;
             obj.std_dTdt_ref_array_h = std_h_mean;
         
         end
         
         %% Create Plots
             function obj = plot_h_array(obj,T_range,vSensor,vCycle,vCase,nMonochrome)

                    objArray_ = obj.objArray; 
                    [NS_, NC_, NHC_]= size(objArray_);

                     if (1 < nargin) && (nargin < 5) % Add default cases
                         vSensor = 1:NS_;
                         vCase = 1:NHC_;
                         vCycle = 1:NC_;
                         nMonochrome = 0;
                     end

                     figure
                     xlabel('Temperature in degC')
                        dT = 1; % display x-axis in dT intervalls
%                         set(gca,'XLim',[obj.Tamb_min_array obj.Tamb_max_array])
%                         set(gca,'XTick',(obj.Tamb_min_array:dT:obj.Tamb_max_array))
                        set(gca,'XLim',T_range)
                        set(gca,'XTick',(min(T_range):dT:max(T_range)))
                     ylabel('Enthalpy in kJ/kg') % Note: h_eval * 1e-3
                     grid on
                     grid minor
                     box on
                     hold on
                         for k = vCase    % cooling or heating
                             for i = vSensor          % Sensor
                                 for j = vCycle       % Cycle
                                  

                                     % Case statements
                                     if nMonochrome == 0
                                         if i==1
%                                              nColour = 'b.';
                                             nColour = 'b';
                                         elseif i==2
%                                              nColour = 'r.';
                                             nColour = 'r';
                                         elseif i==3
%                                              nColour = 'g.';
                                             nColour = 'g';
                                         else
                                             nColour = 'k.';
                                         end
                                     elseif nMonochrome == 1
                                         nColour = 'k.';
                                     end

                                     if k==1
                                         nType = '-';   % Cooling
                                     else
                                         nType = '--';  % Heating
                                     end

                                     T_eval_ = []; h_eval_=[];

                                     T_eval_ = [objArray_(i,j,k).T_eval];
                                     h_eval_ = [objArray_(i,j,k).h_eval]*1e-3; % (J/kg)*1e-3 = (kJ/kg)

                                     nLegend = objArray_(i,j,k).nObjID; % String: 'Name of object': 'nSample_ExpID_nSensor_nCycle_nType' (e.g. RT28HC_B2I_t_1_c)
                                     
%                                      plot(T_eval_,h_eval_,[nType nColour],'DisplayName',nLegend)
                                     plot(T_eval_,h_eval_,[nType nColour])
                                 end
                             end
                         end
                    
                     
                                 % Custom legend
                                 l(1)=plot(NaN,NaN,'-b');
                                 l(2)=plot(NaN,NaN,'-r');
                                 l(3)=plot(NaN,NaN,'-g');
                                 l(4)=plot(NaN,NaN,'--b');
                                 l(5)=plot(NaN,NaN,'--r');
                                 l(6)=plot(NaN,NaN,'--g');
                                legend(l, 'Top sensor (cooling)', 'Center sensor (cooling)', 'Bottom sensor (cooling)',...
                                    'Top sensor (heating)', 'Center sensor (heating)', 'Bottom sensor (heating)');
                     
                      hold off;
%                      legend('show');
             end
             
             %% Boxplot
             function obj = boxplot_h_p_array(obj,vSensor,vCycle,vCase)
                 
                 objArray_ = obj.objArray;
                 [NS_, NC_, NHC_]= size(objArray_);
                 
                 if nargin < 3 % Add default cases
                     vSensor = 1:NS_;
                     vCase = 1:NHC_;
                     vCycle = 1:NC_;
                 end

                 h_p_ = [objArray_(vSensor,vCycle,vCase).h_p]*1e-3; % (J/kg)*1e-3 = (kJ/kg)

                 figure
                 boxplot(h_p_)
                 ylabel('Enthalpy in kJ/kg') % Note: h_eval * 1e-3
             end
             
             %% Plot dTdt for all
             function obj = plot_dTdt_array_all(obj,vSensor,vCycle,vCase,nMonochrome)
             
                 objArray_ = obj.objArray; 
                 [NS_, NC_, NHC_]= size(objArray_);
                 
                 if nargin < 5 % Add default cases
                         vSensor = 1:NS_;
                         vCase = 1:NHC_;
                         vCycle = 1:NC_;
                         nMonochrome = 0;
                         elseif nargin < 3
                         nMonochrome = 0;
                 end
                 
                 figure
                 xlabel('Temperature in degC')
                 dT = 1; % display x-axis in dT intervalls
                        set(gca,'XLim',[obj.Tamb_min_array obj.Tamb_max_array])
                        set(gca,'XTick',(obj.Tamb_min_array:dT:obj.Tamb_max_array))
                 ylabel('first time derivative dTdt in degC/s')
%                  grid on
%                  grid minor
                 box on
                 hold on
                     for i = vSensor          % Sensor
                         for j = vCycle       % Cycle
                             for k = vCase    % cooling or heating
                                 
                                 % Case statements
                                 if nMonochrome == 0
                                         nPlot = '.:';
                                 elseif nMonochrome == 1
                                     if k == 1 % cooling
                                         nPlot = '.:b';
                                     elseif k==2 % heating
                                         nPlot = '.:r';
                                     end
                                 end
                                 
                                 dTdt_sample_   = [objArray_(i,j,k).dTdt_sample];    % extract values
                                 T_sample_      = [objArray_(i,j,k).T_sample];       % extract values
                                 dTdt_ref_      = [objArray_(i,j,k).dTdt_ref];       % extract values
                                 T_ref_         = [objArray_(i,j,k).T_ref];          % extract values

                                 plot(T_sample_,dTdt_sample_,nPlot);
                                 plot(T_ref_,dTdt_ref_,nPlot);
                             end
                         end
                     end
                 T_s_range_array_ = obj.T_s_range_array;
                 T_l_range_array_ = obj.T_l_range_array;
                 % Plot Vertical lines
                 plot([min(T_s_range_array_) min(T_s_range_array_)],  ylim, 'k','LineWidth',2)   % Plot vertical line
                 plot([max(T_s_range_array_) max(T_s_range_array_)],  ylim, 'k','LineWidth',2)   % Plot vertical line 
                 plot([min(T_l_range_array_) min(T_l_range_array_)],  ylim, 'k','LineWidth',2)   % Plot vertical line 
                 plot([max(T_l_range_array_) max(T_l_range_array_)],  ylim, 'k','LineWidth',2)   % Plot vertical line 
                 hold off
             end
             
             %% Plot dTdt for sample
             function obj = plot_dTdt_array_sample(obj,vSensor,vCycle,vCase,nMonochrome)
             
                 objArray_ = obj.objArray; 
                 [NS_, NC_, NHC_]= size(objArray_);
                 
                 if nargin < 4 % Add default cases
                         vSensor = 1:NS_;
                         vCase = 1:NHC_;
                         vCycle = 1:NC_;
                         nMonochrome = 0;
                         elseif nargin < 3
                         nMonochrome = 0;
                 end
                 
                 figure
                 xlabel('Temperature in degC')
                 dT = 1; % display x-axis in dT intervalls
                        set(gca,'XLim',[obj.Tamb_min_array obj.Tamb_max_array])
                        set(gca,'XTick',(obj.Tamb_min_array:dT:obj.Tamb_max_array))
                 ylabel('first time derivative dTdt in degC/s')
                 grid on
                 grid minor
                 box on
                 hold on
                     for k = vCase    % cooling or heating
                         for i = vSensor          % Sensor
                             for j = vCycle       % Cycle
                             
                                 % Case statements
                                 if nMonochrome == 0
                                         nPlot = '.:';
                                 elseif nMonochrome == 1
                                     if k == 1 % cooling
                                         nPlot = '.:b';
                                     elseif k==2 % heating
                                         nPlot = '.:r';
                                     end
                                 end
                                 
                                 dTdt_sample_   = [objArray_(i,j,k).dTdt_sample];    % extract values
                                 T_sample_      = [objArray_(i,j,k).T_sample];       % extract values
                                 
                                 nLegend = objArray_(i,j,k).nObjID; % String: 'Name of object': 'nSample_ExpID_nSensor_nCycle_nType' (e.g. RT28HC_B2I_t_1_c)
%                                  plot(T_sample_(1:end-1),dTdt_sample_,nPlot,'DisplayName',nLegend);
                                 
                                 % Custom legend
                                    plot(T_sample_(1:end-1),dTdt_sample_,nPlot)
                                    l(1)=plot(NaN,NaN,'.:b');
                                    l(2)=plot(NaN,NaN,'.:r');
                                    legend(l, 'Cooling', 'Heating')
                                 
                                 %% Plot Evaluated
%                                  dTdt_sample_eval_   = [objArray_(i,j,k).dTdt_sample_eval];    % extract values
%                                  T_eval_             = [objArray_(i,j,k).T_eval];       % extract values
%                                  
%                                  nLegend_eval = [nLegend '\_eval'];
%                                  plot(T_eval_,dTdt_sample_eval_,'.:g','DisplayName',nLegend_eval);

                             end
                         end
                     end
                 T_s_range_array_ = obj.T_s_range_array;
                 T_l_range_array_ = obj.T_l_range_array;
                 % Plot Vertical lines
                 l1=plot([min(T_s_range_array_) min(T_s_range_array_)],  ylim, 'k','LineWidth',2);   % Plot vertical line
                 l2=plot([max(T_s_range_array_) max(T_s_range_array_)],  ylim, 'k','LineWidth',2);   % Plot vertical line
                 l3=plot([min(T_l_range_array_) min(T_l_range_array_)],  ylim, 'k','LineWidth',2);   % Plot vertical line
                 l4=plot([max(T_l_range_array_) max(T_l_range_array_)],  ylim, 'k','LineWidth',2);   % Plot vertical line
                 
                 set(get(get(l1,'Annotation'),'LegendInformation'),...
                     'IconDisplayStyle','off'); % Exclude line from legend
                 set(get(get(l2,'Annotation'),'LegendInformation'),...
                     'IconDisplayStyle','off'); % Exclude line from legend
                 set(get(get(l3,'Annotation'),'LegendInformation'),...
                     'IconDisplayStyle','off'); % Exclude line from legend
                 set(get(get(l4,'Annotation'),'LegendInformation'),...
                     'IconDisplayStyle','off'); % Exclude line from legend
                     
                 hold off
                 legend('show')
             end
             %% Plot dTdt for ref
             function obj = plot_dTdt_array_ref(obj,vSensor,vCycle,vCase,nMonochrome)
             
                 objArray_ = obj.objArray; 
                 [NS_, NC_, NHC_]= size(objArray_);
                 
                 if nargin < 4 % Add default cases
                         vSensor = 1:NS_;
                         vCase = 1:NHC_;
                         vCycle = 1:NC_;
                         nMonochrome = 0;
                 elseif nargin < 3
                         nMonochrome = 0;
                 end
                 
                 figure
                 xlabel('Temperature in degC')
                 dT = 1; % display x-axis in dT intervalls
                        set(gca,'XLim',[obj.Tamb_min_array obj.Tamb_max_array])
                        set(gca,'XTick',(obj.Tamb_min_array:dT:obj.Tamb_max_array))
                 ylabel('first time derivative dTdt in degC/s')
                 grid on
                 grid minor
                 box on
                 hold on
                     for k = vCase    % cooling or heating
                         for i = vSensor          % Sensor
                             for j = vCycle       % Cycle
                             
%                                  % Case statements
%                                  if nMonochrome == 1
%                                          nPlot = '.:k';
%                                  elseif nMonochrome == 0
%                                      if k == 1 % cooling
%                                          nPlot = '.:b';
% %                                          nPlot = '-b';
%                                      elseif k==2 % heating
%                                          nPlot = '.:r';
% %                                          nPlot = '-r';
%                                      end
%                                  end

                                     % Case statements
%                                      if nMonochrome == 0
%                                          if i==1
% %                                              nColour = 'b.';
%                                              nColour = 'b';
%                                          elseif i==2
% %                                              nColour = 'r.';
%                                              nColour = 'r';
%                                          elseif i==3
% %                                              nColour = 'g.';
%                                              nColour = 'g';
%                                          else
%                                              nColour = 'k.';
%                                          end
%                                      elseif nMonochrome == 1
%                                          nColour = 'k.';
%                                      end
                                     
                                     % Case statements
                                     if nMonochrome == 0
                                             nPlot = '.:';
                                     elseif nMonochrome == 1
                                         if k == 1 % cooling
                                             nPlot = '.:b';
                                         elseif k==2 % heating
                                             nPlot = '.:r';
                                         end
                                     end

%                                      if k==1
%                                          nType = '-';   % Cooling
%                                      else
%                                          nType = '--';  % Heating
%                                      end
%                                      nPlot = [nType nColour];
                                 
                                 
                                 nLegend = objArray_(i,j,k).nObjID; % String: 'Name of object': 'nSample_ExpID_nSensor_nCycle_nType' (e.g. RT28HC_B2I_t_1_c)
     
                                 dTdt_ref_      = [objArray_(i,j,k).dTdt_ref];       % extract values
                                 T_ref_         = [objArray_(i,j,k).T_ref];          % extract values
%                                  plot(T_ref_,dTdt_ref_,nPlot,'DisplayName',nLegend);
                                 
                                 
%                                  % Custom Legend
                                 hold on
                                 plot(T_ref_(1:end-1),dTdt_ref_,nPlot)
%                                  l(1)=plot(NaN,NaN,'-b');
%                                  l(2)=plot(NaN,NaN,'-r');
%                                  l(3)=plot(NaN,NaN,'-g');
%                                  l(4)=plot(NaN,NaN,'--b');
%                                  l(5)=plot(NaN,NaN,'--r');
%                                  l(6)=plot(NaN,NaN,'--g');
%                                 legend(l, 'Top sensor (cooling)', 'Center sensor (cooling)', 'Bottom sensor (cooling)',...
%                                     'Top sensor (heating)', 'Center sensor (heating)', 'Bottom sensor (heating)');
                                 
                                   l(1)=plot(NaN,NaN,'.:b');
                                   l(2)=plot(NaN,NaN,'.:r');
                                   legend(l, 'Cooling', 'Heating')
                                 
                                 %% Plot Evaluated
%                                  dTdt_ref_eval_   = [objArray_(i,j,k).dTdt_ref_eval];    % extract values
%                                  T_eval_      = [objArray_(i,j,k).T_eval];       % extract values
%                                  
%                                  plot(T_eval_,dTdt_ref_eval_,'-g');
                                 
                             end
                         end
                     end
                     T_s_range_array_ = obj.T_s_range_array;
                     T_l_range_array_ = obj.T_l_range_array;
                     % Plot Vertical lines
                     l1=plot([min(T_s_range_array_) min(T_s_range_array_)],  ylim, 'k','LineWidth',2);   % Plot vertical line
                     l2=plot([max(T_s_range_array_) max(T_s_range_array_)],  ylim, 'k','LineWidth',2);   % Plot vertical line
                     l3=plot([min(T_l_range_array_) min(T_l_range_array_)],  ylim, 'k','LineWidth',2);   % Plot vertical line
                     l4=plot([max(T_l_range_array_) max(T_l_range_array_)],  ylim, 'k','LineWidth',2);   % Plot vertical line
                     
                     
                     set(get(get(l1,'Annotation'),'LegendInformation'),...
                         'IconDisplayStyle','off'); % Exclude line from legend
                     set(get(get(l2,'Annotation'),'LegendInformation'),...
                         'IconDisplayStyle','off'); % Exclude line from legend
                     set(get(get(l3,'Annotation'),'LegendInformation'),...
                         'IconDisplayStyle','off'); % Exclude line from legend
                     set(get(get(l4,'Annotation'),'LegendInformation'),...
                         'IconDisplayStyle','off'); % Exclude line from legend
                 
                 hold off
                 legend('show')
             end
             
              %% Plot dTdt for ref
             function obj = plot_Tt_array(obj,iSensor,iCycle,iCase)
                 
                 objArray_ = obj.objArray; 
                 
                 figure
                 ylabel('Temperature in degC')
                 xlabel('time in s')
%                  grid on
%                  grid minor
                 box on
                 hold on
                 
                  T_sample_      = [objArray_(iSensor,iCycle,iCase).T_sample];       % extract values
                  T_ref_         = [objArray_(iSensor,iCycle,iCase).T_ref];          % extract values
                  dt_            = [objArray_(iSensor,iCycle,iCase).dt];          % extract values
                  
                  nLegend = objArray_(iSensor,iCycle,iCase).nObjID; % String: 'Name of object': 'nSample_ExpID_nSensor_nCycle_nType' (e.g. RT28HC_B2I_t_1_c)
                    nLegend_sample = [nLegend '\_PCM'];
                    nLegend_ref = [nLegend '\_ref'];
                  plot(0:dt_:length(T_sample_)*dt_-1,T_sample_,'.:','DisplayName',nLegend_sample)
                  plot(0:dt_:length(T_sample_)*dt_-1,T_ref_,'.:','DisplayName',nLegend_ref)
                  
                  
                  hold off
                  legend('show')
                  
             end
             
     end
         
         
         
         
    
end