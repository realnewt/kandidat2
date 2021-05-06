classdef THistoryExperiment
    
    
    properties (SetAccess = private, GetAccess = public)
        
        % Strings
        nExpID;     % 'Experimental Setup'
        nSample;    % 'PCM name'
        nRef;       % 'Reference name'
        
        nSensor;    % 'Sensor position'
        nCycle;     % 'Number of Cycle'
        nType;      % 'Heating' or 'Cooling'
        
        nObjID;     % 'Name of object': 'nSample_ExpID_nSensor_nCycle_nType' (e.g. RT28HC_B2I_t_1_c)
        
        % Experimental Parameters
        dt;         % Data recording rate (s)
        
        msample;    % PCM sample mass (kg)
        mref;       % Reference sample mass (kg)
        mtsample;   % PCM sample holder mass (kg)
        mtref;      % Reference sample holder mass (kg)
        
        Tamb_max;   % Highest ambient temperature (degC)
        Tamb_min;   % Lowest ambient temperature (degC)
        
        % Vectors
        T_sample;       % Experimental Data: PCM (degC)
        T_ref;          % Experimental Data: Reference (degC)
        
        i_eval_s;       % Index vector for T_sample: Tmin < T_sample(i_eval_s) < Tmax
        i_eval_r;       % Index vector for T_ref: Tmin < T_ref(i_eval_r) < Tmax
        
        i_eval_c_part_I=[];         % Index vector for Part I of the cooling curve (until recalescence, for cooling dTdt std estimate)
        i_eval_c_part_II=[];        % Index vector for Part II of the cooling curve (after recalescence, for cooling dTdt std estimate)
        
        dTdt_sample;    % first time derivative from experimental data: PCM (K/s)
        dTdt_ref;       % first time derivative from experimental data: Reference (K/s)
        dTdt_sample_eval;    % first time derivative from experimental data: PCM (K/s)
        dTdt_ref_eval;       % first time derivative from experimental data: Reference (K/s)
        
        % Evaluation Parameters
        dT = [];
        
        dT_adj = 1e-2;  % Adjusted min. temp. scanning step for enthalpy calculation (degC) (to avoid errors during supercooling)
        Tmax_sc;        % Max. temp. for supercooling region (degC) (used for dT_adj)
        Tmin_sc;        % Min. temp. for supercooling region (degC) (used for dT_adj)
        
        Tmax;           % Max. temp. for enthalpy calculation (degC)
        Tmin;           % Min. temp. for enthalpy calculation (degC)
        Tnorm;          % Norm. temperature for enthalpy calculation (degC)
        
        isSmSample;         % Turn on/off smoothing
        isSmRef;            % Turn on/off smoothing
        
        Tsm_sample_max_h;    % Max. temp. for applying smoothing (degC) (Avoids over smoothing in sensible region)
        Tsm_sample_min_h;    % Min. temp. for applying smoothing (degC) (Avoids over smoothing in sensible region)
        Tsm_sample_max_c;    % Max. temp. for applying smoothing (degC) (Avoids over smoothing in sensible region)
        Tsm_sample_min_c;    % Min. temp. for applying smoothing (degC) (Avoids over smoothing in sensible region)
        
        % Evaluation Results
        T_eval=[];      % Vector: Evaluated temperature (degC)
        t_eval=[];
        h_eval=[];      % Vector: Evaluated specific enthalpy assigned to T_eval and T_norm (J/kg)
        cp_eff_avg=[];  % Vector: Evaluated effective specific heat capacity (J/(kg*K))
        
        T_p=[];        % Point Evaluation
        h_p=[];        % Point Evaluation
        
        T_s_range=[];  % Solid region temperature intervall [T_s1 Ts2] (degC)
        T_l_range=[];  % Liquid region temperature intervall [T_s1 Ts2] (degC)
        cp_s=[];       % Averaged specific heat capacity over temperature intervall (J/(kg*K)
        cp_l=[];       % Average specific heat capacity over temperature intervall (J/(kg*K)
        
        % Cell arrays for MC trials results
        T_eval_MC;
        h_eval_MC;
        h_p_MC;        % Point Evaluation
        
    end
    
    methods
        
        %% Class constructor for 'obj'
        function obj =THistoryExperiment(nExpID_,nSample_,nRef_,msample_,mref_,mtsample_,mtref_,dt_,... % Superclass properties
                nSensor_,nCycle_,nType_,... % Subclass properties
                T_sample_,T_ref_,Tmax_,Tmin_,Tnorm_,isSmSample_,isSmRef_,Tsm_sample_max_h_,Tsm_sample_min_h_,...
                Tsm_sample_max_c_,Tsm_sample_min_c_,Tmax_sc_,Tmin_sc_,Tamb_max_,Tamb_min_,dT_) % Subclass properties
            
            if nargin > 0
                
                obj.nSample = nSample_;
                obj.nRef = nRef_;
                obj.nExpID = nExpID_;
                
                obj.nSensor = nSensor_;
                obj.nCycle = nCycle_;
                obj.nType = nType_;
                
                obj.nObjID = [obj.nSample '\_' obj.nExpID '\_' obj.nSensor '\_' obj.nCycle '\_' obj.nType];
                
                obj.dt = dt_;
                obj.msample = msample_;
                obj.mref = mref_;
                obj.mtsample = mtsample_;
                obj.mtref = mtref_;
                
                obj.Tamb_max = Tamb_max_;
                obj.Tamb_min = Tamb_min_;
                
                obj.T_sample = T_sample_;
                obj.T_ref = T_ref_;
                
                obj.Tmax = Tmax_;
                obj.Tmin = Tmin_;
                obj.Tnorm = Tnorm_;
                
                obj.isSmSample = isSmSample_;
                obj.isSmRef = isSmRef_;
                obj.Tsm_sample_max_h = Tsm_sample_max_h_;
                obj.Tsm_sample_min_h = Tsm_sample_min_h_;
                obj.Tsm_sample_max_c = Tsm_sample_max_c_;
                obj.Tsm_sample_min_c = Tsm_sample_min_c_;
                obj.Tmax_sc = Tmax_sc_;
                obj.Tmin_sc = Tmin_sc_;
                
                obj.dT = dT_;
                
                % Calculate first time derivative
                % Determine heating or cooling case
                if mean(diff(T_ref_))>0        % Heating
                    type = 1;
                elseif mean(diff(T_ref_))<0    % Cooling
                    type = 2;
                end
                
                % Determine evaluation indices
                if type == 1 % Heating
                    idx_high_r = find(T_ref_>=Tmax_,1);
                    idx_low_r  = find(T_ref_>=Tmin_,1);
                    idx_high_s = find(T_sample_>=Tmax_,1);
                    idx_low_s  = find(T_sample_>=Tmin_,1);
                    i_eval_s_ = (idx_low_s:idx_high_s);
                    i_eval_r_ = (idx_low_r:idx_high_r);
                else % cooling
                    idx_high_r = find(T_ref_<=Tmax_,1);
                    idx_low_r  = find(T_ref_<=Tmin_,1);
                    idx_high_s = find(T_sample_<=Tmax_,1);
                    idx_low_s  = find(T_sample_<=Tmin_,1);
                    i_eval_s_ = (idx_high_s:idx_low_s);
                    i_eval_r_ = (idx_high_r:idx_low_r);
                end
                
                obj.i_eval_s = i_eval_s_;
                obj.i_eval_r = i_eval_r_;
                
                % First order forward differentiation
                obj.dTdt_sample = diff(T_sample_)/dt_;
                obj.dTdt_ref = diff(T_ref_)/dt_;
            end
        end
        
    end
    
    methods(Static)% Object independent functions
        %% Evaluation: Calculate reference specific heat capacity
        function [ cp_H2O ] = cp_H20_Patek_et_al_2009(t, a_cpref )
            
            % t in degC
            % T in K
            % cp in J/kg/K
            % expanded relative uncertainty k=2: +/- 0.1%
            % relative standard uncertainty = 0.05%
            % Information available: standard uncertainty -> Gaussian
            
            %             if nargin<2
            %                 a_cpref=0;  % relative standard uncertainty
            %             end
            %             a_cpref=0.05e-2;
            
            
            T = t+273.15;
            
            tau = T/10;
            alpha = 10./(593-T);
            beta = 10./(T-232);
            
            R = 461.51805;
            
            c3 = -8.983025854;
            
            n_i = [4 5 7];
            a_i = [-1.661470539e5 2.708781640e6 -1.557191544e8];
            b_i = [-8.237426256e-1 1.908956353 -2.017597384 8.546361348e-1];
            m_i = [2 3 4 5];
            
            cp_H2O_0 = -R * (c3 + tau.*( ...
                n_i(1)*(n_i(1)+1)*a_i(1)*alpha.^(n_i(1)+2)...
                + n_i(2)*(n_i(2)+1)*a_i(2)*alpha.^(n_i(2)+2)...
                + n_i(3)*(n_i(3)+1)*a_i(3)*alpha.^(n_i(3)+2) )...
                + tau.*(...
                m_i(1)*(m_i(1)+1)*b_i(1)*beta.^(m_i(1)+2)...
                + m_i(2)*(m_i(2)+1)*b_i(2)*beta.^(m_i(2)+2)...
                + m_i(3)*(m_i(3)+1)*b_i(3)*beta.^(m_i(3)+2)...
                + m_i(4)*(m_i(4)+1)*b_i(4)*beta.^(m_i(4)+2))...
                );
            
            %             cp_H2O_0 = -R * (c3 + tau.*( ...
            %                             4*(4+1)*-1.661470539e5*alpha.^(4+2)...
            %                             + 5*(5+1)*2.708781640e6*alpha.^(5+2)...
            %                             + 7*(7+1)*-1.557191544e8*alpha.^(7+2) )...
            %                     + tau.*(...
            %                             2*(2+1)*-8.237426256e-1*beta.^(2+2)...
            %                             + 3*(3+1)*1.908956353*beta.^(3+2)...
            %                             + 4*(4+1)*-2.017597384*beta.^(4+2)...
            %                             + 5*(5+1)*8.546361348e-1*beta.^(5+2))...
            %                 );
            
            
            % Random number drawn for all T (default)
            cp_H2O = cp_H2O_0 + cp_H2O_0*a_cpref; % value + value*rel.uncertainty
            
            % Random number drawn for each T
            %             cp_H2O = cp_H2O_0 + cp_H2O_0.*a_cpref.*randn(size(cp_H2O_0));
            
            
        end
        
        %% Evaluation: Calculate copper tube specific heat capacity
        function [ cp_Cu ] = cp_copper_Arblaster_2015(t, a_cpt )
            
            % t in degC
            % T in K
            % cp in J/kg/K
            % relative standard uncertainty = 0.1%
            % Information available: standard uncertainty -> Gaussian
            
            %             if nargin<2
            %                 a_cpt=0;  % relative standard uncertainty
            %             end
            
            %             a_cpt = 0.1e-2;
            
            Mw_Cu_0 = 63.546e-3; % +/- 0.003e-3 [kg/mol]
            
            T = t+273.15;
            
            cp_Cu_0 = zeros(size(T));   % Create zero vector
            
            % -73.15 C <= t < 25 C // 200 K <= T < 298.15 K
            T1 = T(T>=200 & T<298.15);
            cp_Cu_0(T>=200 & T<298.15) = 6.33481...
                + 0.162424*T1...
                - (5.78862e-4)*T1.^2 ...
                + (9.95052e-7)*T1.^3 ...
                - (6.62868e-10)*T1.^4; % J/(molK)
            
            % 25 C <= t < 1084.62 C // 298.15 K <= T < 1357.77 K
            T2 = T(T>=298.15 & T<1357.77);
            cp_Cu_0(T>=298.15 & T<1357.77) = 23.55055...
                + (6.89498e-3)*T2...
                - (2.95229e-6)*T2.^2 ...
                + (1.78088e-9)*T2.^3 ...
                - 84616.4./(T2.^2); % J/(molK)
            
            % J/(molK) -> J/(kgK)
            cp_Cu_0 = cp_Cu_0/Mw_Cu_0;
            
            % Random number drawn for all T (default)
            cp_Cu = cp_Cu_0 + cp_Cu_0*a_cpt; % value + value*rel.uncertainty
            
            % Random number drawn for each T
            %             cp_Cu = cp_Cu_0 + cp_Cu_0.*a_cpt.*randn(size(cp_Cu_0));
            cp_Cu=900;
        end
        
        %% Evaluation: Calculate copper tube specific heat capacity
        function [ cp_Cu ] = cp_copper_Sabbah_et_al_1999(t, a_cpt )
            
            % t in degC
            % T in K
            % cp in J/kg/K
            % expanded relative uncertainty k=2: +/- 0.3% for 50-300K
            % relative standard uncertainty = 0.15%
            % Information available: standard uncertainty -> Gaussian
            
            if nargin<2
                a_cpt=0;  % relative standard uncertainty
            end
            
            %             n = size(t);
            
            T = t+273.15;
            
            A_i = [-0.1285753818e1,...      % A0
                0.3098967121,...        % A1
                -0.2924985792e-1,...    % A2
                0.1418586260e-2,...     % A3
                -0.3370489513e-4,...    % A4
                0.4856675621e-6,...     % A5
                -0.4646773402e-8,...    % A6
                0.3070527023e-10,...    % A7
                -0.1419198886e-12,...   % A8
                0.455751904e-15,...     % A9
                -0.9894731263e-18,...   % A10
                0.1370529662e-20,...    % A11
                -0.1074497377e-23,...   % A12
                0.3517161374e-27];      % A13
            
            i = [0 1 2 3 4 5 6 7 8 9 10 11 12 13];
            
            %             cp_Cu_0 = sum(A_i.*T.^i)  / (63.546e-3);   % (J/(mol*K)) / (kg/mol) = J/(kg*K)
            cp_Cu_0 = 1/(63.546e-3) * (...
                A_i(1)*T.^i(1)...
                +A_i(2)*T.^i(2)...
                +A_i(3)*T.^i(3)...
                +A_i(4)*T.^i(4)...
                +A_i(5)*T.^i(5)...
                +A_i(6)*T.^i(6)...
                +A_i(7)*T.^i(7)...
                +A_i(8)*T.^i(8)...
                +A_i(9)*T.^i(9)...
                +A_i(10)*T.^i(10)...
                +A_i(11)*T.^i(11)...
                +A_i(12)*T.^i(12)...
                +A_i(13)*T.^i(13)...
                +A_i(14)*T.^i(14));
            
            cp_Cu = normrnd(cp_Cu_0,cp_Cu_0.*a_cpt);
        end
        
        z = cumtrapz_adj(x,y,dim)
        
    end
    
end
