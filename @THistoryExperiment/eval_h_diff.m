function obj = eval_h_diff(obj, a_m, a_cpt, a_cpref)


if nargin==1
    a_m = 0; a_cpt = 0; a_cpref = 0;
end

% Load local data
dt = obj.dt;
Ts_0 = obj.T_sample;
Tr_0 = obj.T_ref;
ts_0 = (0:dt:dt*length(Ts_0)-1)';  % [N 1] vector
tr_0 = (0:dt:dt*length(Tr_0)-1)';  % [N 1] vector
Tmax = obj.Tmax;
Tmin = obj.Tmin;
Tnorm = obj.Tnorm;
dT = obj.dT;
Tsm_sample_max_h = obj.Tsm_sample_max_h;
Tsm_sample_min_h = obj.Tsm_sample_min_h;
Tsm_sample_max_c = obj.Tsm_sample_max_c;
Tsm_sample_min_c = obj.Tsm_sample_min_c;

isSmSample = obj.isSmSample;
isSmRef = obj.isSmRef;

mref = obj.mref;
msample = obj.msample;
mtref = obj.mtref;
mtsample = obj.mtsample;

if mean(diff(Tr_0))>0        % Heating
    %%  (1) Heating
    
    %%  (2) Find indices: Tmin -> Tmax
    
    idx_high_r = find(Tr_0>Tmax,1);
    idx_low_r  = find(Tr_0>Tmin,1);
    idx_high_s = find(Ts_0>Tmax,1);
    idx_low_s  = find(Ts_0>Tmin,1);
    idx_eval_s = (idx_low_s:idx_high_s);
    idx_eval_r = (idx_low_r:idx_high_r);
    
    Ts = Ts_0(idx_eval_s);
    Tr = Tr_0(idx_eval_r);
    ts = ts_0(idx_eval_s);
    tr = tr_0(idx_eval_r);
    
    % - Smoothing: Ref
    switch isSmRef
        case 1
            slmr = slmengine(tr,Tr,'knots',10,'increasing','on');
            Tr_smooth =  slmeval(tr,slmr,0);
            Tr = Tr_smooth;
            
            %             figure
            %             hold on; plot(tr_0/3600,Tr_0);plot(tr/3600,Tr_smooth); hold off;
            %             ylabel('Temperature in degC')
            %             xlabel('time in h')
            
        otherwise
            Tr = Tr;
    end
    
    % - Smoothing: Sample
    switch isSmSample
        case 1
            
            idx_high_sm_sample = find(Ts>=Tsm_sample_max_h,1);
            idx_low_sm_sample = find(Ts>Tsm_sample_min_h,1);
            idx_smooth = (idx_low_sm_sample:idx_high_sm_sample);
            
            slmr = slmengine(ts(idx_smooth),Ts(idx_smooth),'knots',10,'increasing','on');
            Ts_smooth =  slmeval(ts(idx_smooth),slmr,0);
            Ts(idx_smooth) = Ts_smooth;
            
            %             figure
            %             hold on; plot(ts_0/3600,Ts_0);plot(ts/3600,Ts); hold off;
            %             ylabel('Temperature in degC')
            %             xlabel('time in h')
            
        otherwise
            Ts = Ts;
    end
    
    
    % Skip minimum dT
    mask = [false; abs(diff(Ts))<dT];
    while sum(mask)>0
        i=find(mask==1,1);
        Ts(i)=[];
        ts(i)=[];
        mask = [false; abs(diff(Ts))<dT];
    end
    
    %%  (4) Interpolation
    F_r = griddedInterpolant(Tr',tr'); % transpose to row vector [1 n]
    tr_ip = F_r(Ts); % [n 1]
    
    dt_s = diff(ts);     % dt_s > 0
    dt_r_ip_abs = abs(diff(tr_ip)); % |dt_r_ip_abs| > 0
    
    %%  (5) Calculate cpeff
    [ cpref ] = obj.cp_H20_Patek_et_al_2009( Ts(1:end-1),a_cpref );  % [N 1] vector
    [ cpt ] = obj.cp_copper_Arblaster_2015( Ts(1:end-1),a_cpt );  % [N 1] vector
    
    a = (mref.*cpref + mtref.*cpt)./msample .*dt_s./dt_r_ip_abs;
    b = (mtsample/msample).*cpt;
    cp_eff= a-b;                                     % [J/(kg K)]
    cp_eff(cp_eff<0)=0;
    
    %%  (6) Calculate h
    % Note that '1:end-1' is a forward integration -> INT_i to INT_i+1 cp_eff_i dT
    h = cumtrapz_adj(Ts(1:end-1),cp_eff); % [J/(kg)]
    
    %%  (7) Shift h=0 at Tnorm
    i_n_high=find(Ts>=Tnorm,1);
    i_n_low=i_n_high-1;
    h_ip = h(i_n_low:i_n_high)'; % [1 N] vector
    T_ip = Ts(i_n_low:i_n_high)'; % [1 N] vector
    
    F_h = griddedInterpolant(T_ip,h_ip); % ([1 N], [1 N])
    h_n =  F_h(Tnorm);
    
    h = h-h_n; % [N 1] vector
    
    Ts_save = Ts(1:end-1);
    dTdt_sample_save = diff(Ts)./dt_s;
    
    
else  % Cooling
    
    %%  (1) Cooling
    
    %%  (2) Find indices: Tmax -> Tmin
    Tmax
    Tmin
    Tr_0
    idx_high_r = find(Tr_0<Tmax,1)
    idx_low_r  = find(Tr_0<Tmin,1)
    idx_high_s = find(Ts_0<Tmax,1);
    idx_low_s  = find(Ts_0<Tmin,1);
    idx_eval_s = (idx_high_s:idx_low_s);
    idx_eval_r = (idx_high_r:idx_low_r);
    
    Ts = Ts_0(idx_eval_s);
    Tr = Tr_0(idx_eval_r);
    ts = ts_0(idx_eval_s);
    tr = tr_0(idx_eval_r);
    
    % - Smoothing: Ref
    switch isSmRef
        case 1
            slmr = slmengine(tr,Tr,'knots',10,'decreasing','on');
            Tr_smooth =  slmeval(tr,slmr,0);
            Tr = Tr_smooth;
            
            %                         figure
            %                         hold on; plot(tr_0/3600,Tr_0);plot(tr/3600,Tr_smooth); hold off;
            %                         ylabel('Temperature in degC')
            %                         xlabel('time in h')
            
        otherwise
            Tr = Tr;
    end
    
    % - Smoothing: Sample
    switch isSmSample
        case 'never'
            
            idx_high_sm_sample = find(Ts<Tsm_sample_max_c,1);
            idx_low_sm_sample = find(Ts<Tsm_sample_min_c,1);
            idx_smooth = (idx_high_sm_sample:idx_low_sm_sample);
            
            slmr = slmengine(ts(idx_smooth),Ts(idx_smooth),'knots',10,'decreasing','on');
            Ts_smooth =  slmeval(ts(idx_smooth),slmr,0);
            Ts(idx_smooth) = Ts_smooth;
            
            %                         figure
            %                         hold on; plot(ts_0/3600,Ts_0);plot(ts/3600,Ts); hold off;
            %                         ylabel('Temperature in degC')
            %                         xlabel('time in h')
            
        otherwise
            Ts = Ts;
    end
    
    % Skip the sc recalesence onset to avoid apparant high thermal capacity due to sensor lag
    i_n = find(diff(Ts)>0,1);
    i_sc = [i_n; i_n+1;i_n+2];
    Ts(i_sc)=[];
    ts(i_sc)=[];
    
    % Skip minimum dT
    mask = [false; abs(diff(Ts))<dT];
    while sum(mask)>0
        i=find(mask==1,1);
        Ts(i)=[];
        ts(i)=[];
        mask = [false; abs(diff(Ts))<dT];
    end
    
    %%  (4) Interpolation
    F_r = griddedInterpolant(fliplr(Tr'),fliplr(tr')); % transpose to row vector [1 n]
    tr_ip = F_r(Ts); % [n 1]
    
    dt_s = diff(ts);     % dt_s > 0
    dt_r_ip_abs = abs(diff(tr_ip)); % |dt_r_ip_abs| > 0
    
    %%  (5) Calculate cp_eff
    [ cpref ] = obj.cp_H20_Patek_et_al_2009( Ts(1:end-1),a_cpref );  % [N 1] vector
    [ cpt ] = obj.cp_copper_Arblaster_2015( Ts(1:end-1),a_cpt );  % [N 1] vector
    
    a = (mref.*cpref + mtref.*cpt)./msample .*dt_s./dt_r_ip_abs;
    b = (mtsample/msample).*cpt;
    cp_eff= a-b;                                     % [J/(kg K)]
    cp_eff(cp_eff<0)=0;
    
    %  (6) Calculate enthalpy h
    % Note that '1:end-1' is a forward integration -> INT_i to INT_i+1 cp_eff_i dT
    h = -1*cumtrapz_adj(Ts(1:end-1),cp_eff); % [J/(kg)]; '-1' for cooling
    
    %  (7) Shift h=0 at Tnorm
    i_n_low = find(Ts<=Tnorm,1)
    Ts
    i_n_high = i_n_low-1
    h_ip = h(i_n_high:i_n_low)' % [1 N] vector
    T_ip = Ts(i_n_high:i_n_low)' % [1 N] vector
    
    F_h = griddedInterpolant(fliplr(T_ip),fliplr(h_ip));
    h_n =  F_h(Tnorm);
    
    h = h-h_n; % [N 1] vector
    
    
    Ts_save = Ts(1:end-1);
    dTdt_sample_save = diff(Ts)./dt_s;
    
    
end % Heating or Cooling

%% Store results in object
obj.h_eval = h; % [N 1] vector
obj.T_eval = Ts_save; % [N 1] vector
obj.t_eval = ts(1:end-1);
obj.cp_eff_avg = cp_eff;
obj.dTdt_sample_eval = dTdt_sample_save;

% figure;plot(Ts_save,h);
end

