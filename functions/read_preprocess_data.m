function dataStruct = read_preprocess_data(filesPath,dataFile,POS_OP,POS_OP_TrialsIndx,lpFilter,subjID,varargin)
%% This function reads all trials with fixed ankle position within an FLB file of an experiment
options={{'decimation' 10 'decimation factor'} ...
         {'norm_method' 'max' 'EMG normalization method. Possible options are: self, max, mvc'} ...
         {'plot_raw_data' 1 'flag to whether plot raw data'} ...
         {'plot_decimated_data' 1 'flag to whether plot decimated data'} ...
         {'pca_thresh' 95 'threshold in percentage on the input variance expalined by PCs to keep PCs'} ... % {'pca_thresh' 0.1 'threshold on singular values of PCA to keep features'} ...
         {'data_file_format' 'flb' 'file formats for data: either flb or csv'} ...
         {'lp_filter_tf' NaN 'EMG envelop low-pass filter in transfer function format'} ...
         {'transient_time' 0 'transient time to be discarded due to EMG envelope filter transient'} ...
         {'bp_filter' [] 'EMG envelop band-pass filter'}
         };

if arg_parse(options,varargin)
    return
end

%-- We enforce the sampling rate of data acquisition to 0.001s since: 
%-- (1) the csv file has no sampling time, and 
%-- (2) the recorded flb file has a sampling time that is not precisely 0.001s - it is 1.000000047497451e-03 seconds 
Ts = 0.001;

TRIALS = POS_OP_TrialsIndx;

nOP = length(POS_OP);
[POS_OP_SORTED, INDX] = sort(POS_OP);

%-- set the bpFilter for EMG envelope
bpFilter = bp_filter;

%-- Signal indices in the data
switch subjID
    case {'IES01'}
        POS = 1;
        TQ = 2;
        GM = 4;   %-- Gstrocnemius Medial EMG
        GL = 6;   %-- Gastrcnemius Lateral EMG
        Sol = 5;  %-- Soleus
        TA = 3;   %-- Tibialis Anterior
    otherwise
        POS = 1;
        TQ = 2;
        GM = 3;   %-- Gstrocnemius Medial EMG
        GL = 4;   %-- Gastrcnemius Lateral EMG
        Sol = 5;  %-- Soleus
        TA = 6;   %-- Tibialis Anterior
end

switch data_file_format
    case 'flb'
        data_all_trials = flb2mat([filesPath,dataFile],'read_all');
    case 'csv'
        switch subjID
            case 'ES04'
                t1 = 40;    
                t2 = 127;   
            case 'ES05'
                t1 = 1;    
                t2 = 89;   
        end
        startRow = 11;
        startCol = 0;
        data_all_trials = csv2flb(filesPath,dataFile,startRow,startCol,t1,t2,Ts);
end

dataStruct = cell(1,nOP);

for i = 1:nOP
    %-- Assign the position
    dataStruct{1,i}.POS_OP = POS_OP_SORTED(i);
    
    %-- Assign the low-pass filter used in IRF and TF format
    dataStruct{1,i}.emgFilter.irf = lpFilter;
    dataStruct{1,i}.emgFilter.tf = lp_filter_tf;

    j = INDX(i);
    trial_data = data_all_trials{1,TRIALS(j)};
    trial_data.domainIncr = Ts;

    %-- Read in signals of the trial
    pos = trial_data(:,POS); 
    tq = trial_data(:,TQ);
    gm = demean(trial_data(:,GM));  %-- demean is required for EMG signals acquired by the wireless EMG sensors.
    gl = demean(trial_data(:,GL));
    sol = demean(trial_data(:,Sol));
    ta = demean(trial_data(:,TA));
    
    %-- Plot raw data
    if plot_raw_data == 1
        figure;
        subplot(6,1,1)
        plot(pos); ylabel('Pos (rad)'); xlabel(''); title('Raw Data'); % title('Unperturbed Trial, Quasi-Stationary PF Position, TV Torque'); 
        subplot(6,1,2)
        plot(tq); ylabel('TQ (Nm)'); xlabel(''); title('')
        subplot(6,1,3)
        plot(abs(gm)); ylabel('GM EMG (mV)'); xlabel(''); title('')
        subplot(6,1,4)
        plot(abs(gl)); ylabel('GL EMG (mV)'); xlabel(''); title('')
        subplot(6,1,5)
        plot(abs(sol)); ylabel('Sol EMG (mV)'); xlabel(''); title('')
        subplot(6,1,6)
        plot(abs(ta)); ylabel('TA EMG (mV)'); xlabel('Time (s)'); title('')
        xAxisPanZoom
    end

    %-- Assigning output structure
    dataStruct{1,i}.raw.pos = pos;
    dataStruct{1,i}.raw.tq = tq;
    dataStruct{1,i}.raw.gm = gm;
    dataStruct{1,i}.raw.gl = gl;
    dataStruct{1,i}.raw.sol = sol;
    dataStruct{1,i}.raw.ta = ta;

    %-- Decimate the data
    pos_d = decimate(pos, decimation);
    tq_d = decimate(tq, decimation);
    gm_d = decimate(gm, decimation);
    gl_d = decimate(gl, decimation);
    sol_d = decimate(sol, decimation);
    ta_d = decimate(ta, decimation);

    if plot_decimated_data == 1
        figure;
        subplot(6,1,1)
        plot(pos_d); ylabel('Pos (rad)'); xlabel(''); title('Decimated Data'); % title('Unperturbed Trial, Quasi-Stationary PF Position, TV Torque'); 
        subplot(6,1,2)
        plot(tq_d); ylabel('TQ (Nm)'); xlabel(''); title('')
        subplot(6,1,3)
        plot(abs(gm_d)); ylabel('GM EMG (mV)'); xlabel(''); title('')
        subplot(6,1,4)
        plot(abs(gl_d)); ylabel('GL EMG (mV)'); xlabel(''); title('')
        subplot(6,1,5)
        plot(abs(sol_d)); ylabel('Sol EMG (mV)'); xlabel(''); title('')
        subplot(6,1,6)
        plot(abs(ta_d)); ylabel('TA EMG (mV)'); xlabel('Time (s)'); title('')
    end

    %-- Assigning output structure
    dataStruct{1,i}.decimated.pos = pos_d;
    dataStruct{1,i}.decimated.tq = tq_d;
    dataStruct{1,i}.decimated.gm = gm_d;
    dataStruct{1,i}.decimated.gl = gl_d;
    dataStruct{1,i}.decimated.sol = sol_d;
    dataStruct{1,i}.decimated.ta = ta_d;

    dataStruct{1,i}.emg_raw = cat(2, dataStruct{1,i}.raw.gm, dataStruct{1,i}.raw.gl, dataStruct{1,i}.raw.sol);
    dataStruct{1,i}.emg_abs = cat(2, abs(dataStruct{1,i}.raw.gm), abs(dataStruct{1,i}.raw.gl), abs(dataStruct{1,i}.raw.sol));
    dataStruct{1,i}.emg_env = emg_envelope(dataStruct{1,i}.emg_raw,bpFilter,lpFilter);

end

%% Normalizing the EMG envelopes
%-- First, let's calculate the maximum of EMG envelopes at various positions
nEMG = 3;
max_emg_env = zeros(nEMG,nOP);
max_emg_abs = zeros(nEMG,nOP);

for i = 1:nOP
    max_emg_env(:,i) = max(dataStruct{1,i}.emg_env).dataSet';
    max_emg_abs(:,i) = max(dataStruct{1,i}.emg_abs).dataSet';
end

global_max_emg_env = zeros(nEMG,1);
global_max_emg_abs = zeros(nEMG,1);
for j = 1:nEMG
    global_max_emg_env(j,1) = max(max_emg_env(j,:));
    global_max_emg_abs(j,1) = max(max_emg_abs(j,:));
end

for i = 1:nOP
    switch norm_method
        case 'self'
            dataStruct{1,i}.emg_abs_norm = dataStruct{1,i}.emg_abs ./ max_emg_abs(:,i)';
            dataStruct{1,i}.emg_env_norm = dataStruct{1,i}.emg_env ./ max_emg_env(:,i)';
            dataStruct{1,i}.emg_env_norm_factor = 1 ./ max_emg_env(:,i)';
        
        case 'max'
            dataStruct{1,i}.emg_abs_norm = dataStruct{1,i}.emg_abs ./ global_max_emg_abs';
            dataStruct{1,i}.emg_env_norm = dataStruct{1,i}.emg_env ./ global_max_emg_env';
            dataStruct{1,i}.emg_env_norm_factor = 1 ./ global_max_emg_env';
        
        case 'mvc'
            disp('This normalization method is not yet implemented.')
            return;
        
        otherwise
            disp('The normalization method is not defined.')
            return;
    end
end

%% Principal component analysis on normalized EMG envelopes
% Ts = pos.domainIncr;
for i = 1:nOP
    [coeff,~,latent,~,explainedVar] = pca(dataStruct{1,i}.emg_env_norm.dataSet,'Algorithm','svd');
    %== The intial method I had used was based on normalizing latent (or singular values w.r.t. to the largest one)
    % latent_n = latent/latent(1);
    % nPC = sum(latent_n > pca_thresh);
    
    %== The new method (June 2024) is based on using explained variance (i.e., how much of the input variance each PC explains) 
    %++ This can either be calculated from the latent as follows: totalVar = sum(latent); explainedVar = 100 * latent / totalVar;
    %++ Or simply be obtained as output argument of pca -> I chose this.
    
    %++ Now, we set nPC based on thresholding variance explained
    cumulativeExplainedVar = cumsum(explainedVar);
    nPC = find(cumulativeExplainedVar >= pca_thresh, 1);

    dataStruct{1,i}.pca.coeff = coeff;
    dataStruct{1,i}.pca.PCs = nldat(dataStruct{1,i}.emg_env_norm.dataSet * dataStruct{1,i}.pca.coeff,'domainIncr',Ts);
    dataStruct{1,i}.pca.SVs = latent; % latent_n;
    dataStruct{1,i}.pca.explainedVar = explainedVar;
    dataStruct{1,i}.pca.nPC = nPC;
    dataStruct{1,i}.pca.weights = coeff(:,1:nPC);
    dataStruct{1,i}.pca.mainPCs = nldat(dataStruct{1,i}.emg_env_norm.dataSet * dataStruct{1,i}.pca.weights,'domainIncr',Ts);
    dataStruct{1,i}.pca.svThresh = pca_thresh;
end

%% Remove the transient phase of the data only if transient time is non-zero
if transient_time ~= 0
    n_transient = fix(transient_time / Ts) + 1;
    n_transient_d = fix(transient_time / (decimation * Ts)) + 1;
    
    %-- Remove transient phase from the raw data
    for i = 1:nOP
        dataStruct{1,i}.raw.pos = dataStruct{1,i}.raw.pos(n_transient:end,1);
        dataStruct{1,i}.raw.tq = dataStruct{1,i}.raw.tq(n_transient:end,1);
        dataStruct{1,i}.raw.gm = dataStruct{1,i}.raw.gm(n_transient:end,1);
        dataStruct{1,i}.raw.gl = dataStruct{1,i}.raw.gl(n_transient:end,1);
        dataStruct{1,i}.raw.sol = dataStruct{1,i}.raw.sol(n_transient:end,1);
        dataStruct{1,i}.raw.ta = dataStruct{1,i}.raw.ta(n_transient:end,1);
    end

    %-- Remove transient phase from the decimated data
    for i = 1:nOP
        dataStruct{1,i}.decimated.pos = dataStruct{1,i}.decimated.pos(n_transient_d:end,1);
        dataStruct{1,i}.decimated.tq = dataStruct{1,i}.decimated.tq(n_transient_d:end,1);
        dataStruct{1,i}.decimated.gm = dataStruct{1,i}.decimated.gm(n_transient_d:end,1);
        dataStruct{1,i}.decimated.gl = dataStruct{1,i}.decimated.gl(n_transient_d:end,1);
        dataStruct{1,i}.decimated.sol = dataStruct{1,i}.decimated.sol(n_transient_d:end,1);
        dataStruct{1,i}.decimated.ta = dataStruct{1,i}.decimated.ta(n_transient_d:end,1);
    end
    
    %-- Remove transient phase from all other signals
    dataStruct{1,1}.emg_raw = dataStruct{1,1}.emg_raw(n_transient:end,:);
    dataStruct{1,1}.emg_abs = dataStruct{1,1}.emg_abs(n_transient:end,:);
    dataStruct{1,1}.emg_env = dataStruct{1,1}.emg_env(n_transient:end,:);
    dataStruct{1,1}.emg_abs_norm = dataStruct{1,1}.emg_abs_norm(n_transient:end,:);
    dataStruct{1,1}.emg_env_norm = dataStruct{1,1}.emg_env_norm(n_transient:end,:);
    
    %-- Remove transient phase from principal components
    dataStruct{1,1}.pca.PCs = dataStruct{1,1}.pca.PCs(n_transient:end,:);
    dataStruct{1,1}.pca.mainPCs = dataStruct{1,1}.pca.mainPCs(n_transient:end,:);

end

end