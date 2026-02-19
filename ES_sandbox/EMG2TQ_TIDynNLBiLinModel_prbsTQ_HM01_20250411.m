%% This is a code for the analysis and system identification of data I had acquired for EMG-TQ modeling of 
%% subject HM. Data was acquired on 2025-Apr-11. The experiment consisted of voluntary torque variation at  
%% eight (8) isometric conditions; i.e., 8 fixed ankle angular positions equally spaced in the range of motion.
%% Subject varied voluntary torque following a visual cue that was a PRBS between 0 and 40% (MVC) with base switching time of 1 second.

%% The index the random sequence of 8 isometric positions (their index) were:  [6, 8, 2, 7, 4, 3, 1, 5];
clear 
close all         
clc

%% Set user path
userPath = 'C:\Users\ehsan'; %'C:\Users\esobha1';

set(0,'DefaultFigureWindowStyle','docked');

save_flag = 0; %1; %0;

subjName = 'HM01';

%% Add path to my personal MATLAB sandbox:
addpath([userPath,'\Dropbox (NRP)\myPostdoc\REKLAB RA 2021\MATLAB sandbox']);

%% Add path to NLID toolbox and its subdirectories
%-- This is where I have cloned reklab_pulic (NLID toolbox) on my machine at REKLAB: 217W-REKLAB12
nlidBasePath = [userPath,'\reklab_public'];
run([nlidBasePath,'\reklabPaths.m'])

%% Add the path to local clone of the "function" on the repo. These are functions I developed for EMG-Torque modeling.
addpath('C:\Users\ehsan\research_repos\EMG_TQ_Modelling\functions');

%% Path to the folder where all experimental data are stored
filePath = 'C:\Users\ehsan\NRP Dropbox\Ehsan Sobhani Tehrani\myPostdoc\REKLAB RA 2021\Experimental_Data\HM01\';

%% Range of motion and random position sequence info for the subject
maxDF = +0.25; % in rad
maxPF = -0.50; % in rad
nPos = 8;
POS_OP = linspace(maxDF,maxPF,8);
pos_rand_sequence = [6, 1, 7, 4, 5, 8, 3, 2]; 
POS_OP = POS_OP(pos_rand_sequence);
nOP = length(POS_OP);
POS_OP_SORTED = sort(POS_OP);

GM = 1;
GL = 2;
SOL = 3;

%% Identification Hyperparameters
%-- Chebychev nonlinearity
tcheb_order = 3; %5;

%-- Spline nonlinearity
spline_order = 5;
n_knots = 2;

%% Parameters of data and analysis 
pca_thresh = 0.02; %0.01; 0.02; 0.05; 0.1;

input_config = 'emg_env_norm_pc';  % 'emg_env_norm'; 'emg_env_norm_pc';
emg_channels = [1,2,3]; %[1,2,3]; %1; %2; %3; %[1,2]; %[1,3]; %[2,3]; %[1,2,3];

%-- Freeze the number of PCs to use in the modelling
switch input_config
    case 'emg_env_norm_pc'
        nPCs = 1; %1; %2; %3
    case 'emg_env_norm'
        nPCs = length(emg_channels);
end

chan_names = num2str(emg_channels);
chan_names = chan_names(~isspace(chan_names)); %-- Removes spaces in the channel name string

Ts = 0.001;
Fs = 1/Ts;

decimation = 10;

sv_orders = [3,3,3]; %[3,3,3]; [5,5,5];

pc_construct_method = 'all_data'; %'per_segment';

muscle_state_det_method = 'activation_thresh_rate'; % 'activation_thresh_adjust'; 

%++ Threshold to use for finding the contraction and relaxation phases.
switch input_config
    case 'emg_env_norm_pc'
        emg_env_thresh_T = [0.30, 0.35, 0.35, 0.4, 0.4, 0.45, 0.35, 0.3]; % 0.3*ones(1,nOP); 
        emg_env_thresh_R = [0.35, 0.30, 0.32, 0.4, 0.4, 0.35, 0.35, 0.4];  % 0.3*ones(1,nOP); 
    case 'emg_env_norm'
        emg_env_thresh = 0.04; %0.06; %0.04; % This number is from ES06, not optimized for IES01
end

deltaT = 0.75; %0.60; %0.75; 
deltaN = fix(deltaT/Ts);

plot_flag = false; %true; % false;

GREY = [0.8, 0.8, 0.8];

LIGHT_BLUE = [0.68, 0.85, 0.90];
FADE_BLUE = [0.80, 0.90, 0.98];

LIGHT_RED = [1, 0.6, 0.6];
FADE_RED = [1, 0.8, 0.8];

%% EMG envelope filter characteristics
%-- bandpass filter
use_bp_filter = false; %false; true;
bp_order = 20;
bp_f1 = 4;     %-- Hz
bp_f2 = 450;    %-- Hz

%-- filter choice
filter_type = 'irf'; %'irf'; % 'iir'

%-- If we use a digital IIR filter, we only need bandwidth and filter order
BW_res = 0.1;
lpBW_vec = 1; %[0.4, 0.6, 0.8, 1, 1.2]; %[0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2]; %[0.4,0.5,0.6,0.7,0.8]; % 0.4:BW_res:2.0;  % 1.8:BW_res:2.0; 
lpOrder = 2;
lpPassbandRipple = 0.2;

%-- If we use a 2nd order discrete-time IRF, we will have 4 parameters to set 
irfGain = 1;                   %-- filter gain
irfZeta = 0.7;                 %-- filter damping
irfWn_vec = 2*pi*lpBW_vec;     %-- filter natural frequency in rad/s
irfLen_vec = 2.5; %[6, 5, 5, 2.5, 2.5]; %[7,6,5,5,2.5*ones(1,length(lpBW_vec)-2)]; %[5,5,2.5,2.5,2.5]; % 2.5*ones(1,3);  % [5, 2.5*ones(1,8)];                 %-- filter length in seconds
irfTransient = 1; %1; %0;   %TBD: To be incorporated into the code to eliminate the transient phase of the filter before any processing.
irfTs = Ts;
nTransient = fix(irfTransient / Ts) + 1;

%-- Continuity filter
cGain = 1;              %-- filter gain
cZeta = 0.8; %0.7;      %-- filter damping
cBW = 20;               %-- filter bandwidth in Hz
cWn = 2*pi*cBW;         %-- filter natural frequency in rad/s
cLen = 1;               %-- filter length in seconds
cFilter = irf;
cFilter = irf2(cFilter,'g',cGain, 'z',cZeta, 'w',cWn, 'delt',irfTs,'irflen',cLen);
cFilter.domainIncr = Ts;
% figure; plot(cFilter); title('IRF of the LPF to enforce torque continuity')

nBWs = length(lpBW_vec);
lpBW_vec_str = cell(1,nBWs);
POS_OP_str = cell(1,nOP);

for jj = 1:nBWs
    lpBW_vec_str{1,jj} = num2str(lpBW_vec(1,jj));
end

for jj = 1:length(POS_OP)
    POS_OP_str{1,jj} = [num2str(POS_OP(1,jj)),'rad'];
end

n_pos_grid = 200;
pos_OP = linspace(min(POS_OP),max(POS_OP),n_pos_grid)';
avg_pos = (max(pos_OP)+min(pos_OP))/2;
rng_pos = max(pos_OP) - min(pos_OP);
pos_OP_n = (pos_OP - avg_pos)*2/rng_pos;

%-- Initialize 
dataStruct_T = cell(nBWs,1);
dataStruct_R = cell(nBWs,1);

%% Start the main loop, iterating on all possible filter BW options
%% Generating All Data
disp('++++++ STEP-1: Generate all data +++++++')
disp('++++++++++++++++++++++++++++++++++++++++')
for j = 1:nBWs    
    disp('===============================================================')
    disp(['Trying with EMG envelope BW = ',num2str(lpBW_vec(j)), 'Hz'])
    
    %+++ Bandpass filter for EMG envelope calculation
    if use_bp_filter
        bpFilter = designfilt('bandpassiir','FilterOrder',bp_order,...
                              'HalfPowerFrequency1',bp_f1,...
                              'HalfPowerFrequency2',bp_f2,...
                              'SampleRate',Fs);
    else
        bpFilter = [];
    end
    
    figure; freqz(bpFilter,[],Fs);

    %+++ Lowpass filter for EMG envelope calculation
    switch filter_type
        case 'iir'
            lpFilter = designfilt('lowpassiir','FilterOrder',lpOrder, ...
                        'PassbandFrequency',lpBW_vec(1,j),'PassbandRipple',lpPassbandRipple, ...
                        'SampleRate',Fs);
            lpFilterIRF = irf;  
            lpFilterIRF.dataSet = impz(lpFilter);
            lpFilterIRF.domainIncr = Ts;
            figure; plot(lpFilterIRF);
            lpFilterTF = [];


        case 'irf'
            lpFilter = irf;
            lpFilter = irf2(lpFilter,'g',irfGain, 'z',irfZeta, 'w',irfWn_vec(1,j), ...
                            'delt',irfTs,'irflen',irfLen_vec(1,j));
            lpFilterTF = tf(irfGain*irfWn_vec(1,j)*irfWn_vec(1,j), [1,2*irfZeta*irfWn_vec(1,j),irfWn_vec(1,j)*irfWn_vec(1,j)]); 
            % figure; plot(lpFilter); title('IRF of the LPF to calculate EMG envelope')

    end
    
    %++ Test Trial: ES05T
    % subjName = 'ES05T';
    ES05T_FixedPos_TRIALS = [10, 12, 14, 16, 18, 20, 22, 24];
    filesName = 'HM_110425.flb';
    
    dataFile = filesName;
    
    dataStruct_T{j,1} = read_preprocess_data(filePath,dataFile,POS_OP,ES05T_FixedPos_TRIALS,lpFilter,'IES01','decimation',decimation, ...
                                               'plot_raw_data',0,'plot_decimated_data',0,'norm_method','max','pca_thresh',pca_thresh,'data_file_format','flb',...
                                               'bp_filter',bpFilter);

    %++ Retest Trial: ES05R 
    % subjName = 'ES05R';
    ES05R_FixedPos_TRIALS = [26, 28, 30, 32, 34, 36, 38, 40]; 
    filesName = 'HM_110425.flb';
    
    dataFile = filesName;
    
    dataStruct_R{j,1} = read_preprocess_data(filePath,dataFile,POS_OP,ES05R_FixedPos_TRIALS,lpFilter,'IES01','decimation',decimation, ...
                                               'plot_raw_data',0,'plot_decimated_data',0,'norm_method','max','pca_thresh',pca_thresh,'data_file_format','flb');

end

%++ Store generated data into a table
emgEnvBW = lpBW_vec';
dataStruct = dataStruct_T;
dataTable_T = table(emgEnvBW,dataStruct);

dataStruct = dataStruct_R;
dataTable_R = table(emgEnvBW,dataStruct);

%% Use data to identify 1) Single models and 2) Bilinear models per EMG envelope BW and per ankle position
%% The results from the models will be stored as matrices or cells with the size equal to nBW x nOP
%= We have two sets of results:
%           - VAF of torque prediction versus measured
%           - models for each BW and ankle position. Each model is a structure cosisting of:
%                   + pnlm: A cell array of parallel nonlinear polynomial models
%                   + yp: prediction from the model

%% Extracting contract and relax indices based on the sign of the EMG envelopes
%% Use EMG envelopes with a low-pass filter to find these indices
SW_BW = 1.0; % 0.4; % 0.2;
dataStruct_T_1Hz = dataTable_T.dataStruct(dataTable_T.emgEnvBW==SW_BW);
dataStruct_R_1Hz = dataTable_R.dataStruct(dataTable_R.emgEnvBW==SW_BW);

deltaT_r = 0.75; deltaN_r = fix(deltaT_r/Ts);  %-- Time to go backwards from the threshold point of contraction to relaxation
deltaT_c = 0.40; deltaN_c = fix(deltaT_c/Ts);  %-- Time to go backwards from the threshold point of relaxation to contraction

%% 1) Identify bilinear models for Test trials
%++ Initialize the results matrices for Test trials
bilinear_c_vaf = zeros(nBWs,nOP);       %-- VAF of predicting contraction torques
bilinear_r_vaf = zeros(nBWs,nOP);       %-- VAF of predicting relaxation torques
bilinear_t_d_vaf = zeros(nBWs,nOP);     %-- VAF of total torque using discontinuous (swithching) predictions of contraction and relaxation torques
bilinear_t_vaf = zeros(nBWs,nOP);       %-- VAF of total torque using continuous predictions of contraction and relaxation torques

bilinear_c_models = cell(nBWs,nOP);
bilinear_r_models = cell(nBWs,nOP);

disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp('++++ STEP-2: Identifying Single Models for Test Trial +++')
poly_order = tcheb_order; %5; %[5,5]; [5,5,5];
for i = 1:nBWs
    BW = lpBW_vec(i);
    disp('===============================================================')
    disp(['Trying with EMG envelope BW = ',num2str(BW), 'Hz'])
    dataStruct = dataTable_T.dataStruct(dataTable_T.emgEnvBW==BW);
    
    for j = 1:nOP
        switch input_config 
            case 'emg_env_norm_pc' 
                switch_signal_T = dataStruct_T_1Hz{1,1}{1,j}.pca.mainPCs;
            case 'emg_env_norm'
                switch_signal_T = dataStruct_T_1Hz{1,1}{1,j}.emg_env_norm(:,GM);
        end 
        
        switch_signal_dot_T = ddt(switch_signal_T);

        %++ Identify the muscle state
        [contr_indx_T, relax_indx_T] = muscle_state_detector(switch_signal_T,'muscle_state_det_method',muscle_state_det_method,...
                                                                                     'activation_thresh',emg_env_thresh_T(1,j));
        
        %++ Visualize the contraction and relaxation phases:
        time = domain(switch_signal_T);
        switch_signal_T_vec = switch_signal_T.dataSet;
        emg_env_thresh_vec = emg_env_thresh_T(1,j)*ones(size(time));
        
        %++ Set the input / output of contraction phase
        u_c = dataStruct{1,1}{1,j}.pca.PCs(contr_indx_T,1:nPCs);
        y_c = dataStruct{1,1}{1,j}.raw.tq(contr_indx_T,:);
        
        %-- Remove the transient phase
        u_c_id = u_c(nTransient:end,:);
        y_c_id = y_c(nTransient:end,:);
        [pnlm_c,yhat_c_id] = multi_input_pnlm(u_c_id,y_c_id,poly_order);
        yhat_c = sim_pnlm(pnlm_c,u_c);
        
        u_r = dataStruct{1,1}{1,j}.pca.PCs(relax_indx_T,1:nPCs);
        y_r = dataStruct{1,1}{1,j}.raw.tq(relax_indx_T,:);
        %-- Remove the transient phase
        u_r_id = u_r(nTransient:end,:);
        y_r_id = y_r(nTransient:end,:);
        [pnlm_r,yhat_r_id] = multi_input_pnlm(u_r_id,y_r_id,poly_order);
        yhat_r = sim_pnlm(pnlm_r,u_r);

        %-- Merging the predicted torques of contraction and relaxation phases
        uT = dataStruct{1,1}{1,j}.pca.PCs(:,1:nPCs);
        yT = dataStruct{1,1}{1,j}.raw.tq;
        yhatT_d = zeros(size(time));
        yhatT_d(contr_indx_T,:) = yhat_c.dataSet;
        yhatT_d(relax_indx_T,:) = yhat_r.dataSet;

        %-- Create the continuous output
        yhatT = nlsim(cFilter,nldat(yhatT_d,'domainIncr',yT.domainIncr));
        
        figure;
        limP = 0.02;
        subplot(2,1,1)
        plot(time(contr_indx_T),switch_signal_T_vec(contr_indx_T),'r.',time(relax_indx_T),switch_signal_T_vec(relax_indx_T),'b.'); hold on
        plot(time,emg_env_thresh_vec,'--k','LineWidth',1); ylabel('PC Amplitude'); xlabel('Time (s)')
        title(['Ankle Pos = ',num2str(POS_OP_SORTED(j)),'rad'])
        grid
        uT_data = uT.dataSet;
        [umin,umax] = ylim_calc(uT_data(nTransient:end),limP);
        ylim([umin,umax]);
        
        subplot(2,1,2)
        plot(time,yT.dataSet,'k--','LineWidth',1); hold on
        plot(time(contr_indx_T),yhatT_d(contr_indx_T),'r.',time(relax_indx_T),yhatT_d(relax_indx_T),'b.'); hold on
        plot(time,yhatT.dataSet,'m','LineWidth',1)
        ylabel('Torque (Nm)'); xlabel('Time (s)')
        yT_data = yT.dataSet;
        [ymin,ymax] = ylim_calc(yT_data(nTransient:end),limP);
        ylim([ymin,ymax]);
        grid
        legend('measured','pred: conract','pred: relax','pred: continuous')
        xAxisPanZoom

        yhatT_d = nldat(yhatT_d,'domainIncr',yT.domainIncr);
              
        VAF_c = vaf(y_c(nTransient:end,:),yhat_c(nTransient:end,:));
        bilinear_c_vaf(i,j) = VAF_c.dataSet;

        VAF_r = vaf(y_r(nTransient:end,:),yhat_r(nTransient:end,:));
        bilinear_r_vaf(i,j) = VAF_r.dataSet; 

        VAF_t_d = vaf(yT(nTransient:end,:),yhatT_d(nTransient:end,:));
        bilinear_t_d_vaf(i,j) = VAF_t_d.dataSet; 

        VAF_t = vaf(yT(nTransient:end,:),yhatT(nTransient:end,:));
        bilinear_t_vaf(i,j) = VAF_t.dataSet;
      
        %++ Populate the CONTRACTION model structure
        model_c.emgFilter = dataStruct{1,1}{1,j}.emgFilter;
        model_c.emgFilterBW = BW;
        model_c.emgNormFactor = dataStruct{1,1}{1,j}.emg_env_norm_factor;
        model_c.pcWeights = dataStruct{1,1}{1,j}.pca.coeff;
        model_c.id_i = u_c;
        model_c.id_o = y_c;
        model_c.pnlm = pnlm_c;
        model_c.yp = yhat_c;
        model_c.stateIndx = contr_indx_T;
        bilinear_c_models{i,j} = model_c;
        
        %++ Populate the RELAXATION model structure
        model_r.emgFilter = dataStruct{1,1}{1,j}.emgFilter;
        model_r.emgFilterBW = BW;
        model_r.emgNormFactor = dataStruct{1,1}{1,j}.emg_env_norm_factor;
        model_r.pcWeights = dataStruct{1,1}{1,j}.pca.coeff;
        model_r.id_i = u_r;
        model_r.id_o = y_r;
        model_r.pnlm = pnlm_r;
        model_r.yp = yhat_r;
        model_r.stateIndx = relax_indx_T;
        bilinear_r_models{i,j} = model_r;

        dataStruct_T{1,1}{1,j}.models.model_c = model_c;
        dataStruct_T{1,1}{1,j}.models.model_r = model_r;
        dataStruct_T{1,1}{1,j}.models.vaf_c = VAF_c.dataSet;
        dataStruct_T{1,1}{1,j}.models.vaf_r = VAF_r.dataSet;
        dataStruct_T{1,1}{1,j}.models.vaf_d = VAF_t_d.dataSet;
        dataStruct_T{1,1}{1,j}.models.vaf = VAF_t.dataSet;

        dataStruct_T{1,1}{1,j}.models.u = uT;
        dataStruct_T{1,1}{1,j}.models.y = yT;
        dataStruct_T{1,1}{1,j}.models.yp_d = yhatT_d;
        dataStruct_T{1,1}{1,j}.models.yp = yhatT;
    end
end

%++ Plot VAF as a function of envelope bandwidth at various ankle positions 
% figure; 
% subplot(3,1,1)
% plot(lpBW_vec,bilinear_c_vaf,'o--','LineWidth',2); legend(POS_OP_str,'Location','southeast')
% ylabel('VAF (%)');
% title('Bilinear model VAF: CONTRACTION')
% subplot(3,1,2)
% plot(lpBW_vec,bilinear_r_vaf,'o--','LineWidth',2); legend(POS_OP_str,'Location','southeast')
% title('Bilinear model VAF: RELAXATION')
% ylabel('VAF (%)');
% subplot(3,1,3)
% plot(lpBW_vec,bilinear_t_vaf,'o--','LineWidth',2); legend(POS_OP_str,'Location','southeast')
% ylabel('VAF (%)');
% title('Bilinear model VAF: TOTAL')
% xlabel('EMG envelope filter BW (Hz)');

%++ Plot the BW with maximum VAF vs. ankle position
[maxVAF_c,maxIndx_c] = max(bilinear_c_vaf,[],1);
best_bilinear_c_bw_op = lpBW_vec(maxIndx_c);

[maxVAF_r,maxIndx_r] = max(bilinear_r_vaf,[],1);
best_bilinear_r_bw_op = lpBW_vec(maxIndx_r);

[maxVAF_t_d,maxIndx_t_d] = max(bilinear_t_d_vaf,[],1);
best_bilinear_t_d_bw_op = lpBW_vec(maxIndx_t_d);

[maxVAF_t,maxIndx_t] = max(bilinear_t_vaf,[],1);
best_bilinear_t_bw_op = lpBW_vec(maxIndx_t);

%-- Plot best BW and VAF for CONTRACTION
% figure;
% yyaxis left
% plot(POS_OP_SORTED,best_bilinear_c_bw_op,'o--','LineWidth',2);
% xlabel('Ankle Position (rad)')
% ylabel('EMG Envelope Filter BW (Hz)')
% title('CONTRACTION: Best EMG envelope BW at each ankle position')
% 
% yyaxis right
% plot(POS_OP_SORTED,maxVAF_c,'o--','LineWidth',2);
% ylabel('VAF (%)')

%-- Plot best BW and VAF for RELAXATION
% figure;
% yyaxis left
% plot(POS_OP_SORTED,best_bilinear_r_bw_op,'o--','LineWidth',2);
% xlabel('Ankle Position (rad)')
% ylabel('EMG Envelope Filter BW (Hz)')
% title('RELAXATION: Best EMG envelope BW at each ankle position')
% 
% yyaxis right
% plot(POS_OP_SORTED,maxVAF_r,'o--','LineWidth',2);
% ylabel('VAF (%)')

%-- Compare best bandwidth between CONTRACTION and RELAXATION
figure;
subplot(2,1,1)
plot(POS_OP,best_bilinear_c_bw_op,'ro--','LineWidth',2); hold on
plot(POS_OP,best_bilinear_r_bw_op,'bo--','LineWidth',2); hold on
% legend('Contraction','Relaxation')
ylabel('EMG Envelope Filter BW (Hz)')
title('Test Trial: Best EMG envelope BW at each ankle position')

subplot(2,1,2)
plot(POS_OP_SORTED,maxVAF_c,'ro--','LineWidth',2); hold on
plot(POS_OP_SORTED,maxVAF_r,'bo--','LineWidth',2); hold on
plot(POS_OP_SORTED,maxVAF_t_d,'k^--','LineWidth',2); hold on
plot(POS_OP_SORTED,maxVAF_t,'ms--','LineWidth',2);
legend('Contraction','Relaxation','Total Switching','Total Continuous')
title('Test Trial: VAF of best EMG envelope BW at each ankle position')
xlabel('Ankle Position (rad)')
ylabel('VAF (%)')

identVAF_T = maxVAF_t;

%% Plot the best nonlinear models
best_bilinear_model_c_op = cell(1,nOP);
best_bilinear_model_r_op = cell(1,nOP);
for j = 1:nOP
    best_bilinear_model_c_op{1,j} = bilinear_c_models{maxIndx_c(j),j};
    best_bilinear_model_r_op{1,j} = bilinear_r_models{maxIndx_r(j),j};
end

%++ Plot systems for 
POS_OP_BW_c_str = cell(1,nOP);
POS_OP_BW_r_str = cell(1,nOP);
line_style_c = {'b','r','k','g','m','b--','r--','k--'};
line_style_r = {'b--','r--','k--','g--','m--','b','r','k'};
for j = 1:nOP
    POS_OP_BW_c_str{1,j} = ['Contraction: POS=',num2str(POS_OP(1,j)),'rad; ','EMG BW=',num2str(best_bilinear_c_bw_op(1,j)),'Hz'];
    POS_OP_BW_r_str{1,j} = ['Relaxation: POS=',num2str(POS_OP(1,j)),'rad; ','EMG BW=',num2str(best_bilinear_r_bw_op(1,j)),'Hz'];
end

for j=1:nOP
    figure;
    
    %-- Models of Contraction
    %-- Chebychev model
    [u_axis_c, y_axis_c] = plot_pnlm(best_bilinear_model_c_op{1,j}.pnlm,'plot_off',true);
    plot(u_axis_c,y_axis_c,'r--','LineWidth',6); hold on; 
    dataStruct_T{1,1}{1,j}.models.u_axis_c = u_axis_c;
    dataStruct_T{1,1}{1,j}.models.y_axis_c = y_axis_c;

    %-- Scatter plot of input/output
    u_id_c = best_bilinear_model_c_op{1,j}.id_i.dataSet;
    y_id_c = best_bilinear_model_c_op{1,j}.id_o.dataSet;
    plot(u_id_c,y_id_c,'r.','MarkerSize',1); hold on;
    xlabel('PC-1 amplitude')
    ylabel('Nm')
    title(['Model at ankle position = ', num2str(POS_OP_SORTED(1,j)),'rad'])
    
    %-- Models of Relaxation
    %-- Chebychev model
    [u_axis_r, y_axis_r] = plot_pnlm(best_bilinear_model_r_op{1,j}.pnlm,'plot_off',true);
    plot(u_axis_r,y_axis_r,'b--','LineWidth',6); hold on;
    dataStruct_T{1,1}{1,j}.models.u_axis_r = u_axis_r;
    dataStruct_T{1,1}{1,j}.models.y_axis_r = y_axis_r;    

    %-- Scatter plot of input/output
    u_id_r = best_bilinear_model_r_op{1,j}.id_i.dataSet;
    y_id_r = best_bilinear_model_r_op{1,j}.id_o.dataSet;
    plot(u_id_r,y_id_r,'b.','MarkerSize',1); hold on;

    legend(['Contraction Chebychev: Order=',num2str(poly_order)], 'Contraction i/o scatter plot',...
           ['Relaxation Chebychev: Order=',num2str(poly_order)], 'Relaxation i/o scatter plot');

end


%% 2) Identify bilinear models for Re-Test trials
%++ Initialize the results matrices for Re-Test trials
bilinear_c_vaf = zeros(nBWs,nOP);       %-- VAF of predicting contraction torques
bilinear_r_vaf = zeros(nBWs,nOP);       %-- VAF of predicting relaxation torques
bilinear_t_d_vaf = zeros(nBWs,nOP);     %-- VAF of total torque using discontinuous (swithching) predictions of contraction and relaxation torques
bilinear_t_vaf = zeros(nBWs,nOP);       %-- VAF of total torque using continuous predictions of contraction and relaxation torques

bilinear_c_models = cell(nBWs,nOP);
bilinear_r_models = cell(nBWs,nOP);

disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp('++++ STEP-2: Identifying Single Models for Re-Test Trial +++')
poly_order = tcheb_order; %5; %[5,5]; [5,5,5];
for i = 1:nBWs
    BW = lpBW_vec(i);
    disp('===============================================================')
    disp(['Trying with EMG envelope BW = ',num2str(BW), 'Hz'])
    dataStruct = dataTable_R.dataStruct(dataTable_R.emgEnvBW==BW);
    
    for j = 1:nOP
        switch input_config 
            case 'emg_env_norm_pc' 
                switch_signal_R = dataStruct_R_1Hz{1,1}{1,j}.pca.mainPCs;
            case 'emg_env_norm'
                switch_signal_R = dataStruct_R_1Hz{1,1}{1,j}.emg_env_norm(:,GM);
        end 
        
        switch_signal_dot_R = ddt(switch_signal_R);

        %++ Identify the muscle state
        [contr_indx_R, relax_indx_R] = muscle_state_detector(switch_signal_R,'muscle_state_det_method',muscle_state_det_method,...
                                                                     'activation_thresh',emg_env_thresh_R(1,j));
        
        %++ Visualize the contraction and relaxation phases:
        time = domain(switch_signal_R);
        switch_signal_R_vec = switch_signal_R.dataSet;
        emg_env_thresh_vec = emg_env_thresh_R(1,j)*ones(size(time));
        
        u_c = dataStruct{1,1}{1,j}.pca.PCs(contr_indx_R,1:nPCs);
        y_c = dataStruct{1,1}{1,j}.raw.tq(contr_indx_R,:);
        %-- Remove the transient phase
        u_c_id = u_c(nTransient:end,:);
        y_c_id = y_c(nTransient:end,:);
        [pnlm_c,yhat_c_id] = multi_input_pnlm(u_c_id,y_c_id,poly_order);
        yhat_c = sim_pnlm(pnlm_c,u_c);
        
        u_r = dataStruct{1,1}{1,j}.pca.PCs(relax_indx_R,1:nPCs);
        y_r = dataStruct{1,1}{1,j}.raw.tq(relax_indx_R,:);
        %-- Remove the transient phase
        u_r_id = u_r(nTransient:end,:);
        y_r_id = y_r(nTransient:end,:);
        [pnlm_r,yhat_r_id] = multi_input_pnlm(u_r_id,y_r_id,poly_order);
        yhat_r = sim_pnlm(pnlm_r,u_r);

        %-- Merging the predicted torques of contraction and relaxation phases
        uT = dataStruct{1,1}{1,j}.pca.PCs(:,1:nPCs);
        yT = dataStruct{1,1}{1,j}.raw.tq;
        yhatT_d = zeros(size(time));
        yhatT_d(contr_indx_R,:) = yhat_c.dataSet;
        yhatT_d(relax_indx_R,:) = yhat_r.dataSet;

        %-- Create the continuous output
        yhatT = nlsim(cFilter,nldat(yhatT_d,'domainIncr',yT.domainIncr));
        
        figure;
        limP = 0.02;
        subplot(2,1,1)
        plot(time(contr_indx_R),switch_signal_R_vec(contr_indx_R),'r.',time(relax_indx_R),switch_signal_R_vec(relax_indx_R),'b.'); hold on
        plot(time,emg_env_thresh_vec,'--k','LineWidth',1); ylabel('PC Amplitude'); xlabel('Time (s)')
        title(['Ankle Pos = ',num2str(POS_OP_SORTED(j)),'rad'])
        grid
        uT_data = uT.dataSet;
        [umin,umax] = ylim_calc(uT_data(nTransient:end),limP);
        ylim([umin,umax]);
        
        subplot(2,1,2)
        plot(time,yT.dataSet,'k--','LineWidth',1); hold on
        plot(time(contr_indx_R),yhatT_d(contr_indx_R),'r.',time(relax_indx_R),yhatT_d(relax_indx_R),'b.'); hold on
        plot(time,yhatT.dataSet,'m','LineWidth',1)
        ylabel('Torque (Nm)'); xlabel('Time (s)')
        yT_data = yT.dataSet;
        [ymin,ymax] = ylim_calc(yT_data(nTransient:end),limP);
        ylim([ymin,ymax]);
        grid
        legend('measured','pred: contract','pred: relax','pred: continuous')
        xAxisPanZoom

        yhatT_d = nldat(yhatT_d,'domainIncr',yT.domainIncr);
              
        VAF_c = vaf(y_c(nTransient:end,:),yhat_c(nTransient:end,:));
        bilinear_c_vaf(i,j) = VAF_c.dataSet;

        VAF_r = vaf(y_r(nTransient:end,:),yhat_r(nTransient:end,:));
        bilinear_r_vaf(i,j) = VAF_r.dataSet; 

        VAF_t_d = vaf(yT(nTransient:end,:),yhatT_d(nTransient:end,:));
        bilinear_t_d_vaf(i,j) = VAF_t_d.dataSet; 

        VAF_t = vaf(yT(nTransient:end,:),yhatT(nTransient:end,:));
        bilinear_t_vaf(i,j) = VAF_t.dataSet;

        %++ Populate the CONTRACTION model structure
        model_c.emgFilter = dataStruct{1,1}{1,j}.emgFilter;
        model_c.emgFilterBW = BW;
        model_c.emgNormFactor = dataStruct{1,1}{1,j}.emg_env_norm_factor;
        model_c.pcWeights = dataStruct{1,1}{1,j}.pca.coeff;
        model_c.id_i = u_c;
        model_c.id_o = y_c;
        model_c.pnlm = pnlm_c;
        model_c.yp = yhat_c;
        model_c.stateIndx = contr_indx_R;
        bilinear_c_models{i,j} = model_c;
        
        %++ Populate the RELAXATION model structure
        model_r.emgFilter = dataStruct{1,1}{1,j}.emgFilter;
        model_r.emgFilterBW = BW;
        model_r.emgNormFactor = dataStruct{1,1}{1,j}.emg_env_norm_factor;
        model_r.pcWeights = dataStruct{1,1}{1,j}.pca.coeff;
        model_r.id_i = u_r;
        model_r.id_o = y_r;
        model_r.pnlm = pnlm_r;
        model_r.yp = yhat_r;
        model_r.stateIndx = relax_indx_R;
        bilinear_r_models{i,j} = model_r;

        dataStruct_R{1,1}{1,j}.models.model_c = model_c;
        dataStruct_R{1,1}{1,j}.models.model_r = model_r;
        dataStruct_R{1,1}{1,j}.models.vaf_c = VAF_c.dataSet;
        dataStruct_R{1,1}{1,j}.models.vaf_r = VAF_r.dataSet;
        dataStruct_R{1,1}{1,j}.models.vaf_d = VAF_t_d.dataSet;
        dataStruct_R{1,1}{1,j}.models.vaf = VAF_t.dataSet;

        dataStruct_R{1,1}{1,j}.models.u = uT;
        dataStruct_R{1,1}{1,j}.models.y = yT;
        dataStruct_R{1,1}{1,j}.models.yp_d = yhatT_d;
        dataStruct_R{1,1}{1,j}.models.yp = yhatT;
    end
end

%++ Plot VAF as a function of envelope bandwidth at various ankle positions 
% figure; 
% subplot(3,1,1)
% plot(lpBW_vec,bilinear_c_vaf,'o--','LineWidth',2); legend(POS_OP_str,'Location','southeast')
% ylabel('VAF (%)');
% title('Bilinear model VAF: CONTRACTION')
% subplot(3,1,2)
% plot(lpBW_vec,bilinear_r_vaf,'o--','LineWidth',2); legend(POS_OP_str,'Location','southeast')
% title('Bilinear model VAF: RELAXATION')
% ylabel('VAF (%)');
% subplot(3,1,3)
% plot(lpBW_vec,bilinear_t_vaf,'o--','LineWidth',2); legend(POS_OP_str,'Location','southeast')
% ylabel('VAF (%)');
% title('Bilinear model VAF: TOTAL')
% xlabel('EMG envelope filter BW (Hz)');

%++ Plot the BW with maximum VAF vs. ankle position
[maxVAF_c,maxIndx_c] = max(bilinear_c_vaf,[],1);
best_bilinear_c_bw_op = lpBW_vec(maxIndx_c);

[maxVAF_r,maxIndx_r] = max(bilinear_r_vaf,[],1);
best_bilinear_r_bw_op = lpBW_vec(maxIndx_r);

[maxVAF_t_d,maxIndx_t_d] = max(bilinear_t_d_vaf,[],1);
best_bilinear_t_d_bw_op = lpBW_vec(maxIndx_t_d);

[maxVAF_t,maxIndx_t] = max(bilinear_t_vaf,[],1);
best_bilinear_t_bw_op = lpBW_vec(maxIndx_t);

%-- Plot best BW and VAF for CONTRACTION
% figure;
% yyaxis left
% plot(POS_OP_SORTED,best_bilinear_c_bw_op,'o--','LineWidth',2);
% xlabel('Ankle Position (rad)')
% ylabel('EMG Envelope Filter BW (Hz)')
% title('CONTRACTION: Best EMG envelope BW at each ankle position')
% 
% yyaxis right
% plot(POS_OP_SORTED,maxVAF_c,'o--','LineWidth',2);
% ylabel('VAF (%)')

%-- Plot best BW and VAF for RELAXATION
% figure;
% yyaxis left
% plot(POS_OP_SORTED,best_bilinear_r_bw_op,'o--','LineWidth',2);
% xlabel('Ankle Position (rad)')
% ylabel('EMG Envelope Filter BW (Hz)')
% title('RELAXATION: Best EMG envelope BW at each ankle position')
% 
% yyaxis right
% plot(POS_OP_SORTED,maxVAF_r,'o--','LineWidth',2);
% ylabel('VAF (%)')

%-- Compare best bandwidth between CONTRACTION and RELAXATION
figure;
subplot(2,1,1)
plot(POS_OP,best_bilinear_c_bw_op,'ro--','LineWidth',2); hold on
plot(POS_OP,best_bilinear_r_bw_op,'bo--','LineWidth',2); hold on
% legend('Contraction','Relaxation')
ylabel('EMG Envelope Filter BW (Hz)')
title('Re-Test Trial: Best EMG envelope BW at each ankle position')

subplot(2,1,2)
plot(POS_OP_SORTED,maxVAF_c,'ro--','LineWidth',2); hold on
plot(POS_OP_SORTED,maxVAF_r,'bo--','LineWidth',2); hold on
plot(POS_OP_SORTED,maxVAF_t_d,'k^--','LineWidth',2); hold on
plot(POS_OP_SORTED,maxVAF_t,'ms--','LineWidth',2);
legend('Contraction','Relaxation','Total Switching','Total Continuous')
title('Re-Test Trial: VAF of best EMG envelope BW at each ankle position')
xlabel('Ankle Position (rad)')
ylabel('VAF (%)')

identVAF_R = maxVAF_t;

%% Plot the best nonlinear models
best_bilinear_model_c_op = cell(1,nOP);
best_bilinear_model_r_op = cell(1,nOP);
for j = 1:nOP
    best_bilinear_model_c_op{1,j} = bilinear_c_models{maxIndx_c(j),j};
    best_bilinear_model_r_op{1,j} = bilinear_r_models{maxIndx_r(j),j};
end

%++ Plot systems for 
POS_OP_BW_c_str = cell(1,nOP);
POS_OP_BW_r_str = cell(1,nOP);
line_style_c = {'b','r','k','g','m','b--','r--','k--'};
line_style_r = {'b--','r--','k--','g--','m--','b','r','k'};
for j = 1:nOP
    POS_OP_BW_c_str{1,j} = ['Contraction: POS=',num2str(POS_OP(1,j)),'rad; ','EMG BW=',num2str(best_bilinear_c_bw_op(1,j)),'Hz'];
    POS_OP_BW_r_str{1,j} = ['Relaxation: POS=',num2str(POS_OP(1,j)),'rad; ','EMG BW=',num2str(best_bilinear_r_bw_op(1,j)),'Hz'];
end

for j=1:nOP
    figure;
    
    %-- Models of Contraction
    %-- Chebychev model
    [u_axis_c, y_axis_c] = plot_pnlm(best_bilinear_model_c_op{1,j}.pnlm,'plot_off',true);
    plot(u_axis_c,y_axis_c,'r--','LineWidth',6); hold on; 
    dataStruct_R{1,1}{1,j}.models.u_axis_c = u_axis_c;
    dataStruct_R{1,1}{1,j}.models.y_axis_c = y_axis_c;

    %-- Scatter plot of input/output
    u_id_c = best_bilinear_model_c_op{1,j}.id_i.dataSet;
    y_id_c = best_bilinear_model_c_op{1,j}.id_o.dataSet;
    plot(u_id_c,y_id_c,'r.','MarkerSize',1); hold on;
    xlabel('PC-1 amplitude')
    ylabel('Nm')
    title(['Model at ankle position = ', num2str(POS_OP_SORTED(1,j)),'rad'])
    
    %-- Models of Relaxation
    %-- Chebychev model
    [u_axis_r, y_axis_r] = plot_pnlm(best_bilinear_model_r_op{1,j}.pnlm,'plot_off',true);
    plot(u_axis_r,y_axis_r,'b--','LineWidth',6); hold on;
    dataStruct_R{1,1}{1,j}.models.u_axis_r = u_axis_r;
    dataStruct_R{1,1}{1,j}.models.y_axis_r = y_axis_r;    

    %-- Scatter plot of input/output
    u_id_r = best_bilinear_model_r_op{1,j}.id_i.dataSet;
    y_id_r = best_bilinear_model_r_op{1,j}.id_o.dataSet;
    plot(u_id_r,y_id_r,'b.','MarkerSize',1); hold on;

    legend(['Contraction Chebychev: Order=',num2str(poly_order)], 'Contraction i/o scatter plot',...
           ['Relaxation Chebychev: Order=',num2str(poly_order)], 'Relaxation i/o scatter plot');

end

%% Plot and compare the identified contraction and relaxation models between Test and Re-Test trials
LineWidth = 3;
u_axis_n = linspace(-1,1,length(dataStruct_T{1,1}{1,j}.models.u_axis_c));

for j = 1:nOP
    figure;
    subplot(2,1,1)
    plot(dataStruct_T{1,1}{1,j}.models.u_axis_c,dataStruct_T{1,1}{1,j}.models.y_axis_c,'r','LineWidth',LineWidth); hold on
    plot(dataStruct_R{1,1}{1,j}.models.u_axis_c,dataStruct_R{1,1}{1,j}.models.y_axis_c,'r--','LineWidth',LineWidth); hold on
    plot(dataStruct_T{1,1}{1,j}.models.u_axis_r,dataStruct_T{1,1}{1,j}.models.y_axis_r,'b','LineWidth',LineWidth); hold on
    plot(dataStruct_R{1,1}{1,j}.models.u_axis_r,dataStruct_R{1,1}{1,j}.models.y_axis_r,'b--','LineWidth',LineWidth); 
    legend('Contraction - Test','Contraction - ReTest','Relaxation - Test','Relaxation - ReTest')
    title(['Model at ankle position = ', num2str(POS_OP_SORTED(1,j)),'rad'])
    xlabel('PC-1 Amplitude')
    ylabel('Nm')

    subplot(2,1,2)
    plot(u_axis_n,dataStruct_T{1,1}{1,j}.models.y_axis_c,'r','LineWidth',LineWidth); hold on
    plot(u_axis_n,dataStruct_R{1,1}{1,j}.models.y_axis_c,'r--','LineWidth',LineWidth); hold on
    plot(u_axis_n,dataStruct_T{1,1}{1,j}.models.y_axis_r,'b','LineWidth',LineWidth); hold on
    plot(u_axis_n,dataStruct_R{1,1}{1,j}.models.y_axis_r,'b--','LineWidth',LineWidth); 
    legend('Contraction - Test','Contraction - ReTest','Relaxation - Test','Relaxation - ReTest')
    % title(['Model at ankle position = ', num2str(POS_OP_SORTED(1,j)),'rad'])
    xlabel('Normalized PC-1 Amplitude')
    ylabel('Nm')
end

%% Using model identified with Test data to predict Retest torque and vice versa
simVAF_R_by_T = zeros(1,nOP);
simVAF_T_by_R = zeros(1,nOP);

%% Predicting Retest trial torque using the PC martix and PNLM model identified with Test trial data
for j = 1:nOP
    %-- First, extract the model (PCA weight and NL model) 
    weightT = dataStruct_T{1,1}{1,j}.pca.coeff(:,1:nPCs);
    modelT = dataStruct_T{1,1}{1,j}.models;
    modelT_c = modelT.model_c;
    modelT_r = modelT.model_r;
    % Then, construct the input for Retest trial by reading the normalized EMG envelopes of Retest multiplied by the Test PC matrix
    switch input_config
        case 'emg_env_norm_pc'
            uR = nldat_matrix_prod(dataStruct_R{1,1}{1,j}.emg_env_norm, weightT);
        case 'emg_env_norm'
            uR = dataStruct_R{1,1}{1,j}.emg_env_norm(:,emg_channels);
    end

    %-- Extract the contraction and relaxation phases indices of the Restest trials
    contr_indx_R = dataStruct_R{1,1}{1,j}.models.model_c.stateIndx;
    relax_indx_R = dataStruct_R{1,1}{1,j}.models.model_r.stateIndx;
    
    %-- Dividing the input into contraction and relaxation phases
    uR_c = uR(contr_indx_R,:);
    uR_r = uR(relax_indx_R,:);

    %-- Reading the output of the Retest trials
    yR = dataStruct_R{1,1}{1,j}.raw.tq;
    yR_c = uR(contr_indx_R,:);
    yR_r = uR(relax_indx_R,:);
    
    %========= To add in the future if required
    %-- Remove the transient phase of both inputs and outputs (relaxation and contraction)
    % u_c_sim = u_c(nTransient:end,:);
    % y_c_sim = y_c(nTransient:end,:);
    %=========
    yhatR_c = sim_pnlm(modelT_c.pnlm,uR_c);
    yhatR_r = sim_pnlm(modelT_r.pnlm,uR_r);
    
    yhatR_d = zeros(size(time));
    yhatR_d(contr_indx_R,:) = yhatR_c.dataSet;
    yhatR_d(relax_indx_R,:) = yhatR_r.dataSet;
    yhatR_d = nldat(yhatR_d,'domainIncr',yR.domainIncr);
    yhatR = nlsim(cFilter,yhatR_d);

    simVAF_R_by_T(1,j) = vaf(yR,yhatR).dataSet;

    figure;
    uR_vec = uR.dataSet;
    subplot(2,1,1)
    plot(time(contr_indx_R),uR_vec(contr_indx_R),'r.',time(relax_indx_R),uR_vec(relax_indx_R),'b.'); hold on
    ylabel('PC-1 Amplitude')
    title(['Ankle Pos = ',num2str(POS_OP_SORTED(j)),'rad'])
    grid
    [umin,umax] = ylim_calc(uR_vec(nTransient:end),limP);
    ylim([umin,umax]);
    
    subplot(2,1,2)
    yhatR_d_vec = yhatR_d.dataSet;
    plot(time,yR.dataSet,'k--','LineWidth',1); hold on
    % plot(time(contr_indx_R),yhatR_d_vec(contr_indx_R),'r.',time(relax_indx_R),yhatR_d_vec(relax_indx_R),'b.'); hold on
    plot(time,yhatR.dataSet,'m','LineWidth',1)
    ylabel('Torque (Nm)'); xlabel('Time (s)')
    title(['VAF of predicting Retest trial torque with the Test trial model = ',sprintf('%.1f',simVAF_R_by_T(1,j)),'%'])
    yR_vec = yR.dataSet;
    [ymin,ymax] = ylim_calc(yR_vec(nTransient:end),limP);
    ylim([ymin,ymax]);
    grid
    % legend('measured','pred: conract','pred: relax','pred: continuous')
    legend('measured','pred: continuous')
    xAxisPanZoom
end

%% Predicting Test trial torque using the PC martix and PNLM model identified with Retest trial data
for j = 1:nOP
    %-- First, extract the model (PCA weight and NL model) 
    weightR = dataStruct_R{1,1}{1,j}.pca.coeff(:,1:nPCs);
    modelR = dataStruct_R{1,1}{1,j}.models;
    modelR_c = modelR.model_c;
    modelR_r = modelR.model_r;
    % Then, construct the input for Test trial by reading the normalized EMG envelopes of Test multiplied by the Retest PC matrix
    switch input_config
        case 'emg_env_norm_pc'
            uT = nldat_matrix_prod(dataStruct_T{1,1}{1,j}.emg_env_norm, weightR);
        case 'emg_env_norm'
            uT = dataStruct_T{1,1}{1,j}.emg_env_norm(:,emg_channels);
    end

    %-- Extract the contraction and relaxation phases indices of the Test trials
    contr_indx_T = dataStruct_T{1,1}{1,j}.models.model_c.stateIndx;
    relax_indx_T = dataStruct_T{1,1}{1,j}.models.model_r.stateIndx;
    
    %-- Dividing the input into contraction and relaxation phases
    uT_c = uT(contr_indx_T,:);
    uT_r = uT(relax_indx_T,:);

    %-- Reading the output of the Test trials
    yT = dataStruct_T{1,1}{1,j}.raw.tq;
    yT_c = uT(contr_indx_T,:);
    yT_r = uT(relax_indx_T,:);
    
    %========= To add in the future if required
    %-- Remove the transient phase of both inputs and outputs (relaxation and contraction)
    % u_c_sim = u_c(nTransient:end,:);
    % y_c_sim = y_c(nTransient:end,:);
    %=========
    yhatT_c = sim_pnlm(modelR_c.pnlm,uT_c);
    yhatT_r = sim_pnlm(modelR_r.pnlm,uT_r);
    
    yhatT_d = zeros(size(time));
    yhatT_d(contr_indx_T,:) = yhatT_c.dataSet;
    yhatT_d(relax_indx_T,:) = yhatT_r.dataSet;
    yhatT_d = nldat(yhatT_d,'domainIncr',yT.domainIncr);
    yhatT = nlsim(cFilter,yhatT_d);

    simVAF_T_by_R(1,j) = vaf(yT,yhatT).dataSet;

    figure;
    uT_vec = uT.dataSet;
    subplot(2,1,1)
    plot(time(contr_indx_T),uT_vec(contr_indx_T),'r.',time(relax_indx_T),uT_vec(relax_indx_T),'b.'); hold on
    ylabel('PC-1 Amplitude')
    title(['Ankle Pos = ',num2str(POS_OP_SORTED(j)),'rad'])
    grid
    [umin,umax] = ylim_calc(uT_vec(nTransient:end),limP);
    ylim([umin,umax]);
    
    subplot(2,1,2)
    yhatT_d_vec = yhatT_d.dataSet;
    plot(time,yT.dataSet,'k--','LineWidth',1); hold on
    % plot(time(contr_indx_T),yhatT_d_vec(contr_indx_T),'r.',time(relax_indx_T),yhatT_d_vec(relax_indx_T),'b.'); hold on
    plot(time,yhatT.dataSet,'m','LineWidth',1)
    ylabel('Torque (Nm)'); xlabel('Time (s)')
    title(['VAF of predicting Test trial torque with the Retest trial model = ',sprintf('%.1f',simVAF_T_by_R(1,j)),'%'])
    yT_vec = yT.dataSet;
    [ymin,ymax] = ylim_calc(yT_vec(nTransient:end),limP);
    ylim([ymin,ymax]);
    grid
    % legend('measured','pred: conract','pred: relax','pred: continuous')
    legend('measured','pred: continuous')
    xAxisPanZoom
end

%-- Replace any negative VAF values with a VAF of zero
simVAF_R_by_T(simVAF_R_by_T<0) = 0;
simVAF_T_by_R(simVAF_T_by_R<0) = 0;

%% Total VAF of Test and Retest trials as a function of ankle position
figure;
subplot(1,2,1)
plot(POS_OP_SORTED,identVAF_T,'bo--','LineWidth',2);
hold on
plot(POS_OP_SORTED,simVAF_T_by_R,'ro--','LineWidth',2);
legend('ident in T trial','pred by R trial model')
title(['Test Trial: Model Order = ',num2str(tcheb_order)])
xlabel('Ankle Position (rad)')
ylabel('VAF (%)')
ylim([0,100])

subplot(1,2,2)
plot(POS_OP_SORTED,identVAF_R,'bs--','LineWidth',2);
hold on
plot(POS_OP_SORTED,simVAF_R_by_T,'rs--','LineWidth',2);
legend('ident in R trial','pred by T trial model')
title(['Re-Test Trial: Model Order = ',num2str(tcheb_order)])
xlabel('Ankle Position (rad)')
ylabel('VAF (%)')
ylim([0,100])

%% Save the identified models
fileName = ['models_',subjName,'_bilinear','.mat'];
if save_flag
    save(fileName,'POS_OP_SORTED','identVAF_T','identVAF_R','simVAF_T_by_R','simVAF_R_by_T','dataStruct_T','dataStruct_R');
end
