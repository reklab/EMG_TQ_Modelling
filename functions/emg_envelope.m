function envEMG = emg_envelope(rawEMG,bpFilter,lpFilter,varargin)
% This function calculates EMG envelope from raw EMG data using the
% following steps:
% 1) Bandpass filtering using bpFilter to eliminate motion artifacts. If it is empty, it means to skip this step. 
%       + If the experiment does not include movement, no need for this.
%       + Some people call this step Highpass filtering, which is equivalent to setting the BF filter higher frequency to Nyquist frequency (Fs/2)
% 2) Full-wave rectification
% 3) Lowpass filtering using lpFilter with two options for the filter: either an iir digital filter (MATLAB's digitalFilter class) or an IRF filter (NLID's irf)

options={{'filter_phase' 'non_zero_phase' 'filter phase. Available options are: zero_phase and non_zero_phase'}
        };
        
if arg_parse(options,varargin)
    return
end

%% Get the number of EMG signals
[nSamples,nSignals] = size(rawEMG);

%% Extract the data part of the raw EMG signals
U = rawEMG.dataSet;
Ts = rawEMG.domainIncr;
Fs = 1/Ts;
bpfU = zeros(nSamples,nSignals);
lpfU = zeros(nSamples,nSignals);

%% Step-1: BP filtering
if ~isempty(bpFilter)
    for i = 1:nSignals
        % bpfU(:,i) = filtfilt(bpFilter,U(:,i)); 
        bpfU(:,i) = filt(bpFilter,U(:,i)); 
    end
else
    bpfU = U;
end

%% Step-2&3: Full-Wave Rectification + Low-Pass filtering
filterType = class(lpFilter);
for i = 1:nSignals
    
    switch filterType
        case 'digitalFilter'
            epsilon = 1e-6;
            if abs(Fs - lpFilter.SampleRate) > epsilon % Fs~=lpFilter.SampleRate
                disp('The sampling rate of the IIR filter does not match the sampling rate of the signal it must filter!')
                envEMG = [];
                return
            else
                switch filter_phase
                    case 'zero_phase'
                        lpfU(:,i) = filtfilt(lpFilter,abs(bpfU(:,i))); 
                    case 'non_zero_phase'
                        lpfU(:,i) = filter(lpFilter,abs(bpfU(:,i)));
                    otherwise
                        disp('This filter_phase option is not defined! The available options are: zero_phase and non_zero_phase.')
                        return;
                end
            end

        case 'irf'
            epsilon = 1e-6;
            if abs(Ts - lpFilter.domainIncr) > epsilon % Ts ~= lpFilter.domainIncr  % Updated because domainIncr in original data was not precisely 0.001s
                disp('The sampling time of the IRF filter does not match the sampling time of the signal it must filter!')
                envEMG = [];
                return
            else
                temp = nlsim(lpFilter,nldat(abs(bpfU(:,i)),'domainIncr',Ts));
                lpfU(:,i) = temp.dataSet;
                
            end
    end

end

%% Assigning the function's output
envEMG = nldat(lpfU,'domainIncr',Ts);

end