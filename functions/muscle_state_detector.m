function [contr_indx, relax_indx] = muscle_state_detector(activation_signal,varargin)
%% This function finds the contraction and relaxation phases muscles 
%% from an appropriate surrogate of muscle activation

%-- Descriptions of the methods 
options={{'muscle_state_det_method' 'activation_thresh_rate' 'Available options are: activation_thresh_rate and activation_thresh_adjust'},...
         {'activation_thresh' [] 'threshold applied to the surrogate activation signal'},...
         {'delta_t_c2r' 0.75 'Time to go backwards from the threshold point of contraction to relaxation in seconds'},...
         {'delta_t_r2c' 0.20 'Time to go backwards from the threshold point of relaxation to contraction in seconds'},...  %-- the old value was 0.4s
         {'n_times_smooth' 4000 'The number of times to smooth a signal'},...
        };
        
if arg_parse(options,varargin)
    return
end

if isempty('activation_thresh')
    disp('The default value of activation threshold is empty! Please set its value in your function call.')
    return
end

%-- Extracting parameters of data
Ts = activation_signal.domainIncr;
deltaN_c2r = fix(delta_t_c2r/Ts);
deltaN_r2c = fix(delta_t_r2c/Ts);

%-- Running the main code
switch muscle_state_det_method 
    case 'activation_thresh_adjust'
        contr_indx_init = activation_signal.dataSet > activation_thresh;
        relax_indx_init = activation_signal.dataSet <= activation_thresh;

        %-- Refine the state change points
        contr_indx = contr_indx_init;
        relax_indx = relax_indx_init;

        contr_indx_diff = [contr_indx(1); diff(contr_indx)];
        
        %-- Expanding the relaxation phase and shortening the contraction phase by going back to the start of EMG envelope / PC drop
        indx = find(contr_indx_diff==-1);
        for jj = 1:length(indx)
            if indx(jj)-deltaN_c2r > 0
                contr_indx(indx(jj)-deltaN_c2r:indx(jj)) = false; 
                relax_indx(indx(jj)-deltaN_c2r:indx(jj)) = true; 
            end
        end
        
        %-- Expanding the contraction phase and shortening the relaxation phase by going back to the start of EMG envelope / PC rising edge 
        indx = find(contr_indx_diff==+1);
        for jj = 1:length(indx)
            if indx(jj)-deltaN_r2c > 0
                contr_indx(indx(jj)-deltaN_r2c:indx(jj)) = true; 
                relax_indx(indx(jj)-deltaN_r2c:indx(jj)) = false; 
            end
        end
    
    case 'activation_thresh_rate'  %-- This methods is the default method that uses a combination of threshold and activation rate
        contr_indx_init = activation_signal.dataSet > activation_thresh;
        relax_indx_init = activation_signal.dataSet <= activation_thresh;

        %-- Refine the state change points
        contr_indx = contr_indx_init;
        relax_indx = relax_indx_init;

        contr_indx_diff = [contr_indx_init(1); diff(contr_indx_init)];
        indx_c2r = find(contr_indx_diff==-1);
        indx_r2c = find(contr_indx_diff==+1);

        %-- Calculating the rate of the activation signal
        activation_signal_dot = ddt(activation_signal);
        x = activation_signal.dataSet;
        activation_signal_dot_smo = smo(activation_signal_dot,n_times_smooth);
        xdot = activation_signal_dot_smo.dataSet;
        time = domain(activation_signal_dot_smo);
        
        [~,~,zcIndices] = zerocrossrate(xdot);
        zcIndices = zcIndices';
        zcPoints = find(zcIndices==1);

        figure;
        subplot(2,1,1)
        plot(time,x); hold on; plot(time(indx_c2r),x(indx_c2r),'rv'); hold on; plot(time(indx_r2c),x(indx_r2c),'r^'); title('Activation Signal')
        subplot(2,1,2)
        plot(time,xdot); hold on; plot(time(zcIndices),xdot(zcIndices),'ro'); title('Smoothed Activation Rate')
        xAxisPanZoom

        %-- Expanding the relax phase backward to the start of relaxation, based on the zero-crossing 
        for jj = 1:length(indx_c2r)
            delta_indx = zcPoints - indx_c2r(jj);
            [~, maxInd] = max(delta_indx(delta_indx<=0));  %-- This is the logic to find the closest backward zero-crossing
            contr_indx(zcPoints(maxInd):indx_c2r(jj)) = false; 
            relax_indx(zcPoints(maxInd):indx_c2r(jj)) = true;
        end

        %-- Expanding the contraction phase backward to the start of contraction, based on the zero-crossing 
        for jj = 1:length(indx_r2c)
            delta_indx = zcPoints - indx_r2c(jj);
            [~, maxInd] = max(delta_indx(delta_indx<=0));  %-- This is the logic to find the closest backward zero-crossing
            contr_indx(zcPoints(maxInd):indx_r2c(jj)) = true; 
            relax_indx(zcPoints(maxInd):indx_r2c(jj)) = false;
        end
        % disp('Pause')
    
    case 'activation_rate_r_thresh_c'  %-- This was an intermediary method that I implemented, which uses rate for relaxation and threshold for contraction.
        contr_indx_init = activation_signal.dataSet > activation_thresh;
        relax_indx_init = activation_signal.dataSet <= activation_thresh;

        %-- Refine the state change points
        contr_indx = contr_indx_init;
        relax_indx = relax_indx_init;

        contr_indx_diff = [contr_indx_init(1); diff(contr_indx_init)];
        indx_c2r = find(contr_indx_diff==-1);
        indx_r2c = find(contr_indx_diff==+1);

        %-- Calculating the rate of the activation signal
        activation_signal_dot = ddt(activation_signal);
        x = activation_signal.dataSet;
        activation_signal_dot_smo = smo(activation_signal_dot,n_times_smooth);
        xdot = activation_signal_dot_smo.dataSet;
        time = domain(activation_signal_dot_smo);
        
        [~,~,zcIndices] = zerocrossrate(xdot);
        zcIndices = zcIndices';
        zcPoints = find(zcIndices==1);

        figure;
        subplot(2,1,1)
        plot(time,x); hold on; plot(time(indx_c2r),x(indx_c2r),'rv'); hold on; plot(time(indx_r2c),x(indx_r2c),'r^'); title('Activation Signal')
        subplot(2,1,2)
        plot(time,xdot); hold on; plot(time(zcIndices),xdot(zcIndices),'ro'); title('Smoothed Activation Rate')
        xAxisPanZoom

        %-- Expanding the relax phase backward to the start of relaxation, based on the zero-crossing , 
        for jj = 1:length(indx_c2r)
            delta_indx = zcPoints - indx_c2r(jj);
            [~, maxInd] = max(delta_indx(delta_indx<=0));  %-- This is the logic to find the closest backward zero-crossing
            contr_indx(zcPoints(maxInd):indx_c2r(jj)) = false; 
            relax_indx(zcPoints(maxInd):indx_c2r(jj)) = true;
        end
        
        %-- Expanding the contraction phase and shortening the relaxation phase by going back to the start of EMG envelope / PC rising edge 
        for jj = 1:length(indx_r2c)
            if indx_r2c(jj)-deltaN_r2c > 0
                contr_indx(indx_r2c(jj)-deltaN_r2c:indx_r2c(jj)) = true; 
                relax_indx(indx_r2c(jj)-deltaN_r2c:indx_r2c(jj)) = false; 
            end
        end
        % disp('Pause')
    
    otherwise
        disp('This method is not supported!'); 
        return;
end  % End of switch statement

end  % End of function