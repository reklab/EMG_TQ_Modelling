function [u_axis,y_axis] = plot_pnlm(model,varargin)
%% This function plots the pnlm model

options = {{'line_width' 2 'The line width of the plotted NL curve'},...
           {'line_style' 'b' 'The line style including color code'},...
           {'plot_type' 'single' 'whether to plot a group of models or a single model'},...
           {'plot_off' false 'turn off the plot and just calculate the function i/o'},...
          };

if arg_parse(options,varargin)
    return
end

nSamples = 400;

switch plot_type
    case 'single'
        nPaths = length(model);
        u_axis = zeros(nSamples,nPaths);
        y_axis = zeros(nSamples,nPaths);
        
        if ~plot_off
            figure;
        end
        for j = 1:nPaths
            u_min = model{j,1}.polyRange(1);
            u_max = model{j,1}.polyRange(2);
            u = nldat(transpose(linspace(u_min, u_max, nSamples)));
        
            y = nlsim(model{j,1},u);

            u_axis(:,j) = u.dataSet; 
            y_axis(:,j) = y.dataSet;
        
            if ~plot_off
                subplot(nPaths,1,j)
                plot(u.dataSet,y.dataSet,line_style,'LineWidth',line_width);
                xlabel(['NL Input-',num2str(j)])
                ylabel('Output')
                title(['NL Model of Input-',num2str(j)])
            end
        end

    case 'group'
        %% To be implemented
        

end

end