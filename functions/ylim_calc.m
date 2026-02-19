function [ymin,ymax] = ylim_calc(y,limP)
%% This function calculates the limits for plots based on a fixed percentage of min and max added.
%% The reason we have a separate function for this is that the min and max may have negative or positive signs.
y_range = max(y) - min(y);

ymin = min(y) - limP*y_range;
ymax = max(y) + limP*y_range;




