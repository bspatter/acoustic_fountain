function [hist_lab_colors, hist_lab_colors_labels] = histology_lab_colors()

%
% Syntax:  [output] = function_name(inputs)
%
% Inputs: NONE
%   1. 
%
% Outputs:
%	1. 
% 
% Example: 
% 
% 
% Other m-files required: 
% 
% Subfunctions: 
% 
% MAT-files required: 
%
% See also: 
% 
% Author: Brandon patterson
% email address: i.am.brandon.patterson@gmail.com
% 
% Copyright (C) 2017, Brandon Patterson
%
% Last revision: 02-13-2018

hist_lab_colors = [0.5267, 0.6138; 10.5938,  -18.7697; 13.4934, -4.1837; 3.5312, -0.7602]; % A and B values for colors in lung histology, based on LAB color space.
hist_lab_colors_labels = {'white','dark_purple','light_purple','transparent_purple'};