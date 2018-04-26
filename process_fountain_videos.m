function [] = process_fountain_videos()
% Syntax:  [output] = function_name(inputs)
%
% Inputs: 
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
% Copyright (C) 2018, Brandon Patterson
%
% Last revision: 03-01-2018

clc; clear; close all;

% fdir='E:/brandon/research/acoustic_fountain/';

strip_extension = @(mystr) mystr(1:(find(mystr=='.',1,'last')-1));

% 2018-04-16 : Fountain images taken on 2018-04-05
% based on matlab tutorial for motion-based multiple object tracking for computer vision toolbox

CaseCode = 'Da';
basedir = 'E:/brandon/research/acoustic_fountain/acoustic_fountain_videos/';

if true
    fdir=sprintf('%s%s/',basedir,CaseCode);
    fnames = ls([fdir '*.avi']);
    fname = fnames(1,:);
    fpath = [fdir,fname];
    
    
    
end


