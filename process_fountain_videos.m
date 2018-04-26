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

CaseCode = 'Db';
basedir = 'E:/brandon/research/acoustic_fountain/acoustic_fountain_videos/';

if true
    fdir=sprintf('%s%s/',basedir,CaseCode);
    fnames = ls([fdir '*.avi']); %All video file names in the relevant directory
    fname = fnames(1,:); %Specific video file name
    fpath = [fdir,fname]; % Path to video file
    
    v =VideoReader(fpath);
    v.CurrentTime = 5.0; % just for debugging for now;
    
    frame = readFrame(v);    
    frame0 = frame;    
    
    frame_index = {89:size(frame0,1), 129:640};

    vidframe0=frame0(frame_index{:}); %extract only the video (lose the border)  %
    vidframe = vidframe0;
    vidframe0b = imbinarize(vidframe0);
    vidframeb = imbinarize(vidframe);
    vidframezeros = zeros(size(vidframe0));
    
    ii = 0;
    
    while hasFrame(v)
        frame = readFrame(v);
        vidframe = frame(frame_index{:});
        diffimage = imsubtract(vidframe,vidframe0);

        [vidframe_ind_x, vidframe_ind_y] = meshgrid(1:size(vidframe,1), 1:size(vidframe,1));
        
        vidframeb = imbinarize(vidframe);
        vidframeb_clean = vidframeb;

        vidframeb_clean = imfill(vidframeb_clean,'holes'); %Fill the holes in the middle

        first_left_hole = find(~vidframeb_clean,1);
        if first_left_hole < size(vidframeb_clean/2)
%             vidframeb_clean = imfill(vidframeb_clean,first_left_hole); 
        end
        
%         vidframeb_clean = imclearborder(vidframeb_clean);
        
        
        
        fountain_im = imsubtract(vidframe0b,vidframeb);
        fountain_im(fountain_im==-1) = 0; % Not sure what the -1's are, but they need to be removed.
        
        
        intf_y = max(max(vidframe_ind_y(vidframeb_clean == 1)));        
        fountain_y = min(min(vidframe_ind_y(~vidframeb_clean == 1)));
        
        % Deal with drops leaving the top registering as part of the fountain b/c imfill can't handle them.
        if fountain_y == 1
            imshow(vidframeb);
            tophole_ind = size(vidframeb,1)*find(flipud(~imrotate(vidframeb,90)),1)+1; % find the location of the hole at the top
            vidframeb_clean = imfill(vidframeb_clean,tophole_ind);
            fountain_y = min(min(vidframe_ind_y(~vidframeb_clean == 1)));            
        end
        
        figure(1)
        imshow(vidframeb);
        hold on
%         if false && exist('p1','var')
%            p1.YData = intf_y*[1,1];
%            p2.YData = fountain_y*[1,1];
%         else
            p1 = plot([0,size(vidframe_ind_x,2)], intf_y*[1,1],'r','linewidth',2);
            p2 = plot([0,size(vidframe_ind_x,2)], fountain_y*[1,1],'r','linewidth',2);
%         end
        
        
%         fountain_ind = find(fountain_im);
%         fountain_base_ind = max(vidframe_ind_y(fountain_ind)); %Works for base, should plot to check;
        
%         fountain_imb = imbinarize(fountain_im);
%         fountain_only_imb = ~imfill(~fountain_imb,'holes'); %may be useful in isolationg fountain from drops
        
        
%         a = diff(fountain_im,1,1);
%         absa=logical([abs(a); zeros(1,size(a,2))]);
        

        
%         b = find(a);
%         nn = [b; size(fountain_im,1)] - [0; b]
        
%         figure; pcolor(vidframe_ind_x, vidframe_ind_y,[a; zeros(1,512)])
        
%           i = find(diff(fountain_im)) 
%           n = [i numel(fountain_im)] - [0 i]
%           c = arrayfun(@(X) X-1:-1:0, n , 'un',0)
%           y = cat(2,c{:})
        pause(0.05)
    end
    
    
end


