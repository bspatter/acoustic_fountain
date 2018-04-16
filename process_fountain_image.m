%
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

% 2018-03-12 : Fountain images taken during camera demo
if true
   fdir='E:/brandon/research/acoustic_fountain/';
   fname = '2018-03-12_fountain_data_camera_demo.xlsx';   
   basename = strip_extension(fname);
   outdir = [fdir 'figs/'];
   fps = 20e3;
    
   % Rotation matrix based on location of the horizon 
   princ_ax_rot = atan((387-345)/649); %rotate clockwise by this much to align x and y with horizontal and vertical
   Rotation_Matrix = [cos(princ_ax_rot), -sin(princ_ax_rot); sin(princ_ax_rot), cos(princ_ax_rot)];

   % Water data fountain 1
   WD1 = xlsread([fdir fname],1); % Get the data      
   frame1 = WD1(:,1);
   t1 = (frame1-frame1(1)) / fps;
   x1 = WD1(:,2);
   y1 = abs(WD1(:,3)- WD1(1,3));
   
   % rotated
   xy1p = (Rotation_Matrix*[x1, y1]')';
   
   % Water data fountain 
   WD2 = xlsread([fdir fname],2); % Get the data   
   frame2 = WD2(:,1);
   x2 = WD2(:,2);
   y2 = abs(WD2(:,3)- WD2(1,3));
   t2 = (frame2-frame2(1)) / fps;

   % Milk data fountain 1
   MD1 = xlsread([fdir fname],4); % Get the data   
   frame_m1 = MD1(:,1);
   t_m1 = (frame_m1-frame_m1(1)) / fps;
   x_m1 = MD1(:,2);
   y_m1 = abs(MD1(:,3)- MD1(1,3));
   
   % Milk data fountain 1
   MD2 = xlsread([fdir fname],5); % Get the data   
   frame_m2 = MD2(:,1);
   t_m2 = (frame_m2-frame_m2(1)) / fps;
   x_m2 = MD2(:,2);
   y_m2 = abs(MD2(:,3)- MD2(1,3));
   

   hh = figure(1);
   axes; hold on;
   plt(1) = plot(t1*1000, y1); 
   plt(2) = plot(t2*1000, y2);
   plt(3) = plot(t_m1*1000, y_m1,'--');
   plt(4) = plot(t_m2*1000, y_m2,'--');
      
   xlabel('time (ms)')
   ylabel('Fountain heigh (pixels)')
   
   mytitle = strrep(basename,'_',' ')
   title(mytitle)
   
   myleg = legend(plt, {'Water', '', 'Milk', ''},'location', 'best')
   
   box on
   spiffyp(hh)
   
   export_fig(hh,[outdir basename],'-png')
   savefig(hh, [outdir basename '.fig'], 'compact')
   
%  Note Drops lost per fountain ([4 5 4 4])
   
   
   
% Original fountain images taken on my phone
elseif false
    fdir='E:/brandon/research/acoustic_fountain/fountain/';
    fps = 240; %frames per second


    nfig=0;

    image_nums=[82,83];
    Nimages = length(image_nums);

    myimages=cell(1,Nimages);

    cropbox_normalized = [0.4, 0.0, 0.2, 1.0]; %[Xmin, Ymin, width, height]

    for n = 1:Nimages
        fname=sprintf('%sframes%04d.png',fdir,image_nums(n));

        % Read image
        Im_n=imread(fname);    
        Im_size = size(Im_n);    


        % crop image
        cropbox = [Im_size(2)*cropbox_normalized(1), Im_size(1)*cropbox_normalized(2), Im_size(2)*cropbox_normalized(3), Im_size(1)*cropbox_normalized(4)];
        Im_n = imcrop(Im_n,cropbox);

        %assign processed image to cell
        myimages{n} = Im_n;


    end

    if true
        nfig=nfig+1;
        figure(nfig)
    %     imshowpair(myimages{:},'diff')
        subplot(1,2,1)
        imshowpair(myimages{:},'montage')

        subplot(1,2,2)    
        imshowpair(myimages{:},'diff')

    end
end