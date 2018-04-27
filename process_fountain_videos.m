function fountain_data = process_fountain_videos()
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

try
    tic
    
    strip_extension = @(mystr) mystr(1:(find(mystr=='.',1,'last')-1));
    
    % 2018-04-16 : Fountain images taken on 2018-04-05
    % based on matlab tutorial for motion-based multiple object tracking for computer vision toolbox
    
    CaseCode = 'Db';
    basedir = 'E:/brandon/research/acoustic_fountain/acoustic_fountain_videos/';
    
    fdir=sprintf('%s%s/',basedir,CaseCode);
    fnames = ls([fdir '*.avi']); %All video file names in the relevant directory
    NN = length(fnames);
    fountain_data=struct('filename','','filepath','','CaseCode',CaseCode,'CaseNumber','','date_of_experiment','','time_of_experiment','',...
        'm_per_pix',0,'fountain_height_pix',[],'fountain_height',[],'y_interface',[],'y_fountain',[],'time_raw',[],'time',[],'frame_number',[],'camera_data',struct([]));
    fountain_data = repmat(fountain_data,[NN,1]);
    
    
    
    for nn = 1:NN
        fname = fnames(nn,:); %Specific video file name
        fpath = [fdir,fname]; % Path to video file
        
        if contains(fname,'Ruler'); continue; end
        
        v =VideoReader(fpath);
        v.CurrentTime = 0; % just for debugging for now;
        
        frame = readFrame(v);
        frame0 = frame;
        
        frame_index = {89:size(frame0,1), 129:640};
        
        vidframe0=frame0(frame_index{:}); %extract only the video (lose the border)  %
        vidframe = vidframe0;
        vidframe0b = imbinarize(vidframe0);
        vidframeb = imbinarize(vidframe);
        vidframezeros = zeros(size(vidframe0));
        
        [fountain_height_pix,frame_number,y_interface,y_fountain]=deal(zeros(50e3,1));
        
        ii = 0;
        
        while hasFrame(v)
            ii = ii+1; % counter
            
            %read a frame and extract the part of it with the image (scrap the photron border)
            frame = readFrame(v);
            vidframe = frame(frame_index{:});
            
            % build a grid of indices the size of the frame
            [vidframe_ind_x, vidframe_ind_y] = meshgrid(1:size(vidframe,1), 1:size(vidframe,1));
            
            % binarize the image, based on pixel intensity
            vidframeb = imbinarize(vidframe);
            
            % clean up the image, to identify only the fountain.
            vidframeb_clean = vidframeb;
            vidframeb_clean = imfill(vidframeb_clean,'holes'); %Fill the holes in the middle
            first_left_hole = find(~vidframeb_clean,1);
            if first_left_hole < size(vidframeb_clean)/2
                vidframeb_clean = imfill(vidframeb_clean,first_left_hole);
            end
            
            % find and clean weird left edge noise
            edge_x_threshold = floor(size(vidframeb_clean,2)/10); %pixles left of this not counted as part of fountain
            noise_bounds = sum(vidframeb(:,1:edge_x_threshold),2);% Look at the left 10% of the frame and sum
            noise_bounds0 = [find(noise_bounds~=max(noise_bounds),1,'first') find(noise_bounds==min(noise_bounds),1,'first') ];
            vidframeb_clean2 = vidframeb_clean;
            intf_maxy_threshold = max(noise_bounds0)-3;
            vidframeb_clean2(intf_maxy_threshold,(1:edge_x_threshold)) = true;
            vidframeb_clean2(1:intf_maxy_threshold,1) = true;
            vidframeb_clean2=imfill(vidframeb_clean2,'holes');
            leftnoise = imsubtract(vidframeb_clean2,vidframeb_clean);
            
            vidframeb_clean = logical((vidframeb_clean-leftnoise));
            
            
            % Determine the y-location of top and bottom of fountain
            y_interface(ii) = max(max(vidframe_ind_y(vidframeb_clean == 1)));
            y_fountain(ii) = min(min(vidframe_ind_y(~vidframeb_clean == 1)));
            
            % Identify location of drops leaving the top frame and prevent them from registering as part
            % of the fountain by manually removing them w/ imfill which doesn't notice edge holes automatically.
            while y_fountain(ii) == 1                
                tophole_ind = size(vidframeb_clean,1)*(find(flipud(~imrotate(vidframeb_clean,90)),1)-1)+1; % find the location of the hole at the top
                vidframeb_clean = imfill(vidframeb_clean,tophole_ind);
                y_fountain(ii) = min(min(vidframe_ind_y(~vidframeb_clean == 1)));
            end
            
            frame_number(ii) = ii;
            fountain_height_pix(ii) = y_interface(ii) - y_fountain(ii);
            
            
            % Plots and stuff
            if false %|| ii > 133
                figure(1)
                try % Commented commands should increase performance, but some internal bug is breaking things after a set number of frames;
                    %                 if exist('hh','var')
                    %                     hh.CData = vidframeb;
                    %                 else
                    hh = imshow(vidframeb);
                    %                 end
                    hold on
                    %                 if exist('p1','var')
                    %                     p1.YData = intf_y*[1,1];
                    %                     p2.YData = fountain_y*[1,1];
                    %                 else
                    p1 = plot([0,size(vidframe_ind_x,2)], y_interface(ii)*[1,1],'r','linewidth',2);
                    p2 = plot([0,size(vidframe_ind_x,2)], y_fountain(ii)*[1,1],'r','linewidth',2);
                    %                 end
                    drawnow
                catch
                    continue
                end
            end
            
            
            
%             if ii > 2
%                 break
%             end
        end
        
        % remove trailing zeros and data from far before the fountain
        flast = find(fountain_height_pix,1,'last');
        
        % remove initial
        fountain_threshold = 0.04*max(fountain_height_pix-fountain_height_pix(1)); %minimum height at which a fountain is initiated)
        f0 = max([1, find((fountain_height_pix-fountain_height_pix(1))>fountain_threshold,1,'first')-10]); % Frame number 10 frames before the height exceeds the threshold height based on fountain max (or first frame)
        
        fountain_height_pix = fountain_height_pix(f0:flast);
        frame_number = frame_number(f0:flast);
        y_interface = y_interface(f0:flast);
        y_fountain = y_fountain(f0:flast);
        time = frame_number/20000;
        
        % Save fountain data into a structure
        experiment_info = strsplit(fname,'_');
        fountain_data(nn).filename = fname;
        fountain_data(nn).filepath = fpath;
        fountain_data(nn).CaseNumber = experiment_info{1};
        fountain_data(nn).date_of_experiment = experiment_info{2};
        fountain_data(nn).time_of_experiment = strtrim(experiment_info{3});
        fountain_data(nn).time_of_experiment = fountain_data(nn).time_of_experiment(isstrprop(fountain_data(nn).time_of_experiment,'digit'));
        fountain_data(nn).fountain_height_pix = fountain_height_pix;
        fountain_data(nn).m_per_pix = getscale(fountain_data(nn).CaseCode,fountain_data(nn).date_of_experiment);
        fountain_data(nn).fountain_height = fountain_data(nn).fountain_height_pix * fountain_data(nn).m_per_pix;
        fountain_data(nn).y_interface = y_interface;
        fountain_data(nn).y_fountain = y_fountain;
        fountain_data(nn).time_raw = time;
        fountain_data(nn).time = time - time(1);
        fountain_data(nn).frame_number = frame_number;
        fountain_data(nn).camera_data = readcih(strrep(fname,'.avi','.cih'));
        
        toc
        nn
    end
    
catch myerror
    rethrow(myerror)
end

if false
    hhh = figure(3); %axes; hold on;
    jj=1; plot(fountain_data(jj).time*1e3,fountain_data(jj).fountain_height*1e3); hold on
    jj=2; plot(fountain_data(jj).time*1e3,fountain_data(jj).fountain_height*1e3)
    jj=3; plot(fountain_data(jj).time*1e3,fountain_data(jj).fountain_height*1e3)
    xlabel('time (ms)')
    ylabel('Fountain height (mm)')
    spiffyp(hhh)
    
elseif true  
    hhh = figure(4); %axes; hold on;
    axes;
    mycolors = get(gca,'ColorOrder')
    for jj = 1:42;        
        pp(jj) = plot(fountain_data(jj).time*1e3,fountain_data(jj).fountain_height*1e3,'Color',mycolors(floor((jj-1)/6)+1,:)); 
        hold on        
    end
    xlabel('time (ms)')
    ylabel('Fountain height (mm)')
    legend(pp(1:6:end),'-00 dB', '-02 dB', '-04 dB', '-06 dB', '-08 dB', '-10 dB', '-12 dB')
    spiffyp(hhh)    
end

save(sprintf('fountain_data_%s_%s.mat',CaseCode,fountain_data(nn).date_of_experiment),'fountain_data')

end

function meter_per_pixel = getscale(casecode,casedate)
if strcmp(casecode, 'Db') && strcmp(casedate,'20180424')
    meter_per_pixel = 0.01 / sqrt((543-188)^2+(395-395)^2);
end
end

% function camera_data = read_cih_file(cih_fpath)

