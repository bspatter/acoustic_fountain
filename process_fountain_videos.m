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

clc; clear; %close all;

<<<<<<< HEAD
options=struct();
options.PlotFountain = true;
options.WriteFountainVideo = false && options.PlotFountain; %dependent on options above
options.max_analysis_flag = true; %analyze fountain at max height (get width)
options.saveFountainData = false;

paperfigdir = 'E:/brandon/Dropbox/Research/Papers,etc/papers/2018_fountain/figs/';
=======
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415

tic

strip_extension = @(mystr) mystr(1:(find(mystr=='.',1,'last')-1));

% 2018-04-16 : Fountain images taken on 2018-04-05
% based on matlab tutorial for motion-based multiple object tracking for computer vision toolbox

<<<<<<< HEAD
CaseCode = 'Db';%'Da';
date_of_experiment='20180424';%'20180417';%'20180424';
=======
CaseCode = 'Da';
date_of_experiment='20180521';
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
basedir = 'E:/brandon/research/acoustic_fountain/acoustic_fountain_videos/';

fdir=sprintf('%s%s_%s/',basedir,CaseCode,date_of_experiment);
fnames = ls([fdir '*.avi']); %All video file names in the relevant directory
fnames = fnames(~arrayfun(@(x) contains(x,'Ruler'), string(fnames)),:); % Remove the calibration ruler files
NN = length(fnames);
fountain_data=struct('filename','','filepath','','CaseCode',CaseCode,'CaseNumber','','date_of_experiment','','time_of_experiment','',...
<<<<<<< HEAD
    'm_per_pix',0,'fountain_height_pix',[],'fountain_height',[],'fountain_heigh_max',[],'y_interface',[],'y_fountain',[],'time_raw',[],'time',[],'frame_number',[],'camera_data',struct([]));
=======
    'm_per_pix',0,'fountain_height_pix',[],'fountain_height',[],'y_interface',[],'y_fountain',[],'time_raw',[],'time',[],'frame_number',[],'camera_data',struct([]));
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
fountain_data = repmat(fountain_data,[NN,1]);

%     if strcmp(CaseCode,'Da'); water_flag = true;else water_flag = false; end

% Get the starting frame
load(sprintf('%sfountain_start_end.mat',basedir),'fountain_start_end');


<<<<<<< HEAD

%used for detecting stopping fountain at max
if options.max_analysis_flag 
=======
max_analysis_flag = true
%used for detecting stopping fountain at max
if max_analysis_flag
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
    load(sprintf('%s/fountain_data_%s_%s.mat',fdir,CaseCode,date_of_experiment),'fountain_data')
    
end


<<<<<<< HEAD



% Loop over videos
% tallest_fountain_per_set_Db20180424 = [2,10,14,19,25,31,37];
for nn = 1:NN%:-1:1%:NN
=======
% Loop over videos
for nn = 3:NN%:-1:1%:NN
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
    
    try
        fname = fnames(nn,:); %Specific video file name
        fpath = [fdir,fname]; % Path to video file
        [~,basefname,~] = fileparts(fpath); %Specific video file name w/o extension
        
        if contains(fname,'Ruler'); continue; end
        fprintf('\nnn = %g/%g,      Processing %s\n',nn,NN,fname)
        
<<<<<<< HEAD
        %What case and date are this.
        CaseLabel = strsplit(basefname,'_');
        baseCaseLabel = [CaseLabel{1} ' ' CaseLabel{2}];

        
=======
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
        
        % read the camera data and get the frames per second
        camera_data = readcih(sprintf('%s%s',fdir,strrep(fname,'.avi','.cih')));
        fps = camera_data.Record_Rate_fps;
        
        % Get the scale
        m_per_pixel = getscale(CaseCode,date_of_experiment);
        
        
<<<<<<< HEAD

        
        %% Process the intial video
        %create video object for reading
=======
        %create video object
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
        v =VideoReader(fpath); %#ok<TNMLP>
        
        if isfield(fountain_start_end,basefname)
            v.CurrentTime = fountain_start_end.(basefname).frame_start/v.FrameRate;
        else
            v.CurrentTime = 0; % just for debugging for now;
        end
        
        frame = readFrame(v);
        frame0 = frame;
        
        % Remove border from each frame
        frame_index = {89:size(frame0,1), 129:640};
        vidframe0=frame0(frame_index{:}); %extract only the video (lose the border)  %
        vidframe = vidframe0;
        vidframe0b = imbinarize(vidframe0);
        vidframeb = imbinarize(vidframe);
        
        mid_graythreshold = graythresh(vidframe);
        %Used to fix certain bad lighting situations, which break imbinarize
        bad_lighting_flag = false;
        if any(all(~vidframe0b)) || any(all(vidframe0b)) || bad_lighting_flag %If there are lines of all black or white
            bad_lighting_flag=true;
            mid_graythreshold = graythresh(vidframe(:,round(0.1*size(vidframe)):0.9*size(vidframe)));
            %%
            %             vidframe0b = imbinarize(vidframe0,mid_graythreshold);
            vidframeb = imbinarize(vidframe,mid_graythreshold);
        end
        
        % Not sure if this actually helped anything
        while sum(vidframeb(:)) < numel(vidframeb)*2/3 % If the image is more than 1/3 black, somethign is wrong, reduce the gray threshold
            mid_graythreshold = mid_graythreshold*0.9;
            bad_lighting_flag=true;
            vidframeb = imbinarize(vidframe,mid_graythreshold);
        end
        
<<<<<<< HEAD



        %% Loop over frames to determine fountain properties
        
        % Initialize video writer if saving the updating figure
        if options.WriteFountainVideo
            fout_videoname = [fdir basefname '_processed.mp4'];
            fountain_video = VideoWriter(fout_videoname , 'MPEG-4');
            fountain_video.FrameRate = 15;
            fountain_video.Quality = 90;
            open(fountain_video)
        end
        
=======
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
        % Initialize variables
        [fountain_height_pix,frame_number,y_interface,y_fountain,swell_width_]=deal(zeros(50e3,1));
        vidframezeros = zeros(size(vidframe0));
        ii = 0;
        
<<<<<<< HEAD
        % main loop over frames                
=======
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
        while hasFrame(v)
            ii = ii+1; % counter
            
            %read a frame and extract the part of it with the image (scrap the photron border)
            frame = readFrame(v);
            vidframe = frame(frame_index{:});
            
            % build a grid of indices the size of the frame
            [vidframe_ind_x, vidframe_ind_y] = meshgrid(1:size(vidframe,1), 1:size(vidframe,1));
            
            % binarize the image, based on pixel intensity
            if ~bad_lighting_flag
                vidframeb = imbinarize(vidframe);
            else
                vidframeb = imbinarize(vidframe,mid_graythreshold); %Right now, this keeps using the same gray threshold for every frame, may need to adjust this to recalculate in the loop
            end
            
            % clean up the image, to identify only the fountain.
            vidframeb_clean = vidframeb;
            vidframeb_clean = imfill(vidframeb_clean,'holes'); %Fill the holes in the middle
            first_left_hole = find(~vidframeb_clean,1);
            if first_left_hole < size(vidframeb_clean)/2
                vidframeb_clean = imfill(vidframeb_clean,first_left_hole);
            end
            
            
<<<<<<< HEAD
            % find and clean weird left edge noise
            edge_x_threshold = floor(size(vidframeb_clean,2)/10); %pixles left of this not counted as part of fountain
            
            %CHANGED ABOVE TO USE vidframeb_clean THIS SEEMS TO WORK
            noise_bounds_R = sum(vidframeb_clean(:,1:edge_x_threshold),2);% Look at the left 10% of the frame and sum
            
=======
            
            % SOMETHING LIKE THIS SHOWS PROMISE FOR FILLING IN THE LIGHT SPOTS ON FOUNTAIN
            %              figure; imshow(~imfill(~(vf2b-(imsubtract(vidframe,vidframe0)>10)),'holes'))
            
            
            % find and clean weird left edge noise
            edge_x_threshold = floor(size(vidframeb_clean,2)/10); %pixles left of this not counted as part of fountain
            %             noise_bounds = sum(vidframeb(:,1:edge_x_threshold),2);% Look at the left 10% of the frame and sum
            %CHANGED ABOVE TO USE vidframeb_clean THIS SEEMS TO WORK
            noise_bounds_R = sum(vidframeb_clean(:,1:edge_x_threshold),2);% Look at the left 10% of the frame and sum
            %
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
            % May not work
            noise_bounds0_R = [find(noise_bounds_R~=max(noise_bounds_R),1,'first') find(noise_bounds_R==min(noise_bounds_R),1,'first') ]; %First row with noise, first column with noise
            %
            vidframeb_clean2 = vidframeb_clean;
            intf_maxy_threshold_R = max([1,(max(noise_bounds0_R)-3)]);
            vidframeb_clean2(intf_maxy_threshold_R,(1:edge_x_threshold)) = true; %Fill a specific column
            vidframeb_clean2(1:intf_maxy_threshold_R,1) = true; %Fill first row
            vidframeb_clean2=imfill(vidframeb_clean2,'holes');
            leftnoise = imsubtract(vidframeb_clean2,vidframeb_clean);
            vidframeb_clean = logical((vidframeb_clean-leftnoise));
            
            % find and clean weird left edge noise
            edge_x_threshold = floor(size(vidframeb_clean,2)/10); %pixles left of this not counted as part of fountain
            noise_bounds_R = sum(vidframeb_clean(:,end:-1:(end-edge_x_threshold+1)),2);% Look at the left 10% of the frame and sum
            noise_bounds0_R = [find(noise_bounds_R~=max(noise_bounds_R),1,'first') find(noise_bounds_R==min(noise_bounds_R),1,'first') ]; %First row with noise, first column with noise
            vidframeb_clean2 = vidframeb_clean;
            intf_maxy_threshold_R = max([1,(max(noise_bounds0_R)-3)]);
            vidframeb_clean2(intf_maxy_threshold_R,(end:-1:(end-edge_x_threshold+1))) = true; %Fill a specific column
            vidframeb_clean2(1:intf_maxy_threshold_R,end) = true; %Fill first row
            vidframeb_clean2=imfill(vidframeb_clean2,'holes');
            rightnoise = imsubtract(vidframeb_clean2,vidframeb_clean);
            vidframeb_clean = logical((vidframeb_clean-rightnoise));
            
<<<<<<< HEAD
            % Special processing to try to better capture water, which has weird lighting effects
            if strcmp(CaseCode, 'Da') % water experiments                
%                 background_thresh_y = mean(vidframe0,2); % find background lighting as a function of y
                vidframeb_clean(vidframe>vidframe0*1.3) = 0; % Any pixel that is much brighter than the original background is assumed to be water (reflecting/transmitting light)
                vidframeb_clean(vidframe<vidframe0/1.2) = 0; % Any pixel that is much darker than the original background is assumed to be water (handles curvature at edges of structures)
                vidframeb_clean = ~imfill(~vidframeb_clean,'holes');
                vidframeb_clean = imfill(vidframeb_clean,'holes'); % Remove bright droplets that were filled in by the above process
                
                % Clear droplets exiting the left and right side
                first_left_hole = find(~vidframeb_clean,1);
                if first_left_hole < size(vidframeb_clean,1)/2
                    vidframeb_clean = imfill(vidframeb_clean,first_left_hole);
                end
                
                first_right_hole = find(~fliplr(vidframeb_clean),1);                
                if first_right_hole < size(vidframeb_clean,1)/2
                    vidframeb_clean = fliplr(imfill(fliplr(vidframeb_clean),first_right_hole));
                end
                
                % Deal with reflection at the bottom, giving a false early bottom of the fountain
                % Find the contour that touches both sides to determine the surface
                vidframeb_clean_contours = contourc(double(vidframeb_clean),1);
                k = 0;
                contour_ind = 1;
                contour_data = struct('level',0,'numpoints',0,'x',0,'y',0);
                while contour_ind < length(vidframeb_clean_contours)                   
                    k = k+1;
                    contour_data(k).level = vidframeb_clean_contours(1,contour_ind);
                    contour_data(k).numpoints = vidframeb_clean_contours(2,contour_ind);
                    contour_data(k).x = vidframeb_clean_contours(1,contour_ind+(1:contour_data(k).numpoints));
                    contour_data(k).y = vidframeb_clean_contours(2,contour_ind+(1:contour_data(k).numpoints));
                    contour_ind = contour_ind + contour_data(k).numpoints + 1;
                    
                    % If the contour has data at both sides, you found the surface. Stop lookign.
                    if any(contour_data(k).x) == 1 && any(contour_data(k).x == size(vidframeb_clean,2))
                        left_surface_boundary_ind = ceil(contour_data(k).y(contour_data(k).x==1));
                        right_surface_boundary_ind = ceil(contour_data(k).y(contour_data(k).x==size(vidframeb_clean,2)));
                        right_hole_ind = numel(vidframeb_clean) -size(vidframeb_clean,1) + (right_surface_boundary_ind:size(vidframeb_clean,2)) ;
                        left_hole_ind = (left_surface_boundary_ind:size(vidframeb_clean,2));
                        vidframeb_clean = ~imfill(~vidframeb_clean,right_hole_ind(:));
                        vidframeb_clean = ~imfill(~vidframeb_clean,left_hole_ind(:));                        
                        break
                    end
                end
                
                
                % Tried looking for holes on the right side, but that didn't work well.
                if false
                    vidframe0b_right = vidframe0b(:,end);   
                    color_flip_ind0 = find(diff([vidframe0b_right(1); vidframe0b_right]));  
                    if length(color_flip_ind0) > 1 && mod(length(color_flip_ind0),2)==1 % Does it change color more than once on the right side & an odd number of times
                       right_hole_ind = numel(vidframeb_clean) -size(vidframeb_clean,1) + (color_flip_ind0(end-1):color_flip_ind0(end));
                       vidframeb_clean = ~imfill(~vidframeb_clean,right_hole_ind(:));
                    end
                end
            
            % Still get the contours, even if not looking at water
            elseif strcmp(CaseCode,'Db')
                % Deal with reflection at the bottom, giving a false early bottom of the fountain
                % Find the contour that touches both sides to determine the surface
                vidframeb_clean_contours = contourc(double(vidframeb_clean),1);
                k = 0;
                contour_ind = 1;
                contour_data = struct('level',0,'numpoints',0,'x',0,'y',0);
                while contour_ind < length(vidframeb_clean_contours)                   
                    k = k+1;
                    contour_data(k).level = vidframeb_clean_contours(1,contour_ind);
                    contour_data(k).numpoints = vidframeb_clean_contours(2,contour_ind);
                    contour_data(k).x = vidframeb_clean_contours(1,contour_ind+(1:contour_data(k).numpoints));
                    contour_data(k).y = vidframeb_clean_contours(2,contour_ind+(1:contour_data(k).numpoints));
                    contour_ind = contour_ind + contour_data(k).numpoints + 1;
                    
                    % If the contour has data at both sides, you found the surface. Stop lookign.
                    if any(contour_data(k).x) == 1 && any(contour_data(k).x == size(vidframeb_clean,2))
                        left_surface_boundary_ind = ceil(contour_data(k).y(contour_data(k).x==1));
                        right_surface_boundary_ind = ceil(contour_data(k).y(contour_data(k).x==size(vidframeb_clean,2)));
                        right_hole_ind = numel(vidframeb_clean) -size(vidframeb_clean,1) + (right_surface_boundary_ind:size(vidframeb_clean,2)) ;
                        left_hole_ind = (left_surface_boundary_ind:size(vidframeb_clean,2));
                        vidframeb_clean = ~imfill(~vidframeb_clean,right_hole_ind(:));
                        vidframeb_clean = ~imfill(~vidframeb_clean,left_hole_ind(:));                        
                        break
                    end
                end
                
            end
            
            
=======
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
            
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
            fountain_height = fountain_height_pix(ii)*m_per_pixel;
            
<<<<<<< HEAD

            % Plots fountain in real time
            if options.PlotFountain %|| ii > 133
                
                gg = figure(1);
=======
            % Plots fountain in real time
            if false %|| ii > 133
                figure(1)
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
                try % Commented commands should increase performance, but some internal bug is breaking things after a set number of frames;
                    if exist('hh','var')
                        hh.CData = [vidframeb_clean*255, vidframe];
                    else
                        hh = imshowpair(vidframeb_clean,vidframe,'montage');
                    end
                    hold on
<<<<<<< HEAD
                                        
                    % plot the red lines to show height of interface & the contour of the surface
                    if exist('p1','var')
                        p3.XData = contour_data(end).x;
                        p3.YData = contour_data(end).y;
                        p1.YData = y_interface(ii)*[1,1];
                        p2.YData = y_fountain(ii)*[1,1];
                    else
                        p3 = plot(contour_data(end).x,contour_data(end).y,'b');
=======
                    if exist('p1','var')
                        p1.YData = y_interface(ii)*[1,1];
                        p2.YData = y_fountain(ii)*[1,1];
                    else
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
                        p1 = plot([0,size(vidframe_ind_x,2)], y_interface(ii)*[1,1],'r','linewidth',2);
                        p2 = plot([0,size(vidframe_ind_x,2)], y_fountain(ii)*[1,1],'r','linewidth',2);
                    end
                    
                    if ~exist('p_scale','var')
                        pix_per_mm = 0.001 / m_per_pixel;
                        xscalebar = size(vidframeb_clean,2)*0.9+[0,pix_per_mm];
                        yscalebar = size(vidframeb_clean,1)*0.05*[1,1];
                        p_scale=plot(xscalebar,yscalebar,'k');
<<<<<<< HEAD
                        tt2 = text(0.45,0.97,'1 mm','units','normalized','FontWeight','bold');
                        

                        tt3 = text(0.225,0.97,sprintf('%s',baseCaseLabel),'units','normalized','FontWeight','bold'); % Label the top
=======
                        tt2 = text(0.45,0.97,'1 mm','units','normalized');
                        
                        %What case and date are this.
                        CaseLabel = strsplit(basefname,'_');
                        baseCaseLabel = [CaseLabel{1} ' ' CaseLabel{2}];
                        tt3 = text(0.225,0.97,sprintf('%s',baseCaseLabel),'units','normalized'); % Label the top
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
                        
                        myax=gca;
                        myax.Box='off';
                        myax.Visible='on';
                        myax.XTick(myax.XTick>myax.XLim(2)/2) = [];
                        myax.XTick=[];
                        ytick0 = y_interface(ii);
                        [myax.YTick,YTickLabelOrder] = unique([ytick0:-pix_per_mm:0, ytick0:pix_per_mm:size(vidframeb_clean,2)]);
                        ytick0_index = find(YTickLabelOrder==1);
                        YTickLabels = YTickLabelOrder;
                        YTickLabels(1:ytick0_index) = YTickLabelOrder(1:ytick0_index) - 1;
                        YTickLabels(ytick0_index:end) = 0:-1:-(length(YTickLabelOrder)-ytick0_index);
                        myax.YTickLabel=cellstr(num2str(YTickLabels));
                        ylabel('mm')
<<<<<<< HEAD
                        myax.FontSize=14;
                    elseif ~strcmp(tt3.String,baseCaseLabel)
                        tt3.String = baseCaseLabel;
=======
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
                    end
                    
                    if exist('tt1','var')
                        tt1.String = sprintf('t=%f',frame_number(ii)/fps);
                    else
<<<<<<< HEAD
                        tt1 = text(0.025,0.975,sprintf('t=%f',frame_number(ii)/fps),'units','normalized','FontWeight','bold');
=======
                        tt1 = text(0.025,0.975,sprintf('t=%f',frame_number(ii)/fps),'units','normalized');
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
                    end
                    
                    hold off
                    
                    drawnow
                    
                    % Don't run for longer than necessary (not yet calibrated)
                    if ii>900
                        break;
                    end
                catch myerror
<<<<<<< HEAD
                    fprintf('Broke at line %d; nn = %d, ii = %d, Case: %s\n',myerror.stack.line,nn,ii,fname)
                    disp( getReport( myerror, 'extended', 'hyperlinks', 'on' ) );
=======
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
                    continue
                end
            end
            %             fprintf('%d, ',ii); if mod(ii,20)==0; fprintf('\n');; end
            
            % Analyze the fountain at its peak
<<<<<<< HEAD
            swell_width = 0;% This is just to appease the structure at the end when this section is skipped
            if options.max_analysis_flag 
                %                 disp(fountain_height_pix(ii)/pix_per_mm*1e-3 / max(fountain_data(nn).fountain_height))
                % THE 15e-3 MAY NEED TO BE CALIBRATED FOR EACH EXPERIMENTAL SET
                if round(max(fountain_data(nn).fountain_height(fountain_data(nn).time<15e-3)),7) == round(fountain_height_pix(ii)/pix_per_mm*1e-3,7) %|| true
=======
            if max_analysis_flag && false
                %                 disp(fountain_height_pix(ii)/pix_per_mm*1e-3 / max(fountain_data(nn).fountain_height))
                % THE 15e-3 MAY NEED TO BE CALIBRATED FOR EACH EXPERIMENTAL SET
                if round(max(fountain_data(nn).fountain_height(fountain_data(nn).time<15e-3)),7) == round(fountain_height_pix(ii)/pix_per_mm*1e-3,7) || true
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
%                     fprintf('Fountain Height: %f mm\n',fountain_height_pix(ii)/pix_per_mm)
                    
                    fountain_height_max = fountain_height_pix(ii)/pix_per_mm*1e3;
                    
                    % surface profile
                    [~,surface_profile]=max(~vidframeb_clean);
                    surface_profile = abs(surface_profile-max(surface_profile));
                    
                    
                    
                    
                    % Find the fountain's center index;
                    fountain_metric_x = sum(~vidframeb_clean(1:y_interface(ii),:),1);
                    [~,fountain_max_x_ind] = max(fountain_metric_x);
                    
                    % Find the Fountain's swell width at hmax
                    %                     Lmin = find(fountain_metric_x(1:fountain_max_x_ind)==min(fountain_metric_x(1:fountain_max_x_ind)),1,'last');
                    %                     [~,Rmin] = min(fountain_metric_x(fountain_max_x_ind:end)); Rmin = Rmin + fountain_max_x_ind;
                    %                     fountain_threshold_height_pix = max([0.02*max(fountain_metric_x(1:fountain_max_x_ind)),0.1/pix_per_mm*1e3]);
                    
                    % Fountain threshold based on fountain metric < 5% of max (on L and R sides)
                    fountain_threshold_L1 = 0.02*(max(fountain_metric_x(1:fountain_max_x_ind))-min(fountain_metric_x(1:fountain_max_x_ind))) > (fountain_metric_x(1:fountain_max_x_ind)-min(fountain_metric_x(1:fountain_max_x_ind)));
                    fountain_threshold_R1 = 0.02*(max(fountain_metric_x(fountain_max_x_ind:end))-min(fountain_metric_x(fountain_max_x_ind:end))) > (fountain_metric_x(fountain_max_x_ind:end)-min(fountain_metric_x(fountain_max_x_ind:end)));
                    
                    % Fountain threshold based on zero slope (on L and R sides)
                    %                     fountain_threshold_L2 = [false; (diff(smooth(surface_profile(1:fountain_max_x_ind))) == 0)];
                    %                     fountain_threshold_R2 = [false; (diff(smooth(surface_profile(fountain_max_x_ind:end))) == 0)];
                    fountain_threshold_L2 = false(size(fountain_threshold_L1));
                    fountain_threshold_R2 = false(size(fountain_threshold_R1));
                    
                    
                    Lmin = find(  fountain_threshold_L1 | fountain_threshold_L2,1,'last');
                    Rmin = find(  fountain_threshold_R1 | fountain_threshold_R2,1,'first'); Rmin = Rmin + fountain_max_x_ind;
                    
                    fountain_min_height_threshold = 2e-4;
                    if ~isempty(Lmin) && ~isempty(Rmin) && (fountain_height > fountain_min_height_threshold)
                        if exist('psLR','var') %&& false
                            psLR(1).XData = [Lmin,Lmin];
                            psLR(2).XData = [Rmin,Rmin];
                        else 
                            hold on;
                            % half width -> find half the maximum height, look left and right of fountain_max_x_ind at that height, until you reach air, get that width
                            psLR(1) = plot([Lmin,Lmin],myax.YLim,'r--');
                            psLR(2) = plot([Rmin,Rmin],myax.YLim,'r--'); hold off
                        end
                    
                    if fountain_height > fountain_min_height_threshold % doesn't count as a fountain until its at least this tall
                        swell_width = (Rmin-Lmin)*m_per_pixel;
                    else
                        swell_width = 0;
                    end
                        swell_width_(ii) = swell_width;
%                         fprintf('Swell Width: %f\n',swell_width);
                    end
                    
                    %delete(psLR);
                    % Break at max fountain height
<<<<<<< HEAD
                    if round(max(fountain_data(nn).fountain_height(fountain_data(nn).time<15e-3)),7) == round(fountain_height_pix(ii)/pix_per_mm*1e-3,7) && false
                        break
                    end
                    
                    contour_xmax = contour_data(end).x;
                    contour_ymax = contour_data(end).y;
                    
                    % Used to plot the fountain at its maximum
                    true;
                    if false
                        if strcmp(CaseCode,'Db') && strcmp(date_of_experiment,'20180424')
                            g2 = figure(2); clf(2); imshow(vidframe(108:468,176:410)); hold on; plot([190,190+range(xscalebar)],[24,24],'w'); text(0.8,0.97,'1 mm','units','normalized','FontWeight','bold','Color','w');
                        elseif strcmp(CaseCode,'Da') && strcmp(date_of_experiment,'20180417')                            
                            g2 = figure(2); clf(2); imshow(vidframe(1:468,176:350)); hold on; plot([190,190+range(xscalebar)],[24,24],'w'); text(0.8,0.97,'1 mm','units','normalized','FontWeight','bold','Color','w');
                        end
%                         export_fig(g2,sprintf('%s%s_max',paperfigdir,baseCaseLabel),'-r400','-png');
%                         continue
                    end
                    
                    psLR
                    
=======
                    if round(max(fountain_data(nn).fountain_height(fountain_data(nn).time<15e-3)),7) == round(fountain_height_pix(ii)/pix_per_mm*1e-3,7)
                        break
                    end
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
                end
                
                
            end
            
<<<<<<< HEAD
            if options.WriteFountainVideo
               writeVideo(fountain_video,getframe(gg))
            end
            
            
            
        end
        if options.WriteFountainVideo
            close(fountain_video)
        end
=======
        end
        
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
        % remove trailing zeros and data from far before the fountain
        flast = find(fountain_height_pix,1,'last');
        
        % remove initial
        fountain_threshold = 0.04*max(fountain_height_pix-fountain_height_pix(1)); %minimum height at which a fountain is initiated)
        f0 = max([1, find((fountain_height_pix-fountain_height_pix(1))>fountain_threshold,1,'first')-10]); % Frame number 10 frames before the height exceeds the threshold height based on fountain max (or first frame)
        
        fountain_height_pix = fountain_height_pix(f0:flast);
        frame_number = frame_number(f0:flast);
        y_interface = y_interface(f0:flast);
        y_fountain = y_fountain(f0:flast);
        swell_width_ = swell_width_(f0:flast);
        time = frame_number/fps;
        
        
        
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
<<<<<<< HEAD
        fountain_data(nn).fountain_height_max = max(fountain_height);
=======
        fountain_data(nn).fountain_height_max = fountain_height_max;
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
        fountain_data(nn).swell_width = swell_width;
        fountain_data(nn).y_interface = y_interface;
        fountain_data(nn).y_fountain = y_fountain;
        fountain_data(nn).time_raw = time;
        fountain_data(nn).time = time - time(1);
        fountain_data(nn).frame_number = frame_number;
        fountain_data(nn).camera_data = camera_data;
        fountain_data(nn).swell_width = swell_width_;
<<<<<<< HEAD
        fountain_data(nn).contour_x = contour_xmax;
        fountain_data(nn).contour_y = contour_ymax;
=======
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
        
        % Save the frames to analyze next time (used for faster analysis, only look at frames w/ fountain)
        if f0~=1
            fountain_start_end.(basefname).frame_start=f0;
            fountain_start_end.(basefname).frame_end=flast;
            save(sprintf('%sfountain_start_end.mat',basedir),'fountain_start_end','-append');
        end
        
        toc
        
    catch myerror
<<<<<<< HEAD
        fprintf('Broke at line %d; nn = %d, ii = %d, Case: %s\n',myerror.stack.line,nn,ii,fname)
        disp( getReport( myerror, 'extended', 'hyperlinks', 'on' ) );
        if options.WriteFountainVideo; close(fountain_video); end
%         rethrow(myerror)
=======
        fprintf('Broke at nn = %d, ii = %d, Case: %s\n',nn,ii,fname)
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
        toc
        continue
        rethrow(myerror)
    end
end

try
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
        clf(4)
        axes;
        mycolors = get(gca,'ColorOrder');
<<<<<<< HEAD
        jj0=1;
=======
        jj0=5;
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
        for jj = jj0:42
            pp(jj) = plot(fountain_data(jj).time(1:length(fountain_data(jj).fountain_height))*1e3,fountain_data(jj).fountain_height*1e3,'Color',mycolors(floor((jj-1)/6)+1,:));
            hold on
        end
        xlabel('time (ms)')
        ylabel('Fountain height (mm)')
        legend(pp(jj0:6:end),'-00 dB', '-02 dB', '-04 dB', '-06 dB', '-08 dB', '-10 dB', '-12 dB')
        spiffyp(hhh)
<<<<<<< HEAD
        xlim([0,13.5])
=======
        xlim([0,35])
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
    end
catch
end

% save the data (overwrites existing data, so false by default)
<<<<<<< HEAD
if options.saveFountainData
=======
if false
>>>>>>> db6917b0c2c82bb70368f5903bf65bfb972ba415
    save(sprintf('%s/fountain_data_%s_%s.mat',fdir,CaseCode,fountain_data(nn).date_of_experiment),'fountain_data')
end

end

function meter_per_pixel = getscale(casecode,casedate)
if strcmp(casecode, 'Db') && strcmp(casedate,'20180424')
    meter_per_pixel = 0.01 / sqrt((543-188)^2+(395-395)^2);
elseif strcmp(casecode, 'Da') && strcmp(casedate,'20180521')
    meter_per_pixel = 0.01 / sqrt((603-185)^2+(352-355)^2);
elseif strcmp(casecode, 'Da') && strcmp(casedate,'20180518')
    meter_per_pixel = 0.01 / sqrt((138-561)^2+(340-338)^2);
elseif strcmp(casecode, 'Da') && strcmp(casedate,'20180417')
    meter_per_pixel = 0.01 / sqrt((43-421)^2+(178-177)^2);
end

end

% function camera_data = read_cih_file(cih_fpath)

