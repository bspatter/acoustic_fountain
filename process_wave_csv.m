function []=process_wave_csv()
%   process_wave_csv.m
%
%
%   %
%   Remarks
%   -------
%
%   Author: Brandon Patterson                   Creation Date: May 11,2018
%
%   Examples
%   --------
close all; clearvars; clc;

us_machine = 'SSI_ELASTO_cal';%'SSI_ELASTO_cal';% 'ACUSONS3000_ELS'
fpath = sprintf('E:/brandon/research/acoustic_fountain/wavedata/%s',us_machine);
csvnames = ls([fpath '/tek*.csv']);
NN = size(csvnames,1);

% [resolved_waves] = [3, 8:14, 20:24]

for n = [3, 8:14, 20:24] %:40%38:45%NN-1;%[6,20:24]%0:24
    close all
    try
        filename = sprintf('%s/tek%04.0fCH1.csv',fpath,n);
        [~,fname,fext] = fileparts(filename);
        outfname = sprintf('%s/%s.xlsx',fpath,fname);
        
        fid=fopen(sprintf('%s',filename),'r');
        if fid<1; fprintf('%s could not be found',filename); return; end
        
        tline = fgetl(fid);
        
        ii = 0;
        
        header_flag = true;
        header = '';
        
        google_green = [60, 186, 84]/255;
        ln = 0;
        while ischar(tline)
            %     sprintf('%s',tline);
            
            ln = ln+1;
            disp(ln)
            
            tline = fgetl(fid);
            
            
            if strcmp(tline,'TIME,CH1')
                header_flag = false;
                sprintf(header);
                %         ii=1;
                %         tline = fgetl(fid);
                %         break
            end
            
            
            if header_flag
                header = sprintf('%s\n%s',header,tline);
            else
                %         dataline=strsplit(
                time_voltage_data=textscan(fid,'%f,%f');
                time = time_voltage_data{1};
                ch1 = time_voltage_data{2};
                break
            end
            
            
            
            ii = ii+1;
            
        end
        
        header_flag = false;
        
        %remove unreal values
        time(isinf(ch1))=[];
        ch1(isinf(ch1))=[];
        
        voltage_mean = mean(ch1);
        voltage_relative = ch1 - voltage_mean;
        
        calibration = 0.022; %V/MPa
        
        time_micro = time*1e6; % Time in microseconds
        
        Pressure_MPa = voltage_relative/calibration;
        Pressure = Pressure_MPa*1e6;
        rho = 1e3;
        c=1.5e3;
        Velocity = Pressure / (rho*c);
        Instant_Intensity = Pressure.*Velocity;
        Instant_Intensity_Watts_per_square_cm = Instant_Intensity * (0.01)^2;
        Intensity = Pressure.^2/(rho*c);
        
        PII = trapz(time, Instant_Intensity)/range(time); %Pulse intensity integral based on Kinsler Eq. (5.9.1)
        PII_cm = trapz(time, Instant_Intensity)/range(time)/1e4; %Pulse intensity integral based on Kinsler Eq. (5.9.1) per centimeter
        
        % Match Doug's calculations
        DH = diff(time).*(Instant_Intensity_Watts_per_square_cm(1:end-1) + Instant_Intensity_Watts_per_square_cm(2:end))/2; %Doug's column H
        DI = cumsum(DH); %Doug's column I
        DPII = sum(DH(1:9900));% Doug's Pulse average intensity value
        DK = DI/DPII; %Something normalized by pulse average intensity
        
        lt10percent = sum(DK<0.1);
        lt90percent = sum(DK<0.9);
        
        
        date_rec = '02-May-2018';
        
        tekfile = fname;
        
        sample_duration_micro = range(time_micro);
        % pulse_duration_micro =
        
        Pmax_MPa = max(Pressure_MPa);
        Pmin_MPa = min(Pressure_MPa);
        Pabs_mid_MPa = (Pmax_MPa - Pmin_MPa)/2;
        
        
        % PSEUDOCODE
        % Find pulse nearest to zero
        % Measure pulse duration using intensity / absolute value
        % Average over pulse
        % Plot both partial pulse and whole waveform w/ zoom region identified
        % Put everything into an excel sheet
        
        % 1. Find pulse nearest to zero
        [~,t0_index] = min(abs(time)); % find t=0
        
        % Pressure_signal_MPa = Pressure_MPa = mode(Pressure_MPa.^2)
        
        Psquared = Pressure.^2;
        
        [pulse_interval,pulse_notes] = tekfile_info(fname,us_machine);
        if ~all(pulse_interval==0)
            pulse_left_boundary_index = find(time < pulse_interval(1),1,'last');
            pulse_right_boundary_index = find(time > pulse_interval(2),1,'first');
            
        elseif false
            % Determine the pulse boundaries based on where the width of the intensity peak at 0.
            % This works for pulses where the intensity has a single clear peak for the pulse,
            % but fails if the intensity returns to the base noise level during the pulse
            threshold_sq_pressure = max(Psquared)*0.01;
            [intensity_peaks,intensity_peak_locs]=findpeaks(Psquared,'Threshold',threshold_sq_pressure);
            % Find the location of the peak at 0
            [~,t0_peak_loc0] = min(abs(t0_index - intensity_peak_locs));
            t0_peak_loc = intensity_peak_locs(t0_peak_loc0);
            
            % Determine what is in the pulse
            pulse_amp_condition = (Psquared > 0.02*Psquared(t0_peak_loc));
            [pulse_left_condition,pulse_right_condition] = deal(false(size(Psquared)));
            pulse_left_condition(1:t0_peak_loc) = true;
            pulse_right_condition(t0_peak_loc:end) = true;
            pulse_left_boundary_index = find(~pulse_amp_condition & pulse_left_condition,1,'last');
            pulse_right_boundary_index = find(~pulse_amp_condition & pulse_right_condition,1,'first');
        elseif true
            
            % Determine the pulse boundaries based on the distance between peaks in the signal
            if false
                threshold_sq_pressure = max(Psquared)*0.01;
                [intensity_peaks,intensity_peak_locs]=findpeaks(Psquared,'MinPeakProminence',threshold_sq_pressure); % This one works-ish
                %
            elseif false % THIS CAPTURES SOME THINGS THE OTHER DOESN'T, TRY DOING THIS WITH union(peaks(P) and peaks(-P)) instead of P squared, and maybe you will get somethign
                % ALT DEF OF PEAK
                [intensity_peaks,intensity_peak_locs]=findpeaks(Psquared);
                noise_peak_locs = intensity_peaks< mean(intensity_peaks);
                intensity_peak_locs(noise_peak_locs)=[];
                intensity_peaks(noise_peak_locs)=[];
            elseif true
                [pos_pressure_peaks,pos_pressure_peak_locs]=findpeaks(Pressure);
                [neg_pressure_peaks,neg_pressure_peak_locs]=findpeaks(-Pressure);
                
                pressure_peak_locs = sort([pos_pressure_peak_locs; neg_pressure_peak_locs]);
                pressure_peaks = Pressure(pressure_peak_locs);
                
                noise_peak_locs = (pressure_peaks.^2)< mean(pressure_peaks.^2); %May want to separate into pos and neg parts if useful spurious points are capture
                
                pressure_peak_locs(noise_peak_locs)=[];
                pressure_peaks(noise_peak_locs)=[];
                
                intensity_peak_locs = pressure_peak_locs; %Ties this with the rest of the code
                
            end
            %
            % Find the location of the peak at 0
            [~,t0_peak_loc0] = min(abs(t0_index - intensity_peak_locs));
            t0_peak_loc = intensity_peak_locs(t0_peak_loc0);
            
            % Determine the pulse edges based on the distance between peaks, which is greater between pulses than within them
            % NOTE: MAY HAVE TO MAKE THE DIFF CIRCULAR TO HANDLE SIGNALS WITH ONE PULSE
            dtpeaks = diff(time(intensity_peak_locs));
            pulse_right_boundary_indecies = [intensity_peak_locs(dtpeaks>mean(dtpeaks)); intensity_peak_locs(end)];
            pulse_left_boundary_indecies = [intensity_peak_locs(1); intensity_peak_locs(circshift(dtpeaks>mean(dtpeaks),1,1))];
            
            [~,pulse_left_boundary_index0] = min(abs(t0_peak_loc - pulse_left_boundary_indecies));
            pulse_left_boundary_index = pulse_left_boundary_indecies(pulse_left_boundary_index0)-1; %-1 added w/o checking
            
            [~,pulse_right_boundary_index0] = min(abs(t0_peak_loc - pulse_right_boundary_indecies));
            pulse_right_boundary_index = pulse_right_boundary_indecies(pulse_right_boundary_index0)+1; %+1 added w/o checking
            
            
            % Used for testing pulse edge detection code above
            if  false
                figure;
                axes; hold on;
                plot(time_micro,Psquared);
                scatter(time_micro(pulse_left_boundary_indecies),Psquared(pulse_left_boundary_indecies))
                scatter(time_micro(pulse_right_boundary_indecies),Psquared(pulse_right_boundary_indecies))
            end
            
            %     scatter(time_micro(intensity_peak_locs(dtpeaks>mean(dtpeaks))),Psquared(intensity_peak_locs(dtpeaks>mean(dtpeaks)))) %This sort of gets the right boundaries of the pulse
            
            %     pulse_amp_condition = (Psquared > 0.02*Psquared(t0_peak_loc));
            %     [pulse_left_condition,pulse_right_condition] = deal(false(size(Psquared)));
            %     pulse_left_condition(1:t0_peak_loc) = true;
            %     pulse_right_condition(t0_peak_loc:end) = true;
            %     pulse_left_boundary_index = find(~pulse_amp_condition & pulse_left_condition,1,'last');
            %     pulse_right_boundary_index = find(~pulse_amp_condition & pulse_right_condition,1,'first');
            
            
        end
        
        pulse_range_indeces = [pulse_left_boundary_index, pulse_right_boundary_index];
        
        pulse_duration = range(time(pulse_range_indeces));
        pulse_duration_micro = pulse_duration*1e6;
        
        pulse_pressure = Pressure(pulse_range_indeces(1):pulse_range_indeces(2));
        pulse_pressure_MPa = pulse_pressure/1e6;
        pulse_intensity = Intensity(pulse_range_indeces(1):pulse_range_indeces(2));
        pulse_time = time(pulse_range_indeces(1):pulse_range_indeces(2));
        Half_Max_Intensity = 0.5*max(pulse_intensity);
        pulse_half_times = intersections(pulse_time,(pulse_intensity-Half_Max_Intensity),[min(pulse_time),max(pulse_time)],[0,0]);
        pulse_half_times_micro = pulse_half_times * 1e6;
        pulse_half_intensities = interp1q(pulse_time,pulse_intensity,pulse_half_times);
        pulse_half_intensities_cm = pulse_half_intensities/1e4;
        FWHM = range(pulse_half_times);
        
        pulse_average_intensity = trapz(pulse_time,pulse_intensity)/range(pulse_time);
        pulse_average_intensity_cm = pulse_average_intensity / 1e4;
        
        
        
        YP = fft(pulse_pressure);
        L=length(pulse_pressure);
        YP2 = abs(YP/L);
        YP1 = YP2(1:(floor(L/2)+1));
        YP1(2:end-1) = 2*YP1(2:end-1);
        
        sampling_frequency = 1/mean(diff(time));
        
        frequency_spectrum = sampling_frequency*(0:(L/2))/L;
        
        
        pulse_frequency = frequency_spectrum(YP1 == max(YP1));
        
        if ~(max(frequency_spectrum)==pulse_frequency)
            MI = abs(min(pulse_pressure_MPa))/sqrt(pulse_frequency/1e6);
            UnderResolved_flag = false;
        else
            UnderResolved_flag = true;
            MI = 0;
        end
        
        
        
        % FFT plot of pulse
        if false
            plot(frequency_spectrum,YP1)
            title('Single-Sided Amplitude Spectrum of p_{pulse}(t)')
            xlabel('f (Hz)')
            ylabel('|P(f)|')
        end
        
        
        
        
        
        mycolors=mycolororder();
        
        g1 = figure(4); spiffyp(g1);
        h1 = subplot(2,2,1);
        plot(time_micro,Pressure_MPa);hold on
        ppr(1) =plot(time_micro(pulse_left_boundary_index)*[1,1],h1.YLim,'linewidth',1,'Color','r','linestyle','--');
        ppr(2) =plot(time_micro(pulse_right_boundary_index)*[1,1],h1.YLim,'linewidth',1,'Color','r','linestyle','--');
        xlim([min(time_micro),max(time_micro)])
        
        xlabel('Time (\mus)')
        ylabel('Pressure (MPa)')
        
        h2 = subplot(2,2,2);
        plot(time_micro,Pressure_MPa);hold on
        plot(time_micro(pulse_left_boundary_index)*[1,1],h1.YLim,'linewidth',1,'Color','r','linestyle','--')
        plot(time_micro(pulse_right_boundary_index)*[1,1],h1.YLim,'linewidth',1,'Color','r','linestyle','--')
        %FIX THIS
        xlimmin = time_micro(max([(pulse_range_indeces(1) + round(range(pulse_range_indeces)*-0.15)),1]));
        xlimmax = time_micro(min([(pulse_range_indeces(2) + round(range(pulse_range_indeces)*+0.15)), length(time_micro)]));
        %     xlim(time_micro(pulse_range_indeces + range(pulse_range_indeces)*2*[-1, +1]))
        %         xlim(time_micro(pulse_range_indeces + round(range(pulse_range_indeces)*0.15*[-1, +1])))
        xlim([xlimmin,xlimmax])
        %     tt(3) = text(time_micro(pulse_left_boundary_index)+round(range(pulse_range_indeces)*0.975), max(h1.YLim)*0.975,sprintf('Pulse\nboundary'),'Color','r','VerticalAlignment','top','HorizontalAlignment','right');
        tt(3) = text(min(xlim()+0.025*range(xlim())), max(h1.YLim)*0.975,sprintf('Pulse\nboundary'),'Color','r','VerticalAlignment','top','HorizontalAlignment','left');
        xlabel('Time (\mus)')
        ylabel('Pressure (MPa)')
        
        h3 = subplot(2,2,3);
        plot(time_micro,Instant_Intensity_Watts_per_square_cm); hold on
        ppr(3) = plot(time_micro(pulse_left_boundary_index)*[1,1],h3.YLim,'linewidth',1,'Color','r','linestyle','--');
        ppr(4) = plot(time_micro(pulse_right_boundary_index)*[1,1],h3.YLim,'linewidth',1,'Color','r','linestyle','--');
        xlim([min(time_micro),max(time_micro)])
        xlabel('Time (\mus)')
        ylabel('Intensity (W / cm^2)')
        
        
        h4 = subplot(2,2,4);
        plot(time_micro,Instant_Intensity_Watts_per_square_cm); hold on
        plot(time_micro(pulse_left_boundary_index)*[1,1],h4.YLim,'linewidth',1,'Color','r','linestyle','--')
        plot(time_micro(pulse_right_boundary_index)*[1,1],h4.YLim,'linewidth',1,'Color','r','linestyle','--')
        plot(pulse_half_times_micro,pulse_half_intensities_cm,'Color',google_green)
        tt(2)=text(pulse_half_times_micro(end)+pulse_duration_micro*0.015,pulse_half_intensities_cm(end),'FWHM','Color',google_green);
        plot(time_micro([pulse_left_boundary_index, pulse_right_boundary_index]),pulse_average_intensity_cm*[1,1],'Color',mycolors(4,:))
        tt(1)=text(time_micro(pulse_left_boundary_index)-pulse_duration_micro*0.015,pulse_average_intensity_cm,'<I_{pa}>','Color',mycolors(4,:),'HorizontalAlignment','right');
        %     xlim(time_micro(pulse_range_indeces + round(range(pulse_range_indeces)*0.15*[-1, +1])))
        xlim([xlimmin,xlimmax])
        if max(pulse_intensity)/max(Intensity) < 0.25; ylim([0, max(pulse_intensity)/1e4]); end
        xlabel('Time (\mus)')
        ylabel('Intensity (W / cm^2)')
        
        % spiffyp(g1)%,'aspect_ratio',2)
        g1fonts = findall(gcf,'type','text');
        FontSize0 = get(g1fonts,'FontSize');
        g1FontSize = mat2cell(cell2mat(FontSize0)+10,ones(size(FontSize0)));
        [g1fonts.FontSize]=deal(g1FontSize{:});
        [h1.FontSize,h2.FontSize,h3.FontSize, h4.FontSize] = deal(15);
        %     g1.Position(3)=1080;
        screensizes = get(groot,'MonitorPositions');
        g1.Position= screensizes(2,:);
                
        if contains(lower(pulse_notes),'whole') || contains(lower(pulse_notes),'---') | strcmp(lower(pulse_notes),'?')
            delete(h2);
            delete(h4);
            delete(ppr);
            h1.Position([1,3])=[0.1,0.875];
            h3.Position([1,3])=[0.1,0.875];
            g1.Position(3)=g1.Position(3)/2;
            pulse_notes =strrep(pulse_notes,'---','');
        end
        drawnow
        
        
        % derated wave stuff
        if ~UnderResolved_flag
            depth_cm = 0.5;
            attenuation_dB_cmMHz = 1.1;
            
            [Pressure_derated, time_derated] = derate_wave(time,Pressure,depth_cm, attenuation_dB_cmMHz);
            
            g2 = figure(7);
            rated_plot = plot(time_micro,Pressure_MPa);hold on            
            plot(time_micro(pulse_left_boundary_index)*[1,1],h1.YLim,'linewidth',1,'Color','r','linestyle','--')
            plot(time_micro(pulse_right_boundary_index)*[1,1],h1.YLim,'linewidth',1,'Color','r','linestyle','--')
            
            derated_plot = plot(time_derated*1e6,Pressure_derated/1e6,'Color',colors('bright pink'),'LineStyle','-.');
            
            %FIX THIS
            xlimmin = time_micro(max([(pulse_range_indeces(1) + round(range(pulse_range_indeces)*-0.15)),1]));
            xlimmax = time_micro(min([(pulse_range_indeces(2) + round(range(pulse_range_indeces)*+0.15)), length(time_micro)]));
            %     xlim(time_micro(pulse_range_indeces + range(pulse_range_indeces)*2*[-1, +1]))
            %         xlim(time_micro(pulse_range_indeces + round(range(pulse_range_indeces)*0.15*[-1, +1])))
            xlim([xlimmin,xlimmax])
            %     tt(3) = text(time_micro(pulse_left_boundary_index)+round(range(pulse_range_indeces)*0.975), max(h1.YLim)*0.975,sprintf('Pulse\nboundary'),'Color','r','VerticalAlignment','top','HorizontalAlignment','right');
            tt(3) = text(min(xlim()+0.025*range(xlim())), max(h1.YLim)*0.975,sprintf('Pulse\nboundary'),'Color','r','VerticalAlignment','top','HorizontalAlignment','left');
            xlabel('Time (\mus)')
            ylabel('Pressure (MPa)')
            
            ll=legend([rated_plot, derated_plot],'Original Wave', sprintf('Derated Wave: d=%3.1gcm, \\alpha=%2.1f dB/cm/MHz',depth_cm,attenuation_dB_cmMHz));
%             ll.Interpreter = 'latex';
            spiffyp(g2)
%             figure(g2)
            if false
                xlswritefig(g2,outfname,'Sheet1','M69')
            end
            
        end
        
        
        if false
            delete(outfname)
            xlswrite(outfname,{'Time (s)', 'Voltage (V)', '<v>', 'volts-<v>','Calibration (V/MPa)','Pressure (MPa)','Intensity (W/cm^2)','Pulse Intensity Integral','<10%','<90%', 'Time (us)'},1,'A1')
            xlswrite(outfname,{time,ch1,voltage_mean,voltage_relative,calibration,Pressure_MPa,Instant_Intensity_Watts_per_square_cm,PII,lt10percent,lt90percent,time_micro},1,'A2')
            xlswrite(outfname,[time,ch1],1,'A2')
            xlswrite(outfname,voltage_relative,1,'D2')
            xlswrite(outfname,[Pressure_MPa, Instant_Intensity_Watts_per_square_cm],1,'F2')
            xlswrite(outfname,time_micro,1,'K2')
            xlswrite(outfname,{sprintf('Pulse Analysis : %s',pulse_notes)},1,'M7')
            xlswrite(outfname,{'TEKFILE:', 'Pulse Durtion (us)', 'P+ (MPa)', 'P- (MPa)','<p>','Ipa (W/cm2)','MI'},1,'M9')
            xlswrite(outfname,{fname, pulse_duration_micro, round(max(pulse_pressure_MPa),3), round(min(pulse_pressure_MPa),3), round((max(pulse_pressure_MPa) - min(pulse_pressure_MPa))/2,3), round(pulse_average_intensity_cm,3),MI},1,'M10')
            xlswritefig(g1,outfname,'Sheet1','M14')
            xlswrite(outfname,{'PROJECT: ACOUSTIC FOUNTAIN'},1,'M1')
            xlswrite(outfname,{'DATE RECORDED: 2018-05-07'},1,'M3')
            if strcmp(us_machine,'SSI_ELASTO_cal'); xlswrite(outfname,{'Trandsucer: 6.7 MHz,   SSI,   SL15-4'},1,'M4')
            elseif strcmp(us_machine,'ACUSONS3000_ELS'); xlswrite(outfname,{'Trandsucer:  9.0 MHz,   ACUSON,   9L4'},1,'M4')
            end
        end
        disp(fname)
        
    catch myerror
        fprintf('Case #%d failed for %s',n,fname)
        disp(myerror.stack)
        fprintf('Line(s): %d\n', [myerror.stack.line])
        fprintf('\n')
        disp(myerror)
        
        
        continue
        
    end
end




function [pulse_interval, notes]= tekfile_info(fname,us_machine)

% Set default values
pulse_interval = [0, 0];
notes = '';

teknum = str2num(fname(4:7))
switch us_machine
    case 'SSI_ELASTO_cal'
        
        switch teknum
            case 0 %'tek0000CH1'
                notes = 'Bmode, -18dB, 100 ms/div';
            case 1 %'tek0001CH1'
                notes = 'B mode, -18dB, 10 ms/div';
            case 2 %'tek0002CH1'
                notes = 'B mode, -18dB, 200 us/div';
                pulse_interval = [-0.2032, +0.4132]*1e-6;
            case 3 %'tek0003CH1'
                notes = 'B mode, -18dB, 200 ns/div';
                pulse_interval = [-0.2032, +0.4132]*1e-6;
            case 4%'tek0004CH1'
                notes = 'SWE at push pulse, -0 dB, 100 ms/div';
                pulse_interval = [3500, 4400]*1e-6;
            case 5%'tek0005CH1'
                notes = 'SWE at push pulse, -0 dB, 10 ms/div';
                pulse_interval = [3550, 4220]*1e-6;
            case 6%'tek0006CH1'
                notes = 'SWE at push pulse, -0 dB, 10 ms/div';
            case 7%'tek0007CH1'
                notes = 'SWE at push pulse, -0 dB, 200 us/div';
                pulse_interval = [-327, +317]*1e-6;
            case 8%'tek0008CH1'
                notes = 'SWE at push pulse, -0 dB, 200 ns/div';
                pulse_interval = [-0.256, +14.95]*1e-6;
            case 9%'tek0009CH1'
                notes = 'SWE at push pulse, -0 dB, 200 ns/div';
                pulse_interval = [-0.31, +2.995]*1e-6;
            case 10%'tek0010CH1'
                notes = 'SWE at push pulse, -0/-2? dB, 200 ns/div';
                pulse_interval = [-0.30, +2.995]*1e-6;
            case 11%'tek0011CH1'
                notes = 'SWE at push pulse, -2 dB, 200 ns/div';
                pulse_interval = [-0.265, +2.995]*1e-6;
            case 12%'tek0012CH1'
                notes = 'SWE at push pulse, -4 dB, 200 ns/div';
                pulse_interval = [-0.30, +2.995]*1e-6;
            case 13%'tek0013CH1'
                notes = 'SWE at push pulse, -6 dB, 200 ns/div';
                pulse_interval = [-0.26, +2.995]*1e-6;
            case 14%'tek0014CH1'
                notes = 'SWE at push pulse, -8 dB, 200 ns/div';
                pulse_interval = [-0.3, +2.995]*1e-6;
            case 15
                notes = 'B mode, SWE box moved to top, -0 dB, ? ms/div';
                pulse_interval = [-4400, -2700]*1e-6;
            case 16
                notes = 'B mode, SWE box moved to left, -0 dB, ? ms/div';
            case 17
                notes = 'B mode, focus and box moved, -0 dB, ? ms/div';
            case 18
                notes = 'Push pulse, -0 dB, 400 us/div';
                pulse_interval = [-164, +486]*1e-6;
            case 19%'tek0019CH1'
                notes = 'SSI Image pulse between push 1 and 2, -0 dB, 400 us/div';
                pulse_interval = [-0.3, +115]*1e-6;
            case 20
                notes = 'SSI Image pulse between push 1 and 2, -0 dB, 200 ns/div';
                pulse_interval = [-0.25, +0.77]*1e-6;
            case 21
                notes = 'SSI Image pulse between push 1 and 2, -2 dB, 200 ns/div ---';
                pulse_interval = [-0.25, +0.77]*1e-6;
            case 22
                notes = 'SSI Image pulse between push 1 and 2, -4 dB, 200 ns/div';
                pulse_interval = [-0.25, +0.77]*1e-6;
            case 23
                notes = 'SSI Image pulse between push 1 and 2, -6 dB, 200 ns/div';
                pulse_interval = [-0.25, +0.77]*1e-6;
            case 24
                notes = 'SSI Image pulse between push 1 and 2, -8 dB, 200 ns/div';
                pulse_interval = [-0.25, +0.77]*1e-6;
                
            otherwise
                pulse_interval = [0, 0];
                warning('Warning: No manual pulse interval prescribed for this case')
        end
        
    case 'ACUSONS3000_ELS'
        switch teknum
            case 0 %'tek0000CH1'
                notes = 'Wide view of whole wave; 0 dB; 200 ms/div';
                
            case 1 %'tek0001CH1'
                notes = 'One set; 0 dB; 10? ms/div';
                pulse_interval = [-0.5, +1.127]*1e5*1e-6;
            case 2 %'tek0002CH1'
                notes = 'Middle of one set - Bmode like pulse; 0 dB; 10? ms/div';
                pulse_interval = [-40, +210]*1e-6;
            case 3 %'tek0003CH1'
                notes = '2nd push; 0 dB; 200 us/div';
                
            case 4%'tek0004CH1'
                notes = '2nd push; 0 dB; 10 us/div';
                pulse_interval = [-1, +79]*1e-6;
            case 5%'tek0005CH1'
                notes = '2nd push; 0 dB; 200 ns/div';
                pulse_interval = [-0.78, +3.15]*1e-6;
            case 6%'tek0006CH1'
                notes = 'Whole thing, without THI; 0 dB; 1 s/div';
            case 7%'tek0007CH1'
                notes = 'Whole thing, without THI; 0 dB; 1 s/div';
            case 8%'tek0008CH1'
                notes = 'Whole thing, with THI, B mode 100%; 0 dB';
            case 9%'tek0009CH1'
                notes = 'Whole thing, without THI, B mode 5%, 9 MHz; 0 dB';
            case 10%'tek0010CH1'
                notes = 'Whole thing, without THI, B mode 5%, 9 MHz; 0 dB; 400 ms/div';
            case 11%'tek0011CH1'
                notes = 'Whole thing, without THI, B mode 5%, 9 MHz; 0 dB; 10 ms/div';
                
            case 12%'tek0012CH1'
                notes = 'Whole thing, without THI, B mode 5%, 9 MHz; 0 dB; 1 ms/div';
                
            case 13%'tek0013CH1'
                notes = 'Whole thing, without THI, B mode 5%, 9 MHz; 0 dB; 100 ms/div';
                
            case 14%'tek0014CH1'
                notes = 'Whole thing, without THI, B mode 5%, 9 MHz; 0 dB; 200 ns/div';
                
            case 15
                notes = '10 ms interior push + imaging?';
                pulse_interval = [3.525, +3.545]*1e4*1e-6;
            case 16
                notes = '10 ms interior push + imaging (actually 15)?';
            case 17
                notes = 'Single set of push pulses + imaging sequence; 40 ms/div';
            case 18
                notes = 'All 4 sets of push pulses + imaging; 200 ms/div';
            case 19%'tek0019CH1'
                notes = 'All 4 sets of push pulses + imaging; 200 ms/div';
            case 20
                notes = 'Big, far right set of #7 of 8 scans, 20 ms/div?';
                pulse_interval = [+8.698e5, +8.805e5]*1e-6;
            case 21
                notes = 'Big, far right set of #7 of 8 scans, 20 ms/div?---';
                pulse_interval = [+8.698e5, +8.782e5]*1e-6;
            case 22
                notes = 'Big, far right set of #7 of 8 scans "18 pre-push, 95 post-push; 400 ns/div?';
                pulse_interval = [+8.732325121495327e5, +8.732332161993769e5]*1e-6;
            case 23
                notes = '"whole box"; 2 m/div';
                pulse_interval = [+8.6968e5, +8.8050e5]*1e-6;%               
            case 24
                notes = 'Push pulse pair; 20 us/div---';                
            case 25
                notes = 'Push pulse pair; 800 ns/div---';
                pulse_interval = [8.480295e5, 8.4803695e5]*1e-6;
            case 26
                notes = 'Push pulse pair; -8 dB; 200 ns/div---';
                pulse_interval = [8.480293e5, 8.4803695e5]*1e-6;
                
            case mat2cell(27:40,1,ones(size(27:40)))
                notes = 'Attempts to measure beam diameter - Not so good---'; 
            otherwise
                pulse_interval = [0, 0];
                warning('Warning: No manual pulse interval prescribed for this case')
        end
        
        
end


