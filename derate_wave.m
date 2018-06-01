function [Pressure_out] = derate_wave(time,Pressure,depth, attenuation_dB_cmMHz)
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
% clear functions; clearvars;
% Define constants 
f0 = 1e6; % Reference frequency (Hz), at which all reference quantities are defined
c0 = 1540; % Reference speed of sound (m/s) at f0

alpha0_dB_cm = attenuation_dB_cmMHz; % 0.3 dB/cm/MHz = 3e-5 dB/m/Hz (reference attenuation coefficient)
alpha0_dB = alpha0_dB_cmMHz*100/1e6; % 0.3 dB/m/Hz (reference attenuation coefficient)
alpha0 = alpha0_dB*log(10)/10;% 1 Np = 20 log10(e) dB ~= 8.685889638 dB % For intensity derating
alpha0_p = alpha0_dB*log(10)/20;% 1 Np = 20 log10(e) dB ~= 8.685889638 dB % For pressure derating
d = depth;%0.0125; % Distance from measurement to transducer (I went with the depth of the interface)


% Bring the wave into the frequency domain
time = time-time(1); % start time from 0;
T=range(time); % Total time
fs = 1/mean(diff(time)); % sampling frequency
N = round(T*fs); %number of data points
N = round(N/2)*2; %force it to be even
T = N/fs; %new total time domain length, if N changed
t = (0:(N-1))/fs; % time vector (max(t) is 1 bin short of T)
f = (0:(N/2))/N*fs; %frequency vector (min(f)=0, the DC component. max(f) = fs/2 exactly)

% Match pressure to output variables
if length(t)==length(Pressure)
    y = Pressure;%function(t); %hypothetical time domain vector, length(y)=N
else
    y = interp1(time,Pressure, t);%function(t); %hypothetical time domain vector, length(y)=N
end
Y = bfft(y); %Get pressure in frequency domain, P(omega), using Brian's fft function. This Uses exp(i omega t) convention for forward transform. length(Y)=N/2+1

% Acoustic velocity as a function of frequency - used in the derivation of attenuation
% cf= c0+(c0/pi)^2*alpha0*log(frequency_spectrum/f0);

% Get the dominant frequency and location
[~,dominant_frequency_index]=max(abs(Y));
dominant_frequency = f(dominant_frequency_index);

% Define the transfer function
betaprime=@(ff) alpha0_p.*ff.*(1-((2*pi*1j*log(ff./f0))./(pi^2+c0*alpha0_p*log(ff./f0))) );
AmplitudeCorrection = 10.^(-(betaprime(f+1e-15)).*d./20);%function(f); %this is where 10^(-attenuation(f)/10) would go, except you'd use 10^(-attenuation(f)/20) for a pressure field and the other one for a power/intensity

Y_derated = Y.*AmplitudeCorrection;
y_derated =  bifft(Y_derated);


%%%%%%%%%%% Try a second form of derating (remove all frequencies above the dominant, then scale by ratio of original energy to "linearized" waveform
Y_linearized = Y;
Y_linearized(f>dominant_frequency*1.5)=0;
y_linearized = bifft(Y_linearized);
energy_ratio = trapz(t,y.^2)/trapz(t,y_linearized.^2);
y_linearized_scaled = y_linearized*sqrt(energy_ratio);
new_energy_ratio = trapz(t,y.^2)/trapz(t,y_linearized_scaled.^2); % should be 1
assert(round(new_energy_ratio,2)==1, 'Linearized wave did not scale correctly. Failed to obey conservation of energy.')

y_linearized_scaled_derated = y_linearized_scaled*exp(-alpha0_p*dominant_frequency*d);


% FFT plot of pulse
if true
    figure(1);clf(1)
    %axes; hold on;
    h1=subplot(1,2,1); hold on
    plot(t*1e6,y/1e6);
    plot(t*1e6,yNew/1e6); hold off;

    ylabel('Pressure (MPa)')
    xlabel('time (\mus)')
    title('Wideband derating')
    spiffyp(1)
    box on

    legend('Push Pulse', sprintf('Derated %g dB/cm/MHz',alpha0_dB_cmMHz))

    text(0.5,9, 'd=1.25cm','FontSize',16)
    
    h2=subplot(1,2,2); hold on
    plot(t*1e6,y/1e6);
    plot(t*1e6,y_linearized_scaled_derated/1e6); hold off;

    ylabel('Pressure (MPa)')
    xlabel('time (\mus)')
    title('Band-pass derating')

    spiffyp(1)
    box on

    legend('Push Pulse', sprintf('Derated %g dB/cm/MHz',alpha0_dB_cmMHz))

    text(0.5,9, 'd=1.25cm','FontSize',16)

elseif false
    figure(5)
    plot(frequency_spectrum,YP1)
    title('Single-Sided Amplitude Spectrum of p_{pulse}(t)')
    xlabel('f (Hz)')
    ylabel('|P(f)|')
    
    % figure(6); clf
    plot(time,ifft([0; Y2])); hold on
    plot(time,ifft([YP_derate]));
end







function [ freqdata ] = bfft( timedata )
    [rows,columns]=size(timedata);
    
    if columns>rows %I want each column to be the complete fourier information
        timedata=transpose(timedata);
    end
    
    L=round(max(size(timedata))/2)*2; %Account for possibility that L is not a power of 2
    NN=min(size(timedata));
    freqdata=zeros(L/2+1,NN);
    
    for i=1:NN
        temptimedata=timedata(:,i); %Assign the ith column
        temp=conj(fft(temptimedata)); %FFT Convention enforced
        freqdata(:,i)=temp(1:L/2+1); %Grab only the information before the folding frequency
    end
    
    if columns>rows %switch it back if necessary
        freqdata=transpose(freqdata);
    end





function [ timedata ] = bifft( freqdata )
    [rows,columns]=size(freqdata);
    
    if columns>rows %I want each column to be the complete fourier information
        freqdata=transpose(freqdata);
    end
    
    L=2*max(size(freqdata))-2;
    NN=min(size(freqdata));
    timedata=zeros(L,NN);
    
    for i=1:NN
        tempfreqdata=freqdata(:,i); %Assign the ith column
        
        tempfreqdata(1)=real(tempfreqdata(1)); %Force the DC bin to be real, otherwise the time domain has a constant imaginary component. As below, a sign change here doesn't change the realness. Not sure what's "most" correct lol...

        L=2*length(tempfreqdata)-2; %Length of time output, such that Length(freq)=L/2+1
        preIFFT=zeros(L,1); %Prepare a full-length vector for IFFT
        preIFFT(1:L/2+1)=tempfreqdata; %Assign the freq data to first half
        temp=conj(flipud(tempfreqdata)); %Flip freq data, conjugate and store
        preIFFT(L/2+1:L)=temp(1:L/2); %Assign the flipped data, skipping the first data point
        preIFFT(L/2+1)=real(preIFFT(L/2+1)); %Force the folding frequency to be real (abs or real part?) - I'm choosing the real part, since having a negative number here doesn't destroy the real-ness of the time data
        timedata(:,i)=ifft(conj(preIFFT)); %Fixed for FFT Convention
    end
    
    if columns>rows %switch it back if necessary
        timedata=transpose(timedata);
    end
        
    
    %Psuedocode for the preIFFT setup:
    %Step 1: 0 0 0 0 0 0 // Create Empty Vector of Zeroes
    %Step 2: 1 2 3 4 0 0 // Assign first half of Freq Data from 1 to L/2+1
    %Step 3: store 4 3 2 1 // FlipLR, Conjugate, and store
    %Step 4: 1 2 3 4 3 2 // Assign this stored info to L/2+1 to L
    %Step 5: 1 2 3 4 3 2 // Force the "4" to become real
    %Step 6: IFFT




