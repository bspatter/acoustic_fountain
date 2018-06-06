load('time_pressure_example.mat','time','Pressure');

YP = fft(Pressure);
L=length(Pressure);
YP2 = abs(YP/L);
YP1 = YP2(1:(floor(L/2)+1));
YP1(2:end-1) = 2*YP1(2:end-1);

sampling_frequency = 1/mean(diff(time));

frequency_spectrum = sampling_frequency*(0:(L/2))/L;


pulse_frequency = frequency_spectrum(YP1 == max(YP1));

% FFT plot of pulse
if true
    figure(5)
    plot(frequency_spectrum,YP1)
    title('Single-Sided Amplitude Spectrum of p_{pulse}(t)')
    xlabel('f (Hz)')
    ylabel('|P(f)|')
end

figure(6)
plot(time,ifft(YP))
