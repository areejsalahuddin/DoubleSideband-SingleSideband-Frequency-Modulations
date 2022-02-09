
%%%%%%%%%%%%%%%%%%%%%%%%%%%(1)%%%%%%%%%%%%%%%%%%%%%%%%%
[y , fs] = audioread('eric.wav');
% sound(y,fs);
y_Spectrum = fftshift(fft(y));
f = linspace(-fs/2,fs/2,length(y));
t = linspace(0, length(y)/fs, length(y));
plot_time(y,fs,'Signal in time domain');
plot_frequency(abs(y_Spectrum),fs,'Signal in Frequency Domain');


%%%%%%%%%%%%%%%%%%%%%%%%%%%(2)%%%%%%%%%%%%%%%%%%%%%%%%%
filter1 = ones(1,length(y_Spectrum)).*(f >= -4000 & f <= 4000);



%%%%%%%%%%%%%%%%%%%%%%%%%%%(3)%%%%%%%%%%%%%%%%%%%%%%%%%
filteredSpectrum = y_Spectrum.' .* filter1 ;
ytime_filtered = real(ifft(ifftshift(filteredSpectrum)));
plot_time(ytime_filtered,fs,'Filtered Signal in Time Domain');
plot_frequency(abs(filteredSpectrum),fs,'Filtered Signal in Frequency Domain');
% sound(ytime_filtered,fs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%(4)%%%%%%%%%%%%%%%%%%%%%%%%%
fc = 100000;
FS = 5*fc;
resampledSignal = resample(ytime_filtered,FS,fs); 
f_new = linspace(-FS/2, FS/2,length(resampledSignal));
t_new = linspace(0, length(resampledSignal)/FS, length(resampledSignal));
carrierSignal = cos(2 * pi * fc * t_new); 
modulated_DSB_SC_time = carrierSignal.* resampledSignal;
modulated_DSB_SC_Spectrum = fftshift(fft(modulated_DSB_SC_time));
plot_time(modulated_DSB_SC_time,FS,'DSB-SC modulated Signal in Time Domain');
plot_frequency(abs(modulated_DSB_SC_Spectrum),FS,'DSB-SC modulated Signal in Frequency Domain');


%%%%%%%%%%%%%%%%%%%%%%%%%%%(5)%%%%%%%%%%%%%%%%%%%%%%%%%
filter2 = ones(1,length(modulated_DSB_SC_Spectrum)).*(f_new >= -fc & f_new <= fc); 
modulated_SSB_Spectrum = modulated_DSB_SC_Spectrum.*filter2;
modulated_SSB_time = real(ifft(ifftshift(modulated_SSB_Spectrum)));
plot_time(modulated_SSB_time,FS,'SSB-SC modulated Signal in Time Domain (LSB) Ideal Filter');
plot_frequency(abs(modulated_SSB_Spectrum),FS,'SSB-SC modulated Signal in Frequency Domain (LSB) Ideal Filter');


%%%%%%%%%%%%%%%%%%%%%%%%%%%(6)%%%%%%%%%%%%%%%%%%%%%%%%%
demodulated_time = modulated_SSB_time.*carrierSignal;
filter3 = designfilt('lowpassfir', 'FilterOrder', 8000, 'CutoffFrequency', 4000, 'SampleRate', FS);
coherent_demodulated_time  = filter(filter3, demodulated_time);
coherent_demodulated_freq = fft(fftshift(coherent_demodulated_time));
plot_time(real(coherent_demodulated_time),FS,'Demodulated Signal in Time Domain (Coherent Detection Ideal Filter)');
plot_frequency(abs(real(coherent_demodulated_freq)),FS,'Demodulated Signal in Frequency Domain (Coherent Detection Ideal Filter)');
resampled_time = resample(coherent_demodulated_time, fs, FS);
resampled_freq = fftshift(fft(resampled_time));
plot_time(real(resampled_time),fs,'Demodulated Signal in Time Domain (Coherent Detection Ideal Filter)');
plot_frequency(abs(real(resampled_freq)),fs,'Demodulated Signal in Frequency Domain (Coherent Detection Ideal Filter)');
%sound(real(resampled_time),fs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%(7)%%%%%%%%%%%%%%%%%%%%%%%%%
%%% step (5) %%%
fcn = FS/2;
fc1 = (fc - 4000)/fcn;
fc2 = fc/fcn;
[b1, a1] = butter(4,[fc1 fc2]);
SSB_practical_time = filter(b1, a1, modulated_DSB_SC_time);
SSB_practical_Spectrum = real(fftshift(fft(SSB_practical_time)));
plot_time(SSB_practical_time,FS,'Butter-worth SSB-SC modulated Signal in Time Domain (LSB) 4th order');
plot_frequency(abs(SSB_practical_Spectrum),FS,'Butter-worth SSB-SC modulated Signal in Frequency Domain (LSB) 4th order');
%%% step (6) %%%
demodulated_time = SSB_practical_time .*carrierSignal;
demodulated_spectrum = fft(fftshift(demodulated_time));
[b2, a2] = butter(4,4000/fcn);
coherent_demodulated_time = filter(b2, a2, demodulated_time);
coherent_demodulated_Spectrum = fftshift(fft(coherent_demodulated_time));
plot_time(coherent_demodulated_time,FS,'Butter-worth Demodulated Signal in Time Domain (Coherent Detection)');
plot_frequency(abs(coherent_demodulated_Spectrum),FS,'Butter-worth Demodulated Signal in Frequency Domain (Coherent Detection)');
resampled_time = resample(coherent_demodulated_time, fs, FS);
resampled_freq = fftshift(fft(resampled_time));
plot_time(resampled_time,fs,'Butter-worth Demodulated Signal in Time Domain (Coherent Detection)');
plot_frequency(abs(real(resampled_freq)),fs,'Butter-worth Demodulated Signal in Frequency Domain (Coherent Detection)');
% sound(real(resampled_time),fs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%(8)%%%%%%%%%%%%%%%%%%%%%%%%%
for i = [0 10 30]
demodulated_coherent_SSB_SC = awgn(modulated_SSB_time,i);
time = linspace(0, length(demodulated_coherent_SSB_SC)/FS, length(demodulated_coherent_SSB_SC));
demodulated_signal_noise_time = demodulated_coherent_SSB_SC.*cos(2*pi*fc*time);
filter4 = designfilt('lowpassfir', 'FilterOrder', 8000, 'CutoffFrequency', 4000, 'SampleRate', FS);
demodulatedTime  = filter(filter4, demodulated_signal_noise_time);
demodulatedSpectrum = fftshift(fft(demodulatedTime));
plot_time(real(demodulatedTime),FS,['Demodulated SSB-SC signal in Time Domain:',num2str(i),' db']);
plot_frequency(abs(demodulatedSpectrum),FS,['Demodulated SSB-SC signal in Frequency Domain:',num2str(i),' db']);
resampled_demodulatedTime = resample(real(demodulatedTime),fs,FS);
resampled_demodulatedSpectrum = fftshift(fft(resampled_demodulatedTime));
plot_time(resampled_demodulatedTime,fs,['Demodulated SSB-SC signal in Time Domain with noise: ',num2str(i),' db']);
plot_frequency(real(abs(resampled_demodulatedSpectrum)), fs, ['Demodulated SSB-SC signal in Frequency Domain with noise: ',num2str(i),' db']);
audiowrite(['Signal_SSB_SC_SNR',num2str(i),'.wav'],resampled_demodulatedTime,fs);
%sound(resampled_demodulatedTime,fs);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%(9)%%%%%%%%%%%%%%%%%%%%%%%%%
Ac = 2 * max(abs(resampledSignal));
modulated_DSB_TC_time = (Ac + resampledSignal).* carrierSignal ;
filter5 = designfilt('lowpassfir', 'FilterOrder', 8000, 'CutoffFrequency', fc, 'SampleRate', FS);
modulated_SSB_TC_Time  = filter(filter5, modulated_DSB_TC_time);
envelope_SSB_TC = abs(hilbert(modulated_SSB_TC_Time));
plot(envelope_SSB_TC);
title('Envelope of SSB-TC signal');
grid;
envelopeResampled_SSB_TC_time = resample(envelope_SSB_TC,fs,FS);
plot_time(real(envelopeResampled_SSB_TC_time), fs, 'Demodulated SSB-TC signal using Envelope Detection in Time Domain');
% sound(envelopeResampled_SSB_TC,fs);




