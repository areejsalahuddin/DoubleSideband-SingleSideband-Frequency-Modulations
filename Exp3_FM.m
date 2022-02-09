%%%%%%%%%question1%%%%%%%%%%

[y , fs] = audioread('eric.wav');
% sound(y,fs);
plot_time(y,fs,'Signal time domain');
y_Spectrum = (fftshift(fft(y)));
plot_frequency(abs(y_Spectrum),fs,'Signal Spectrum');
f = linspace(-fs/2,fs/2,length(y)); 

filter = generate_filter(length(y_Spectrum),fs);
filteredSpectrum = y_Spectrum.* filter ;
ytime_filtered= real(ifft(ifftshift(filteredSpectrum)));
plot_frequency(real(filter),fs,'Filter');
plot_time(ytime_filtered,fs,'Filtered Signal in Time Domain');
plot_frequency(abs(filteredSpectrum),fs,'Filtered Signal in Frequency Domain');
%sound(ytime_filtered,fs);

%%%%%%%%%question2%%%%%%%%%%

fc = 100000;
FS = 5*fc;
resampledSignal = resample(ytime_filtered,FS,fs);
t = linspace(0,length(resampledSignal)/FS,length(resampledSignal));
A = max(abs(resampledSignal));
Kf = 0.2/(2*pi* max(cumsum(resampledSignal)*(1/FS)));
delta = Kf.*cumsum(resampledSignal).';
NBFM_time = A.*cos(2* pi*fc*t  + delta);
plot_time(NBFM_time, FS, 'NBFM modulated signal in Time Domain');
NBFM_frequency=abs(fftshift(fft(NBFM_time)));
plot_frequency(NBFM_frequency, FS, 'NBFM modulated signal in Frequency Domain');



%%%%%%%%%question3%%%%%%%%%%
%beta<<1


%%%%%%%%%question4%%%%%%%%%%

signalDiff_AM = diff(NBFM_time);
signal_envelope=abs(hilbert(signalDiff_AM)); 
signal_envelopeResampled=resample(signal_envelope,fs,FS);
plot_time(signal_envelopeResampled, fs, 'NBFM demodulated signal in Time Domain using envelope');
signal_envelopeResampledFreq = abs(fftshift(fft(signal_envelopeResampled)));
plot_frequency(signal_envelopeResampledFreq, fs, 'NBFM demodulated signal in Frequency Domain using envelope');
% sound(signal_envelopeResampled,fs);











