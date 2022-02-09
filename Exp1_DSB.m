%%%%%%%%%%%%%%% question 1%%%%%%%%%%%%%%

% filename = input('Enter file name: ','s');
% file_name = [filename '.wav'];
[y , fs] = audioread('eric.wav');
% sound(y,fs);
plot_time(y,fs,'Signal time domain');
Yspectrum = (fftshift(fft(y)));
f = linspace(-fs/2,fs/2,length(y)); 
plot_frequency(abs(Yspectrum),fs,'Signal Spectrum');


%%%%%%%%%%%%%%% question 2%%%%%%%%%%%%%%

 
filter = generate_filter(length(Yspectrum),fs);
plot_frequency(filter,fs,'filter');
filteredSpectrum = Yspectrum .* filter ;
ytime_filtered= real(ifft(ifftshift(filteredSpectrum))); % for sounding


%%%%%%%%%%%%%%% question 3%%%%%%%%%%%%%%

plot_time(ytime_filtered,fs,'Filtered Signal in Time Domain');
plot_frequency(abs(filteredSpectrum),fs,'Filtered Signal in Frequency Domain');


%%%%%%%%%%%%%%% question 4%%%%%%%%%%%%%%

% sound(ytime_filtered,fs);

%%%%%%%%%%%%%%% question 5%%%%%%%%%%%%%%

fc = 100000;
% [p,q] = rat(5*fc/fs);
FS = 5*fc;
resampledSignal = resample(ytime_filtered,FS,fs);    
AC = 2 * max(abs(resampledSignal));
grid;

t = linspace(0,length(resampledSignal)/FS,length(resampledSignal));

%DSB-SC 
carrierSignal = cos(2 * pi * fc * t);
modulated_time_DSB_SC= carrierSignal.'.* resampledSignal;
modulated_DSB_SC = fftshift(fft(modulated_time_DSB_SC));
f_new = linspace(-FS/2, FS/2,length(resampledSignal)); 
figure;
plot(f_new,abs(modulated_DSB_SC));
title('DSB-SC modulated signal in Frequency domain');
grid;

%DSB_TC
modulated_time_DSB_TC= (AC + resampledSignal).* carrierSignal.' ;
modulated_DSB_TC = fftshift(fft(modulated_time_DSB_TC));
figure;
plot(f_new,abs(modulated_DSB_TC));
title('DSB-TC modulated signal in Frequency domain');
grid;

plot_time(modulated_time_DSB_SC,FS,'DSB-SC modulated signal in time domain');
plot_time(modulated_time_DSB_TC,FS,'DSB-TC modulated signal in time domain');

%%%%%%%%%%%%%%% question 6%%%%%%%%%%%%%%

%DSB_SC 
envelope_SC = abs(hilbert(modulated_time_DSB_SC));
%DSB-TC
envelope_TC = abs(hilbert(modulated_time_DSB_TC));

%%%%%%%%%%%%%%% question 7%%%%%%%%%%%%%%

%DSB-SC
figure;
plot(envelope_SC);
title('Envelope of DSB-SC signal');
grid;
envelopeResampled_SC=resample(envelope_SC,fs,FS);
%sound(envelopeResampled_SC,fs);

%DSB-TC
figure;
plot(envelope_TC);
title('Envelope of DSB-TC signal');
grid;
envelopeResampled_TC=resample(envelope_TC,fs,FS);
% sound(envelopeResampled_TC,fs); 
% envelope detector should be used with TC only 
 
%%%%%%%%%%%%%%% question 8%%%%%%%%%%%%%%

% DSB-SC
for i = [0 10 30]
demodCoherent_SC = awgn(modulated_time_DSB_SC,i); 
t = linspace(0,length(demodCoherent_SC)/FS,length(demodCoherent_SC));
demodCoherent_SC= demodCoherent_SC.*cos(2*pi*fc*t).';  %Coherent detection
Yr = fftshift(fft(demodCoherent_SC));  %Fourier transform
Yr_filter= generate_filter(length(modulated_time_DSB_SC),FS);
Yr_msg = Yr.*Yr_filter;  %Obtained audio message after coherent detection
plot_frequency(abs(Yr_msg),FS,'DSB-SC audio signal spectrum');
time_demodSC = ifft(ifftshift(Yr_msg));   %Inverse fourier transform
time_demodSC = resample(real(time_demodSC),fs,FS);
plot_time(time_demodSC,fs,['DSB-SC audio signal in time domain with noise: ',num2str(i),' db']);
freq_demodSC = fftshift(fft(time_demodSC));
plot_frequency(real(abs(freq_demodSC)), fs, ['DSB-SC audio signal in freq domain with noise: ',num2str(i),' db']);

% audiowrite(['Signal_DSB_SC_SNR',num2str(i),'.wav'],time_demodSC,fs);
end

%%%%%%%%%%%%%%% question 9%%%%%%%%%%%%%% 
fc_new = 100100;

t = linspace(0,length(modulated_time_DSB_SC)/FS,length(modulated_time_DSB_SC));
demod_freqError= modulated_time_DSB_SC.*cos(2*pi*fc_new*t).';  %Coherent detection
demod_freqError = fftshift(fft(demod_freqError)); 
demod_freqError = demod_freqError.*Yr_filter;  
plot_frequency(abs(demod_freqError),FS,'DSB-SC audio signal with frequency error spectrum');
time_demodSC = ifft(ifftshift(demod_freqError));   %Inverse fourier transform
time_demodSC = resample(real(time_demodSC),fs,FS);
plot_time(time_demodSC,fs,'DSB-SC audio signal in time domain with freq error');
freq_demodSC = fftshift(fft(time_demodSC));
plot_frequency(real(abs(freq_demodSC)), fs, 'DSB-SC audio signal in freq domain with freq error');
%sound(time_demodSC,fs);
% audiowrite(['Signal_with_frequencyError.wav'],time_demodSC,fs);

%%%%%%%%%%%%%%% question 10%%%%%%%%%%%%%%
phaseshift = degtorad(20);
t = linspace(0,length(modulated_time_DSB_SC)/FS,length(modulated_time_DSB_SC));
demod_phaseError= modulated_time_DSB_SC.*cos(2*pi*fc*t + phaseshift).';  %Coherent detection
demod_phaseError = fftshift(fft(demod_phaseError)); 
demod_phaseError = demod_phaseError.*Yr_filter;  
plot_frequency(abs(demod_phaseError),FS,'DSB-SC audio signal with phase error spectrum');
time_demodSC = ifft(ifftshift(demod_phaseError));   %Inverse fourier transform
time_demodSC = resample(real(time_demodSC),fs,FS);
plot_time(time_demodSC,fs,'DSB-SC audio signal in time domain with phase error');
freq_demodSC = fftshift(fft(time_demodSC));
plot_frequency(real(abs(freq_demodSC)), fs, 'DSB-SC audio signal in freq domain with phase error');
%sound(time_demodSC,fs);
% audiowrite(['Signal_with_frequencyError.wav'],time_demodSC,fs);
