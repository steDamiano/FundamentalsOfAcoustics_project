clear;
close;
clc;

%% Steinway B2 sample C1

[y,Fs] = audioread("./SteinwayB2samples/Piano.mf.C1.aiff");
%[y,Fs] = audioread("./YamahaU3samples/C3.wav");
y_mono = sum(y,2)/size(y,2);
Nfft = 2^17;

y_fft = fft(y_mono,Nfft);
y_fft(1:ceil(20*Nfft/Fs))=0;
y_fft=abs(y_fft(1:Nfft/2));

f = Fs/2 * linspace(0,1,Nfft/2);

% Parameters initialization
B = 0.0001;
delta = 1;
counter = 0;
old_trend = 1;
f1 = 32.323;
deltaF = 0.4 * f1;

%% Peak reduction step
% deltaFreq = 5*f1;
% [peaks, locs] = findpeaks(y_fft(1:floor(deltaFreq*Nfft/Fs)));
% 
% [MAXpeaks, MAXidx] = maxk(peaks,10);
% MAXlocs = locs(MAXidx);
% 
% figure();
% plot(f, y_fft);
% xlim([0,1000]);
% hold on;
% plot(f(MAXlocs),MAXpeaks, 'o');

% Try to cycle through frequencies
bin = 1;
deltaFreq = 5*f1;
deltaBin = floor(deltaFreq * Nfft/Fs);
reducedPeaks = [];
reducedFreqs = [];
while(bin < length(y_fft) - deltaBin)
    [peaks, locs] = findpeaks(y_fft(bin:bin + deltaBin));

    [MAXpeaks, MAXidx] = maxk(peaks,10);
    MAXlocs = locs(MAXidx) + bin;
    
    reducedPeaks = [reducedPeaks; MAXpeaks];
    reducedFreqs = [reducedFreqs; MAXlocs];
    
    bin = bin + deltaBin;
end

disp("Loop end");

% % Plot to check selected peaks
% figure();
% plot(f, y_fft);
% xlim([0,1000]);
% hold on;
% plot(f(reducedFreqs),reducedPeaks, 'o');

%% Iteration Loop
peaks = zeros(1,25);
f_peaks = zeros(1,25);
while(counter <100 && abs(delta) > 10e-4)
    fk_estimation = (1:25) .* f1 / (1+B)^(-0.5) .* (1+B.*(1:25).^2).^0.5;
    
    for i=1:length(fk_estimation)
    % set interval for peak estimation
    lowLimit = ceil(fk_estimation(i) * Nfft/Fs - deltaF);
    upLimit = ceil(fk_estimation(i)*Nfft/Fs + deltaF);
    % find highest peak in selected interval
    peaks(i) = max(y_fft(lowLimit:upLimit));
    f_peaks(i) = find(y_fft == peaks(i));
    end
    
    % compute error and trend
    Dk = f(f_peaks) - fk_estimation;
    trend = sum(sign(Dk(2:end) - Dk(1:end-1)));
    
    %evaluate delta
    if(trend >0)
            delta = delta * sign(delta);
    end
    if(trend<0)
        delta = delta * sign(delta);
        delta = delta * sign(trend);
    end
    if(sign(trend) ~= sign(old_trend))
        delta = delta/2;
    end
    
    %update parameters
    old_trend = trend;
    B = B*10^delta;
    counter= counter+1;
end

% Theoretical values for inharmonicity coefficient, based on eq. (1)
f_theoretical = (1:25) *f1/(1+B)^(-0.5) .* ((1+B.*(1:25).^2).^(0.5));
f_ideal = (1:25) * f1;
% Linear interoplation of spectral peaks
c = polyfit((1:25).^2, (f(f_peaks) ./(1:25)).^2, 1);
y_est = polyval(c, (1:25).^2);


% Plots
figure();
plot(f,abs(y_fft));
title("Spectrum Steinway B2: C1");
xlabel("f [Hz]");
ylabel("FFT");
xlim([0, ceil(f_peaks(length(f_peaks))/Nfft * Fs /10)*10]);
hold on;
plot(f(f_peaks), peaks, 'or');
plot(f_theoretical, peaks, '.b');
legend('FFT','Spectrum Peaks', 'Computed Peaks');

figure()
plot(f, y_fft);
hold on;
title("Ideal vs Real Partials");
xlabel("f [Hz]");
ylabel("FFT");
xlim([0, ceil(f_peaks(15)/Nfft * Fs /10)*10]);
plot(f(f_peaks), peaks, 'or');
stem(f_ideal, peaks);
legend('FFT','Inharmonic partials','Ideal partials');


figure();
plot((1:25).^2, y_est);
hold on;
plot((1:25).^2, (f(f_peaks)./(1:25)).^2, 'or');
title("Inharmonicity of first 25 partials Steinway B2");
ylabel("(f_{n}/n)^2");
xlabel("n^2");

figure();
M=floor(0.050*Fs);
R=floor(M*0.5);
w=window(@bartlett,M);
N=4096;
spectrogram(y_mono,w,R, N, Fs, 'yaxis');
ylim([0,2]);
title("STFT Steinway B2");
