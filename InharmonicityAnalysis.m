clear;
close;
clc;

% This is the final file for the analysis of the piano inharmonicity. All
% the data are initialized to the C2 of the Steinway B2 piano: to change it
% just open the correct sample from B2 or U3 and set the fundamental
% frequency (looking at the arrays in the section below you can find the
% fundamentals for both pianos). 

% Section I implements the algorithm by Vesa Valimaki (cit. in Readme) and
% makes the plots for data analysis. 
% Section II builds a synthetic spectrum from ideal and real partials and
% allows to listen to the differences in a pure sinusoidal spectrum. To run
% this section Audio Toolbox is required.
%% Steinway B2 sample C2

[y,Fs] = audioread("./SteinwayB2samples/Piano.mf.C2.aiff");
%[y,Fs] = audioread("./YamahaU3samples/C2.wav");
y_mono = sum(y,2)/size(y,2);
Nfft = 2^17;

y_fft = fft(y_mono,Nfft);
y_fft(1:ceil(20*Nfft/Fs))=0;
y_fft = abs(y_fft(1:Nfft/2));
y_fft = y_fft/max(y_fft);

f = Fs/2 * linspace(0,1,Nfft/2);

% Parameters initialization
B = 0.0001;
delta = 1;
counter = 0;
old_trend = 1;
f1 = 65.1;
deltaF = 0.4 * f1;

% DATA TO CYCLE OVER NOTES
B2fundamentals = [32.323 65.1 131.1 262.1 524.9];
U3fundamentals = [31.73 64.94 130.5 261.4 523.9];

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

% Theoretical values of partials for computed inharmonicity coefficient, based on eq. (9)
f_theoretical = (1:25) *f1/(1+B)^(-0.5) .* ((1+B.*(1:25).^2).^(0.5));
f_ideal = (1:25) * f1;
% Linear interoplation of spectral peaks
c = polyfit((1:25).^2, (f(f_peaks) ./(1:25)).^2, 1);
y_est = polyval(c, (1:25).^2);


% Plots
figure();
plot(f,abs(y_fft));
title("Spectrum Steinway B2: C2");
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
title("Ideal vs Real Partials Steinway B2: C2");
xlabel("f [Hz]");
ylabel("FFT");
xlim([0, ceil(f_peaks(15)/Nfft * Fs /10)*10]);
plot(f(f_peaks), peaks, 'or');
stem(f_ideal, peaks);
legend('FFT','Real partials','Ideal partials');


figure();
plot((1:25).^2, y_est);
hold on;
plot((1:25).^2, (f(f_peaks)./(1:25)).^2, 'or');
title("Inharmonicity of first 25 partials Steinway B2: C2");
ylabel("(f_{n}/n)^2");
xlabel("n^2");

figure();
M=floor(0.050*Fs);
R=floor(M*0.5);
w=window(@bartlett,M);
N=4*4096;
spectrogram(y_mono,w,R, N, Fs, 'yaxis');
ylim([0,3]);
title("STFT Steinway B2");

%% Play the sound
oscillators = {};
normalizedAmp = peaks ./ sum(peaks);

deviceWriter = audioDeviceWriter(44100);
deviceWriter.SupportVariableSizeInput = true;
deviceWriter.BufferSize = 64;

%% IDEAL SPECTRUM
for i=1:length(peaks)
    osc = audioOscillator('sine','Frequency', f_ideal(i), 'Amplitude', peaks(i), 'SamplesPerFrame', 44100);
    oscillators(i) = {osc};
end

%Create oscillators for first 15 harmonics
osc1 = oscillators{1};
osc2 = oscillators{2};
osc3 = oscillators{3};
osc4 = oscillators{4};
osc5 = oscillators{5};
   
osc6 = oscillators{6};
osc7 = oscillators{7};
osc8 = oscillators{8};
osc9 = oscillators{9};
osc10 = oscillators{10};
    
osc11 = oscillators{11};
osc12 = oscillators{12};
osc13 = oscillators{13};
osc14 = oscillators{14};
osc15 = oscillators{15};

disp("Press enter to hear ideal spectrum, first 15 harmonics");
pause();
tic
while toc<5
    %deviceWriter(waveform);
    deviceWriter(osc1() + osc2() + osc3() + osc4() + osc5() ...
         + osc6() + osc7() + osc8() + osc9() + osc10() ...
        + osc11() + osc12() + osc13() + osc14() + osc15());
end

%% REAL SPECTRUM

% Uncomment to set manually the B coefficient and hear differences
% B = 0.0001;
% f_theoretical = (1:25) *f1/(1+B)^(-0.5) .* ((1+B.*(1:25).^2).^(0.5));
for i=1:length(peaks)
    osc = audioOscillator('sine','Frequency', f_theoretical(i), 'Amplitude', peaks(i), 'SamplesPerFrame', 44100);
    oscillators(i) = {osc};
end

%Create oscillators for first 15 harmonics
osc1 = oscillators{1};
osc2 = oscillators{2};
osc3 = oscillators{3};
osc4 = oscillators{4};
osc5 = oscillators{5};
   
osc6 = oscillators{6};
osc7 = oscillators{7};
osc8 = oscillators{8};
osc9 = oscillators{9};
osc10 = oscillators{10};
    
osc11 = oscillators{11};
osc12 = oscillators{12};
osc13 = oscillators{13};
osc14 = oscillators{14};
osc15 = oscillators{15};

disp("Press enter to hear real spectrum, first 15 harmonics");
pause();
tic
while toc<5
    deviceWriter(osc1() + osc2() + osc3() + osc4() + osc5() ...
         + osc6() + osc7() + osc8() + osc9() + osc10() ...
        + osc11() + osc12() + osc13() + osc14() + osc15());
end