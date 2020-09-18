clear;
close all;
clc;

%% PIANO C1
%Import audio signal
[y, Fs] = audioread("./SteinwayB2samples/Piano.mf.C1.aiff");
y_mono = sum(y,2) / size(y,2);
%y_mono = y_mono(length(y_mono)/32:length(y_mono)/32 + 50000);

%FFT of signal
Nfft = 2^nextpow2(length(y_mono));
y_fft = fft(y_mono, Nfft);
y_fft = y_fft(1:Nfft/2);
y_fft(1:ceil(20*Nfft/Fs)) = 0;
f = Fs /2 * linspace(0,1,Nfft/2);

%peaks in the sampled note
[peaks, bins] = findpeaks(abs(y_fft), 'MinPeakDistance', 30 * Nfft / Fs, 'MinPeakHeight', 3);

%plot fft and peaks found from audio 
figure();
plot(f(1:ceil(1050*Nfft / Fs)), abs(y_fft(1 : ceil(1050 * Nfft/Fs))));
hold on;
stem(f(bins(1:15)), peaks(1:15), 'or');


%peaks computed with formula f = n* f0 -> f0 = 65
fund = 32;
ideal_harmonics = zeros(1,15);
for i=1:15
    ideal_harmonics(i) = i * fund;
end

%plot of ideal harmonics
stem(f(ideal_harmonics) * Nfft / Fs + 1, peaks(1:15), 'og')
title("C1 Grand Piano Steinway Model B");
ylabel("FFT");
xlabel("f [Hz]");
hold off;


% spectrogram
N=length(y_mono);
M=floor(0.050*Fs);
R=floor(M*0.5);
w=window(@bartlett,M);
Nfft=4096;
% MATLAB WAY
figure();
spectrogram(y_mono,w,R, Nfft, Fs, 'yaxis');
ylim([0,2]);

% Analysis of fn / n (behaviour of partials wrt to ideal case)
figure();
plot((1:10).^2, (f(bins(1:10)) ./ (1:10)).^2, 'or');
hold on;
plot((1:10).^2, (f(bins(1))*ones(10)).^2, '--b');
title("Inharmonicity of partials C1");
ylabel("f_{n} / n");
xlabel("n partial")

c = polyfit((1:10).^2, (f(bins(1:10)) ./(1:10)).^2, 1);
y_est = polyval(c, (1:10).^2);
plot((1:10).^2, y_est);
%% PIANO C2
%Import audio signal
[y, Fs] = audioread("./SteinwayB2samples/Piano.mf.C2.aiff");
y_mono = sum(y,2) / size(y,2);

%FFT of signal
Nfft = 2^nextpow2(length(y_mono));
y_fft = fft(y_mono, Nfft);
y_fft = y_fft(1:Nfft/2);
y_fft(1:20) = 0;
f = Fs /2 * linspace(0,1,Nfft/2);

%peaks in the sampled note
[peaks, bins] = findpeaks(abs(y_fft), 'MinPeakDistance', 60 * Nfft / Fs, 'MinPeakHeight', 20);

%plot fft and peaks found from audio 
plot(f(1:ceil(1050*Nfft / Fs)), abs(y_fft(1 : ceil(1050 * Nfft/Fs))));
hold on;
stem(f(bins(1:10)), peaks(1:10), 'or');


%peaks computed with formula f = n* f0 -> f0 = 65
fund = 65;
ideal_harmonics = zeros(1,10);
for i=1:10
    ideal_harmonics(i) = i * fund;
end

%plot of ideal harmonics
stem(f(ideal_harmonics) * Nfft / Fs, peaks(1:10), 'og')
title("C2 Grand Piano Steinway Model B");
ylabel("FFT");
xlabel("f [Hz]");
hold off;

%% PIANO C3
%Import audio signal
[y, Fs] = audioread("Piano.mf.C3.aiff");
y_mono = sum(y,2) / size(y,2);

%FFT of signal
Nfft = 2^nextpow2(length(y_mono));
y_fft = fft(y_mono, Nfft);
y_fft = y_fft(1:Nfft/2);
y_fft(1:20) = 0;
f = Fs /2 * linspace(0,1,Nfft/2);

%peaks in the sampled note
[peaks, bins] = findpeaks(abs(y_fft), 'MinPeakDistance', 60 * Nfft / Fs, 'MinPeakHeight', 20);

%plot fft and peaks found from audio 
plot(f(1:ceil(2050*Nfft / Fs)), abs(y_fft(1 : ceil(2050 * Nfft/Fs))));
hold on;
stem(f(bins(1:10)), peaks(1:10), 'or');


%peaks computed with formula f = n* f0 -> f0 = 65
fund = 130;
ideal_harmonics = zeros(1,10);
for i=1:10
    ideal_harmonics(i) = i * fund;
end

%plot of ideal harmonics
stem(f(ideal_harmonics) * Nfft / Fs, peaks(1:10), 'og')
title("C3 Grand Piano Steinway Model B");
ylabel("FFT");
xlabel("f [Hz]");
hold off;

%% PIANO C4
%Import audio signal
[y, Fs] = audioread("Piano.mf.C4.aiff");
y_mono = sum(y,2) / size(y,2);

%FFT of signal
%Nfft = 2^nextpow2(length(y_mono));
Nfft = 2^16;
y_fft = fft(y_mono, Nfft);
y_fft = y_fft(1:Nfft/2);
y_fft(1:ceil(20*Nfft/Fs)) = 0;
y_fft(ceil(4800*Nfft/Fs):ceil(4900*Nfft/Fs)) = 0;

f = Fs /2 * linspace(0,1,Nfft/2);

%peaks in the sampled note
[peaks, bins] = findpeaks(abs(y_fft), 'MinPeakDistance', 262 * Nfft / Fs, 'MinPeakHeight', 0.01);
peaks = peaks(1:20);
bins = bins(1:20);
%plot fft and peaks found from audio 
plot(f, abs(y_fft));
hold on;
stem(f(bins), peaks, 'or');


%peaks computed with formula f = n* f0 -> f0 = 65
fund = 262;
ideal_harmonics = zeros(1,20);
for i=1:20
    ideal_harmonics(i) = i * fund;
end

%plot of ideal harmonics
%stem(f(ideal_harmonics) * Nfft / Fs, peaks, 'og')
title("C4 Grand Piano Steinway Model B");
ylabel("FFT");
xlabel("f [Hz]");
hold off;

figure();
plot((1:20).^2, (f(bins(1:20)) ./ (1:20)).^2, 'or');
hold on;
plot((1:20).^2, (f(bins(1))*ones(20)).^2, '--b');
title("Inharmonicity of partials C1");
ylabel("f_{n} / n");
xlabel("n partial")

c = polyfit((1:20).^2, (f(bins(1:20)) ./(1:20)).^2, 1);
y_est = polyval(c, (1:20).^2);
plot((1:20).^2, y_est);


%% Experimental Setup:
%{
    2 cardioid microphones t.bone sc140 -> PIANO WITH FRONT PANEL
    Left Mic: A2, 5cm over top of piano wood
    Right Mic: A5, 5cm top of piano wood
    Recording stereo tracks through logic, cut at the attack of note. 
    Fs = 44100;
%}

%% Upright piano notes

[y,Fs] = audioread("./YamahaU3samples/C1.wav");
y = sum(y,2) / size(y,2);
y=y(1:length(y));
Nfft = 2^nextpow2(length(y));

y_fft = fft(y, Nfft);
y_fft = y_fft(1:Nfft/2);
f = Fs/2 * linspace(0, 1,Nfft/2);
[peaks, bins] = findpeaks(abs(y_fft), 'MinPeakDistance', 30 * Nfft / Fs, 'MinPeakHeight', 30);

% FFT plot with partials 
figure();
plot(f(1:ceil(1050*Nfft/Fs)),abs(y_fft(1:ceil(1050*Nfft/Fs))));
hold on;
%stem(f(bins(1:15)), peaks(1:15), 'or');

title("Yamaha U3 Upright Piano");
xlabel("f [Hz]");
ylabel("FFT");
%peaks computed with formula f = n* f0 -> f0 = 65
fund = 32;
ideal_harmonics = zeros(1,15);
for i=1:15
    ideal_harmonics(i) = i * fund;
end
%stem(f(ideal_harmonics)* Nfft / Fs, peaks(1:15), 'og');
hold off;

% Spectrogram

% MANUAL IMPLEMENTATION
N=length(y);
M=floor(0.050*Fs);
R=floor(M*0.5);
w=window(@bartlett,M);
Nfft=4096;
y(N+M)=0;
X=zeros(floor(Nfft/2+1),floor(N/R)+1);
for m=1:floor(N/R)+1
  mR=(m-1)*R+1; x_m=y(mR:mR+M-1);
  x_wm=x_m.*w; X_wm=fft(x_wm, Nfft);
  X(:,m)=X_wm(1:floor(Nfft/2+1));
end
% X=flipud(X);
time=(1:size(X,2))*0.025;
freq=linspace(Fs/2,0,size(X,1));
figure();
%imagesc(time,freq,20*log10(abs(X)));

% MATLAB WAY
spectrogram(y,w,R, Nfft, Fs, 'yaxis');
ylim([0,2]);
title("STFT C1 Yamaha U3");

%Analysis of fn / n (behaviour of partials wrt to ideal case)
figure();
plot(1:10, f(bins(1:10)) ./ (1:10), 'or');
hold on;
plot(1:10, f(bins(1))*ones(10), '--b');
title("Inharmonicity of partials C1");
ylabel("f_{n} / n");
xlabel("n partial")

figure();
plot((1:10).^2, (f(bins(1:10)) ./ (1:10)).^2, 'or');
hold on;
plot((1:10).^2, (f(bins(1))*ones(10)).^2, '--b');
title("Inharmonicity of partials C1");
ylabel("f_{n} / n");
xlabel("n partial")

c = polyfit((1:10).^2, (f(bins(1:10)) ./(1:10)).^2, 1);
y_est = polyval(c, (1:10).^2);
plot((1:10).^2, y_est);