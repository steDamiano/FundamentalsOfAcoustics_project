clear;
close;
clc;

% This file is a study of the implementation of the Fast Algorithm for
% B coefficient estimation by Vesa Valimaki (reference in Readme). 
%% Steinway B2 sample C1

[y,Fs] = audioread("./SteinwayB2samples/Piano.mf.C1.aiff");
%[y,Fs] = audioread("./YamahaU3samples/C1.wav");
y_mono = sum(y,2)/size(y,2);
%y_mono = y_mono(1:44100);
Nfft = 2^17;

y_fft = fft(y_mono,Nfft);
y_fft(1:ceil(20*Nfft/Fs))=0;

y_fft=abs(y_fft(1:Nfft/2));
f = Fs/2 * linspace(0,1,Nfft/2);
f1 = 32.323;
[p, idx] = findpeaks(y_fft, 'MinPeakHeight', 1.7, 'MinPeakDistance', 31.5 * Nfft / Fs);

figure();
plot(f,abs(y_fft));
hold on; 
p=p(1:25);
idx = idx(1:25);
stem(f(idx),p);
title("Spectrum and Partials Steinway B2: C1");
xlim([0 2000]);

%starting value for B coefficient
B = 10e-4;
delta = 1;
counter = 0;
old_trend = 1;

while(counter <100 && abs(delta) > 10e-4)
    fk_estimation = (1:25) .* f1 / (1+B)^(-0.5) .* (1+B.*(1:25).^2).^0.5;
    Dk = f(idx) - fk_estimation;
    trend = sum(sign(Dk(2:end) - Dk(1:end-1)));
    %delta = delta * sign(trend);
    if(trend >0)
            delta = delta * sign(delta);
    end
    if(trend<0)
        if(delta < 0)
            delta = delta * sign(delta);
        end
        if(delta > 0)
            delta = delta * sign(trend);
        end
    end
    if(sign(trend) ~= sign(old_trend))
        delta = delta/2;
    end
    old_trend = trend;
    B = B*10^delta;
    counter= counter+1;
end
% Theoretical values of harmonics, using f1 as fundamental
f_theoretical = (1:25) *f1/(1+B)^(-0.5) .* ((1+B.*(1:25).^2).^(0.5));

c = polyfit((1:25).^2, (f(idx(1:25)) ./(1:25)).^2, 1);
y_est = polyval(c, (1:25).^2);

figure();
plot((1:25).^2, y_est);
hold on;
plot((1:25).^2, (f(idx)./(1:25)).^2, 'or');
title("Inharmonicity of first 30 partials Steinway B2: C1");
ylabel("(f_{n}/n)^2");
xlabel("n^2");


%% Yamaha U3 sample C1

[y,Fs] = audioread("./YamahaU3samples/C1.wav");
y_mono = sum(y,2)/size(y,2);
Nfft = 2^17;

y_fft = fft(y_mono,Nfft);
y_fft(1:ceil(20*Nfft/Fs))=0;
y_fft=abs(y_fft(1:Nfft/2));
f = Fs/2 * linspace(0,1,Nfft/2);
f1 = 31.73;

[p, idx] = findpeaks(y_fft, 'MinPeakHeight', 18, 'MinPeakDistance', 31.5 * Nfft / Fs);
figure();
plot(f,abs(y_fft));
hold on; 
p=p(1:25);
idx = idx(1:25);
stem(f(idx),p);
title("Spectrum and Partials Yamaha U3: C1");
xlim([0 2000]);

%starting value for B coefficient
B = 10e-4;
delta = 1;
counter = 0;
old_trend = 1;

while(counter <100 && abs(delta) > 10e-4)
    fk_estimation = (1:25) .* f1 / (1+B)^(-0.5) .* (1+B.*(1:25).^2).^0.5;
    Dk = f(idx) - fk_estimation;
    trend = sum(sign(Dk(2:end) - Dk(1:end-1)));
    %delta = delta * sign(trend);
    if(trend >0)
            delta = delta * sign(delta);
    end
    if(trend<0)
        if(delta < 0)
            delta = delta * sign(delta);
        end
        if(delta > 0)
            delta = delta * sign(trend);
        end
    end
    if(sign(trend) ~= sign(old_trend))
        delta = delta/2;
    end
    old_trend = trend;
    B = B*10^delta;
    counter= counter+1;
end

% Theoretical values of harmonics, using f1 as fundamental
f_theoretical = (1:25) *f1/(1+B)^(-0.5) .* ((1+B.*(1:25).^2).^(0.5));
c = polyfit((1:25).^2, (f(idx(1:25)) ./(1:25)).^2, 1);
y_est = polyval(c, (1:25).^2);

figure();
plot((1:25).^2, y_est);
hold on;
plot((1:25).^2, (f(idx)./(1:25)).^2, 'or');
title("Inharmonicity of first 30 partials Yamaha U3: C1");
ylabel("(f_{n}/n)^2");
xlabel("n^2");

%% Steinway B2 sample C2

[y,Fs] = audioread("./SteinwayB2samples/Piano.mf.C2.aiff");
y_mono = sum(y,2)/size(y,2);
Nfft = 2^17;

y_fft = fft(y_mono,Nfft);
y_fft(1:ceil(20*Nfft/Fs))=0;
y_fft(ceil(1715*Nfft/Fs):ceil(1720*Nfft/Fs)) = 0;

y_fft=abs(y_fft(1:Nfft/2));
f = Fs/2 * linspace(0,1,Nfft/2);
f1 = 65.2;

[p, idx] = findpeaks(y_fft, 'MinPeakHeight', 1, 'MinPeakDistance', 65 * Nfft / Fs);
figure();
plot(f,abs(y_fft));
hold on; 
p=p(1:25);
idx = idx(1:25);
stem(f(idx),p);
title("Spectrum and Partials Steinway B2: C2");
xlim([0 2000]);

%starting value for B coefficient
B = 10e-4;
delta = 1;
counter = 0;
old_trend = 1;
% Theoretical values of harmonics, using f1 as fundamental

while(counter <100 && abs(delta) > 10e-4)
    fk_estimation = (1:25) .* f1 / (1+B)^(-0.5) .* (1+B.*(1:25).^2).^0.5;
    Dk = f(idx) - fk_estimation;
    trend = sum(sign(Dk(2:end) - Dk(1:end-1)));
    %delta = delta * sign(trend);
    if(trend >0)
            delta = delta * sign(delta);
    end
    if(trend<0)
        if(delta < 0)
            delta = delta * sign(delta);
        end
        if(delta > 0)
            delta = delta * sign(trend);
        end
    end
    if(sign(trend) ~= sign(old_trend))
        delta = delta/2;
    end
    old_trend = trend;
    B = B*10^delta;
    counter= counter+1;
end
f_theoretical = (1:25) *f1/(1+B)^(-0.5) .* ((1+B.*(1:25).^2).^(0.5));

c = polyfit((1:25).^2, (f(idx(1:25)) ./(1:25)).^2, 1);
y_est = polyval(c, (1:25).^2);

figure();
plot((1:25).^2, y_est);
hold on;
plot((1:25).^2, (f(idx)./(1:25)).^2, 'or');
title("Inharmonicity of first 30 partials Steinway B2: C2");
ylabel("(f_{n}/n)^2");
xlabel("n^2");


%% Yamaha U3 sample C2

[y,Fs] = audioread("./YamahaU3samples/C2.wav");
y_mono = sum(y,2)/size(y,2);
Nfft = 2^17;

y_fft = fft(y_mono,Nfft);
y_fft(ceil(1710*Nfft/Fs):ceil(1715*Nfft/Fs)) = 0;
y_fft(1:ceil(20*Nfft/Fs))=0;

y_fft=abs(y_fft(1:Nfft/2));
f = Fs/2 * linspace(0,1,Nfft/2);
f1 = 64.94;
[p, idx] = findpeaks(y_fft, 'MinPeakHeight', 0.2, 'MinPeakDistance', 64 * Nfft / Fs);

figure();
plot(f,abs(y_fft));
hold on; 
p=p(1:25);
idx = idx(1:25);
stem(f(idx),p);
title("Spectrum and Partials Yamaha U3: C2");
xlim([0 2000]);

%starting value for B coefficient
B = 10e-4;
delta = 1;
counter = 0;
old_trend = 1;

while(counter <100 && abs(delta) > 10e-4)
    fk_estimation = (1:25) .* f1 / (1+B)^(-0.5) .* (1+B.*(1:25).^2).^0.5;
    Dk = f(idx) - fk_estimation;
    trend = sum(sign(Dk(2:end) - Dk(1:end-1)));
    %delta = delta * sign(trend);
    if(trend >0)
            delta = delta * sign(delta);
    end
    if(trend<0)
        if(delta < 0)
            delta = delta * sign(delta);
        end
        if(delta > 0)
            delta = delta * sign(trend);
        end
    end
    if(sign(trend) ~= sign(old_trend))
        delta = delta/2;
    end
    old_trend = trend;
    B = B*10^delta;
    counter= counter+1;
end

% Theoretical values of harmonics, using f1 as fundamental
f_theoretical = (1:25) *f1/(1+B)^(-0.5) .* ((1+B.*(1:25).^2).^(0.5));
c = polyfit((1:25).^2, (f(idx(1:25)) ./(1:25)).^2, 1);
y_est = polyval(c, (1:25).^2);

figure();
plot((1:25).^2, y_est);
hold on;
plot((1:25).^2, (f(idx)./(1:25)).^2, 'or');
title("Inharmonicity of first 30 partials Yamaha U3: C2");
ylabel("(f_{n}/n)^2");
xlabel("n^2");

%% Steinway B2 sample C3

[y,Fs] = audioread("./SteinwayB2samples/Piano.mf.C3.aiff");
y_mono = sum(y,2)/size(y,2);
Nfft = 2^17;

y_fft = fft(y_mono,Nfft);
y_fft(1:ceil(20*Nfft/Fs))=0;

y_fft=abs(y_fft(1:Nfft/2));
f = Fs/2 * linspace(0,1,Nfft/2);
f1 = 131.2;
[p, idx] = findpeaks(y_fft, 'MinPeakHeight', 0.1, 'MinPeakDistance', 130 * Nfft / Fs);

figure();
plot(f,abs(y_fft));
hold on; 
p=p(1:25);
idx = idx(1:25);
stem(f(idx),p);
title("Spectrum and Partials Steinway B2: C3");
xlim([0 4000]);

%starting value for B coefficient
B = 10e-4;
delta = 1;
counter = 0;
old_trend = 1;

while(counter <100 && abs(delta) > 10e-4)
    fk_estimation = (1:25) .* f1 / (1+B)^(-0.5) .* (1+B.*(1:25).^2).^0.5;
    Dk = f(idx) - fk_estimation;
    trend = sum(sign(Dk(2:end) - Dk(1:end-1)));
    %delta = delta * sign(trend);
    if(trend >0)
            delta = delta * sign(delta);
    end
    if(trend<0)
        if(delta < 0)
            delta = delta * sign(delta);
        end
        if(delta > 0)
            delta = delta * sign(trend);
        end
    end
    if(sign(trend) ~= sign(old_trend))
        delta = delta/2;
    end
    old_trend = trend;
    B = B*10^delta;
    counter= counter+1;
end
% Theoretical values of harmonics, using f1 as fundamental
f_theoretical = (1:25) *f1/(1+B)^(-0.5) .* ((1+B.*(1:25).^2).^(0.5));
c = polyfit((1:25).^2, (f(idx(1:25)) ./(1:25)).^2, 1);
y_est = polyval(c, (1:25).^2);

figure();
plot((1:25).^2, y_est);
hold on;
plot((1:25).^2, (f(idx)./(1:25)).^2, 'or');
title("Inharmonicity of first 30 partials Steinway B2: C3");
ylabel("(f_{n}/n)^2");
xlabel("n^2");


%% Yamaha U3 sample C3

[y,Fs] = audioread("./YamahaU3samples/C3.wav");
y_mono = sum(y,2)/size(y,2);
Nfft = 2^17;

y_fft = fft(y_mono,Nfft);
y_fft(ceil(3140*Nfft/Fs):ceil(3160*Nfft/Fs)) = 0;
y_fft(1:ceil(20*Nfft/Fs))=0;

y_fft=abs(y_fft(1:Nfft/2));
f = Fs/2 * linspace(0,1,Nfft/2);
f1 = 130.5;
[p, idx] = findpeaks(y_fft, 'MinPeakHeight', 0.1, 'MinPeakDistance', 130 * Nfft / Fs);

figure();
plot(f,abs(y_fft));
hold on; 
p=p(1:25);
idx = idx(1:25);
stem(f(idx),p);
title("Spectrum and Partials Yamaha U3: C3");
xlim([0 4000]);

%starting value for B coefficient
B = 10e-4;
delta = 1;
counter = 0;
old_trend = 1;

while(counter <100 && abs(delta) > 10e-4)
    fk_estimation = (1:25) .* f1 / (1+B)^(-0.5) .* (1+B.*(1:25).^2).^0.5;
    Dk = f(idx) - fk_estimation;
    trend = sum(sign(Dk(2:end) - Dk(1:end-1)));
    %delta = delta * sign(trend);
    if(trend >0)
            delta = delta * sign(delta);
    end
    if(trend<0)
        if(delta < 0)
            delta = delta * sign(delta);
        end
        if(delta > 0)
            delta = delta * sign(trend);
        end
    end
    if(sign(trend) ~= sign(old_trend))
        delta = delta/2;
    end
    old_trend = trend;
    B = B*10^delta;
    counter= counter+1;
end
% Theoretical values of harmonics, using f1 as fundamental
f_theoretical = (1:25) *f1/(1+B)^(-0.5) .* ((1+B.*(1:25).^2).^(0.5));
c = polyfit((1:25).^2, (f(idx(1:25)) ./(1:25)).^2, 1);
y_est = polyval(c, (1:25).^2);

figure();
plot((1:25).^2, y_est);
hold on;
plot((1:25).^2, (f(idx)./(1:25)).^2, 'or');
title("Inharmonicity of first 30 partials Yamaha U3: C3");
ylabel("(f_{n}/n)^2");
xlabel("n^2");
