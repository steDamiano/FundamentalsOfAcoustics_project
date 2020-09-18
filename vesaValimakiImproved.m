clear;
close;
clc;

%% Steinway B2 sample C1

[y,Fs] = audioread("./SteinwayB2samples/Piano.mf.C3.aiff");
y_mono = sum(y,2)/size(y,2);
Nfft = 2^17;

y_fft = fft(y_mono,Nfft);
y_fft(1:ceil(20*Nfft/Fs))=0;

y_fft=abs(y_fft(1:Nfft/2));
f = Fs/2 * linspace(0,1,Nfft/2);

f1 = 131.2;
deltaF = 0.4 * f1;

figure();
plot(f,abs(y_fft));

%[p, idx] = findpeaks(y_fft, 'MinPeakHeight', 1.7, 'MinPeakDistance', 31.5 * Nfft / Fs);
% hold on; 
% p=p(1:25);
% idx = idx(1:25);
% stem(f(idx),p);
% 
% %peaks evaluation
% deltaF = 0.4 * f1;
% B = 0.00025;
% fk_estimation = (1:25) * f1 .* sqrt(1 + B * (1:25).^2);
% 
% peaks = zeros(1,25);
% f_peaks = zeros(1,25);
% 
% for i=1:length(fk_estimation)
%     %set interval 
%     lowLimit = ceil(fk_estimation(i) * Nfft/Fs - deltaF);
%     upLimit = ceil(fk_estimation(i)*Nfft/Fs + deltaF);
%     
%     peaks(i) = max(y_fft(lowLimit:upLimit));
%     f_peaks(i) = find(y_fft == peaks(i));
%     %stem(f(ceil(i*f1 * Nfft/ Fs)), 200, 'or');
%     %stem(fk_estimation(i), 200, 'ob');
% 
% end
% figure();
% plot(f,abs(y_fft));
% hold on;
% stem(f(f_peaks), peaks, 'or');
% 

%starting value for B coefficient
B = 0.0001;
delta = 1;
counter = 0;
old_trend = 1;
% Theoretical values of harmonics, using f1 as fundamental

while(counter <100 && abs(delta) > 10e-4)
    fk_estimation = (1:25) .* f1 / (1+B)^(-0.5) .* (1+B.*(1:25).^2).^0.5;
    
    for i=1:length(fk_estimation)
    %set interval 
    lowLimit = ceil(fk_estimation(i) * Nfft/Fs - deltaF);
    upLimit = ceil(fk_estimation(i)*Nfft/Fs + deltaF);
    
    peaks(i) = max(y_fft(lowLimit:upLimit));
    f_peaks(i) = find(y_fft == peaks(i));
    %stem(f(ceil(i*f1 * Nfft/ Fs)), 200, 'or');
    %stem(fk_estimation(i), 200, 'ob');

    end
    
    Dk = f(f_peaks) - fk_estimation;
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

figure();

c = polyfit((1:25).^2, (f(f_peaks) ./(1:25)).^2, 1);
y_est = polyval(c, (1:25).^2);
plot((1:25).^2, y_est);
hold on;
plot((1:25).^2, (f(f_peaks)./(1:25)).^2, 'or');
title("Inharmonicity of first 30 partials Steinway B2");
ylabel("(f_{n}/n)^2");
xlabel("n^2");
