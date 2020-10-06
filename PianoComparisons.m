clear;
close;
clc;

% This script makes a comparison between the Inharmonicity coefficients for
% 5 different notes of an Upright and a Grand piano. It makes use of the
% code implemented in the "InharmonicityAnalysis" script and plots the
% comparison results.

%% Array construction
% DATA TO CYCLE OVER NOTES
B2fundamentals = [32.323 65.2 131.1 262.1 524.9];
U3fundamentals = [31.73 64.94 130.5 261.4 523.9];
IdealFundamentals = [32.703 65.406 130.813 261.626 523.251];

B2names = ["Piano.mf.C1.aiff", "Piano.mf.C2.aiff", ...
    "Piano.mf.C3.aiff", "Piano.mf.C4.aiff", "Piano.mf.C5.aiff"];
U3names = ["C1.wav", "C2.wav", "C3.wav", "C4.wav", "C5.wav"];

BcoeffsB2 = zeros(1,5);
BcoeffsU3 = zeros(1,5);

BinterpB2 = zeros(1,5);
BinterpU3 = zeros(1,5);

%% B2 Loop

for j=1:5
    [y,Fs] = audioread(strcat("./SteinwayB2samples/", B2names(j)));
    %[y,Fs] = audioread("./YamahaU3samples/C5.wav");
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
    f1 = B2fundamentals(j);
    deltaF = 0.4 * f1;
    
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

    BcoeffsB2(j) = B;
    BinterpB2(j) = c(1) / f1^2;
%     
%     % Plots
%     figure(1);
%     subplot(2,4,j);
%     plot(f,abs(y_fft));
%     title("Spectrum Steinway B2: C"+j);
%     xlabel("f [Hz]");
%     ylabel("FFT");
%     xlim([0, ceil(f_peaks(length(f_peaks))/Nfft * Fs /10)*10]);
%     hold on;
%     plot(f(f_peaks), peaks, 'or');
%     plot(f_theoretical, peaks, '.b');
%     legend('FFT','Spectrum Peaks', 'Computed Peaks');
% 
%     figure(2);
%     subplot(2,4,j);
%     plot(f, y_fft);
%     hold on;
%     title("Ideal vs Real Partials Steinway B2: C"+j);
%     xlabel("f [Hz]");
%     ylabel("FFT");
%     xlim([0, ceil(f_peaks(15)/Nfft * Fs /10)*10]);
%     plot(f(f_peaks), peaks, 'or');
%     stem(f_ideal, peaks);
%     legend('FFT','Real partials','Ideal partials');
% 
% 
%     figure(3);
%     subplot(2,4,j);
%     plot((1:25).^2, y_est);
%     hold on;
%     plot((1:25).^2, (f(f_peaks)./(1:25)).^2, 'or');
%     title("Partials Steinway B2: C"+j);
%     ylabel("(f_{n}/n)^2");
%     xlabel("n^2");
    
    
end

%% U3 Loop

for j=1:5
    [y,Fs] = audioread(strcat("./YamahaU3samples/", U3names(j)));
    %[y,Fs] = audioread("./YamahaU3samples/C5.wav");
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
    f1 = U3fundamentals(j);
    deltaF = 0.4 * f1;
    
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

    BcoeffsU3(j) = B;
    BinterpU3(j) = c(1) / f1^2;
%     
%     % Plots
%     figure(1);
%     subplot(2,4,4+j);
%     plot(f,abs(y_fft));
%     title("Spectrum Yamaha U3: C"+j);
%     xlabel("f [Hz]");
%     ylabel("FFT");
%     xlim([0, ceil(f_peaks(length(f_peaks))/Nfft * Fs /10)*10]);
%     hold on;
%     plot(f(f_peaks), peaks, 'or');
%     plot(f_theoretical, peaks, '.b');
%     hold off;
%     legend('FFT','Spectrum Peaks', 'Computed Peaks');
% 
%     figure(2)
%     subplot(2,4,4+j);
%     plot(f, y_fft);
%     hold on;
%     title("Ideal vs Real Partials Yamaha U3: C" + j);
%     xlabel("f [Hz]");
%     ylabel("FFT");
%     xlim([0, ceil(f_peaks(15)/Nfft * Fs /10)*10]);
%     plot(f(f_peaks), peaks, 'or');
%     stem(f_ideal, peaks);
%     hold off;
%     legend('FFT','Real partials','Ideal partials');
% 
% 
%     figure(3);
%     subplot(2,4,4+j);
%     plot((1:25).^2, y_est);
%     hold on;
%     plot((1:25).^2, (f(f_peaks)./(1:25)).^2, 'or');
%     title("Partials Yamaha U3: C"+j);
%     ylabel("(f_{n}/n)^2");
%     xlabel("n^2");
%     hold off;
    
end

%% Plot B comparison
figure();
plot(1:5, BcoeffsB2, 'or');
hold on;
plot(1:5, BcoeffsU3, 'ob');
legend('B2 inharmonicity', 'U3 inharmonicity');
title("Inharmonicity coefficient comparison Yamaha U3 - Steinway B2");
ylabel("B");
xlabel("C_{1-4}");

%% Plot Fundamentals comparison
notes = [4 16 28 40 52];
B2shift = B2fundamentals - IdealFundamentals;
U3shift = U3fundamentals - IdealFundamentals;

figure();
plot(notes, zeros(1,length(notes)), 'DisplayName', 'Ideal tuning');
hold on
plot(notes, B2shift, 'DisplayName', 'B2 stretched tuning');
plot(notes,U3shift, 'DisplayName', 'U3 stretched tuning');
legend('show');
xlim([4,52]);
xlabel("Note number");
ylabel("f_{fund} - f_{ideal}");