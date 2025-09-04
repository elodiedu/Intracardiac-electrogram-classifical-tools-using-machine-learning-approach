% Elizabeth Annoni
% Last Modiefied: 02 - 24 - 16

% Function to calculate the Frequency estimation for a 1D signal (in time)
% Algorithm based off of:
%       Local Multiscale Frequency and Bandwidth Estimation
%       Hans Knutsson, Carl-Fredrik Westin, Gosta Granlund
% INPUTS:
%       signal_1D - input time series signal
%       T - period, 1/sampling frequency
% OUTPUTS:
%       FE_w - Frequency estimation in radians per second (Array, same size
%       as signal_1D) 
%       FE_Hz - Frequency estimation in Hertz (Array, same size as 
%       signal_1D) 
%       DF - Dominant Frequency for time series (single value)

% Estimate the instantaneous frequency of a 1D signal using a weighted
% summation of 8 log-Gabor filters. (Equation 29 in the paper reference
% above). 

function [FE_Hz,FE_w] = MSF_1D(signal_1D, Fs, DF)


T = 1/Fs;
L = length(signal_1D);
t = (0:L-1)*T;


y = signal_1D;

NFFT = 2^nextpow2(L);
Y = fft(y,NFFT);
f = Fs*linspace(0,1,NFFT);
fin = 2*pi*f;

%% Remove Harmonics
y = double(y);
if (DF ~= 0)
    maxFreq = 50;
    fsample = 1/T;
    q = 10;
    % Remove 1st harmonic
    df1 = DF*2;
    if (df1 < maxFreq)
        %notch filter centered at df1
        w0 = df1/(fsample/2);
        bw = w0/q;
        [b1,a1] = iirnotch(w0,bw,q);
        y1 = filter(b1,a1,y);
        y = y1;
    end
    % Remove 2nd harmonic
    df2 = DF*3;
    if (df2 < maxFreq)
        %notch filter centered at df2
        w0 = df2/(fsample/2);
        bw = w0/q;
        [b2,a2] = iirnotch(w0,bw,q);
        y2 = filter(b2,a2,y);
        y = y2;
    end
    % Remove 3rd harmonic
    df3 = DF*4;
    if (df3 < maxFreq)
        %notch filter centered at df3
        w0 = df3/(fsample/2);
        bw = w0/q;
        [b3,a3] = iirnotch(w0,bw,q);
        y3 = filter(b3,a3,y);
        y = y3;
    end
    
end

Y = fft(y,NFFT);

%% Apply log-Gabor filters to input signal

p0 = pi/3;
N = 8;
pis = [];
for i = 1:N
    pi_i = 2^i * p0;
    pis = [pis, pi_i];
end
% pis = [pi/128, pi/64, pi/32, pi/16, pi/8, pi/4, pi/2, pi];

filtSigTime = [];
for i  = 1:length(pis)
    B = 2*sqrt(2);
    Cb = 4/(B^2 * log(2));
    R2_i = exp(-Cb*log(fin/pis(i)).^2);
    filtSig_i = (Y).*R2_i;
    filtSigTime_i = ifft(filtSig_i);
    filtSigTime = [filtSigTime; filtSigTime_i];
end

 

%% Frequency Estimation
removeEdge = 50;
for j = removeEdge+1:length(signal_1D)-removeEdge
    numTemp = 0;
    for m = 1:N-1
        numTemp = numTemp + ((2^(m+0.5))*(filtSigTime(m+1,j)));
    end
    denTemp = 0;
    for m = 1:N-1
        denTemp = denTemp + (filtSigTime(m,j));
    end
    
    rho(j) = p0 * (1/denTemp) * numTemp;
end

% See how smoothing affects output
% figure
% plot(abs(rho));
% hold on
% plot(smooth(abs(rho),10));
% plot(smooth(abs(rho),50));
% plot(smooth(abs(rho),100));
% legend('Original','Smooth:10','Smooth:50', 'Smooth: 100');


% figure
% plot(abs(rho)/(2*pi));
% 
% avgFE_w = mean(abs(rho));
% avgFE_Hz = mean(abs(rho))/(2*pi);

% Output the array of frequency estimates
FE_w = mean(abs(rho));
FE_Hz = mean(abs(rho)/(2*pi));
end




