function DF = dominant_freq(g,Fs)
L = length(g);                      % Length of the input signal
Y = abs(fft(g))/L;                       % Compute the Fast Fourier Transform
t = (0:L-1)*(1/Fs);                 % Define time vector for plotting based on Sampling Rate
f = Fs*linspace(0,1,length(Y));% Define the corresponding frequency vector for plotting
Y(1)=0;                            % make the dc value to zero to find dominant frequency peak
[~,v] = max(Y);                    % find max frequency's location 
DF = f(v);                         % find max frequency