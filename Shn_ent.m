
function SE = Shn_ent(g)
L = length(g);
mx = max(g);            % Maximum intensity value of the optical  intensity
mi = min(g);            % Minimum intensity value of the optical  intensity
bin = linspace(mi,mx,1000);      % Create amplitude bins with 1 unit intensity difference for bin size
h = hist(g,bin);        % Create histogram data using specified bin size
l = length(h);          % Calculate the size of histogram matirx
P = zeros(1,l);         % Create array for probability density function
S = zeros(1,l);         % Create array to be used for Shannon Entropy calculations
N = zeros(1,l);         % Create array to be used for Shannon Entropy calculations
for i = 1:l
    P(i) = h(i)/L;              % Probability that a given sample falls within a particular amplitude bin
    S(i) = P(i)*log2(P(i));     % Shannon Entropy Term
    N(i) = S(i);
    if isnan(N(i))             % Check for NAN and replace it with 0
        N(i) = 0;
    end
end
SE = -1*sum(N);                 % Calculate Shannon Entropy
 
