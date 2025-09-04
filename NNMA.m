% Written by Shivaram Poigai Arunachalam
% Subroutine function to calculate NNMA- Nearest-Neighbor Moving
% Averaged time series based on the desired time scale factor 's'
% data - incoming time series data
% s - time scale factor chosen to estimate nearest neighbor moving averageNearest Neighbor moving average procedure

% Last Modified Februray 23 2016

function data = NNMA(data,s)
N = length(data);             % Estimate the length of the incoming time series data
for i = 1:N - s               % for loop to compute nearest neighbor moving average based on the time scale factor
    if i > s                  % check for the presence of backward neighbor
data(i-s) = mean(data((i-s):(i + s)));  % compute the avereage of the corresponding nearest neighbors
    end
    
end

% ---------------- END OF PROGRAM -------------------------------------