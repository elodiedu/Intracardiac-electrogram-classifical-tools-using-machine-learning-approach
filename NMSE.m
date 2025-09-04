% Written by Shivaram Poigai Arunachalam
% Subroutine function to calculate Nearest-Neighbor Multi-Scale Entropy
% (NMSE) for the time series data
% data - incoming time series data
% scale - time scale factor chosen to estimate nearest neighbor moving average
% r - Predefined threshold for template matching chosen to be 20% of the
% standard deviation of the incoming time series data
% MS_E - Returns the multi-scale entropy value of the time series
% This subroutine uses two functions namely NNMA- Nearest-Neighbor Moving
% Average and SampEn which computes the sample entropy of the new time
% series data

% Last Modified: Feb 23 2016

function MS_E = NMSE(data, scale)
r = 0.2*std(data);               % Threshold to define matching template vectors
NT = NNMA(data, scale);              % Compute the new time (NT) series Nearest-Neighbor moving average based on the defined time scale factor
MS_E = SampEn(NT,r,1);        % Compute the Sample Entropy of the new time series


% ---------- END OF PROGRAM -----------------------------------