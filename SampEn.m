% Written by Shivaram Poigai Arunachalam
% Subroutine function to calculate NNMA- Nearest-Neighbor Moving
% Averaged time series based on the desired time scale factor 's'
% data - incoming time series data
% r - Predefined threshold for template matching chosen to be 20% of the
% standard deviation of the incoming time series data
% delay - time delay between sucessive time series points to be considered
% Reference: Wu, S. D., Wu, C. W., Lee, K. Y., & Lin, S. G. (2013). Modified multiscale entropy for short-term time series analysis. Physica A: Statistical Mechanics and its Applications, 392(23), 5865-5873.

% Last Modified Februray 23 2016

function MS_Entropy = SampEn(data,delay)
r = 0.2*std(data);
N = length(data);                        % Find the length of the time series
Nn = 0;                                  % Initial value for the number of matched template vectors for dimension 'm'
Nd = 0;                                  % Initial value for the number of matched template vectors for dimension 'm+1'
for i = 1:N - 3*delay                    % Loop to construct template vectors in 'm' dimension
for j = i + delay:1:N - 2*delay          % Loop to construct template vectors in 'm+1' dimension
if abs(data(i)-data(j))<r && abs(data(i + delay)-data(j + delay))<r  % Calculate the Euclidean distance of the tempalte vetors to be within the threshold 'r' in 'm' dimension
Nn = Nn + 1;                                                         % increase the # of matched template vectors by 1 if a matching template vector is found
if abs(data(i + 2*delay)-data(j + 2*delay))<r                        % Calculate the Euclidean distance of the tempalte vetors to be within the threshold 'r' in 'm+1' dimension
Nd = Nd + 1;                                                         % increase the # of matched template vectors by 1 if a matching template vector is found
end                  % end of if statement for 'Nd'
end                  % end of if statement for 'Nn'
end                  % end of for statement for 'j loop'
end                  % end of for statement for 'i loop'

MS_Entropy = -log(Nd/Nn);                    % Calcuate the Sample entropy as the ratio of# of matched template vectors in 'm+1' dimension to 'm' dimension

% ---------------- END OF PROGRAM -------------------------------------