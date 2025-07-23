function [Dt_r, MSD] = fastMSDfft(X,Y)
% fastMSDfft computes the mean squared displacement using fft.
% At least two orders of magnitude faster than standard O(n^2) method.
%
% Sources:
%   https://medium.com/@pavel.kos/how-to-speed-up-msd-calculations-df3b670016f6
%   https://stackoverflow.com/questions/69738376/how-to-optimize-mean-square-
%   displacement-for-several-particles-in-two-dimensions/69767209#69767209
%
% Input:
%   X, Y - particle tracks where rows = particles, columns = times
%
% Output:
%   Dt_r, MSD - time indices, mean squared displacement
%
% Written by Joe Prisaznuk, used ChatGPT for assistance converting python
%   code to MATLAB syntax on 1/14/2024

tracks = cat(3,X,Y);        % put into one 3D array
nTime = size(tracks, 2);    % Number of time steps, same as size(X, 2)

% Compute the frequency domain using FFT
q = ifft(abs(fft(tracks, 2 * nTime, 2)).^2, [], 2);
S2 = sum(real(q(:, 1:nTime, :)), 3) ./ (nTime - (0:nTime-1));

D = sum(tracks.^2, 3); % Compute the direct squared distance
D = [D, zeros(size(tracks, 1), 1)]; % append zeros to the end

% Compute the first part of MSD
p1 = 2 * sum(D, 2) - cumsum([zeros(size(D, 1), 1), D(:, 1:end-1)] + flip(D, 2), 2);
S1 = p1(:, 1:end-1) ./ (nTime - (0:nTime-1));

MSD = S1 - 2 * S2;  % Final MSD result

Dt_r = 1:nTime;   % Time indices for MSD 
MSD = MSD(:, Dt_r); % Skip the first and last time step
Dt_r = Dt_r - 1;    % shift down such that delay = Dt_r * frame_interval