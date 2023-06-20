% jfendi's MATLAB code for AE 4361, Homework 5, Question 2 - Spring 2022
% Justin Effendi, Dr. Lightsey, AE 4361, February 24, 2022

format long

% PROBLEM 2A
% initial values for latitude (phi), longitude (lambda), and altitude (h)
% for the four satellites
% NOTE: each subvector represents a satellite (i.e. QZS-1, QZS-2, etc.)
QZSS = [-11.28 127.27 34666.42;
    41.68 137.46 38886.86;
    0.05 127.02 35790.31;
    -11.12 149.14 34588.99];

% initialize and allocate space for array of QZSS satellite positions
QZSS_ECEF = zeros([4 3]);

% print '2a)'
disp('2a)');

% find ECEF position for each satellite and print results (2a)
for i = 1:4
    % calculate position in ECEF frame with to_r_ECEF
    r_ECEF = to_r_ECEF(QZSS(i, 1), QZSS(i, 2), QZSS(i, 3));
    % add to QZSS_ECEF array for future usage
    QZSS_ECEF(i, 1:3) = r_ECEF(1:3);
    % print results
    fprintf('r of the QZS-');
    fprintf('%g', i);
    fprintf(' satellite is <');
    fprintf('%.3f, ', r_ECEF(1:end-1));
    fprintf('%.3f> km.\n', r_ECEF(end));
end

% PROBLEM 2B
% print '2b)'
disp('2b)');
r_origin_ECEF = [0 0 0];
c = 299792.000; % in km/s
t_u = 0; % in seconds (1 day)
for i = 1:4
    rho_hat = sqrt((QZSS_ECEF(i, 1))^2 + (QZSS_ECEF(i, 2))^2 + (QZSS_ECEF(i, 3))^2) + c*t_u;
    % print results
    fprintf('pseudorange of the QZS-');
    fprintf('%g', i);
    fprintf(' satellite is ');
    fprintf('%.3f', rho_hat);
    fprintf(' km.');
    disp(' ');
end

% converts to lat, lon, h to ECEF reference frame
function [r_ECEF] = to_r_ECEF(lat, lon, h)
    r_E = 6371; % radius of Earth, in km
    x = (r_E + h) * (cosd(lon)*cosd(lat));
    y = (r_E + h) * (sind(lon)*cosd(lat));
    z = (r_E + h) * sind(lat);
    r_ECEF = [x y z];
end