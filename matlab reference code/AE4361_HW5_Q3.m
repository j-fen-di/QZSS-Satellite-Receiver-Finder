% jfendi's MATLAB code for AE 4361, Homework 5, Question 3 - Spring 2022
% Justin Effendi, Dr. Lightsey, AE 4361, February 24, 2022

format long

% measured pseudoranges for QZS 1, 2, 3, and 4 (in km)
rhos_act = [36536.1926 39061.5413 36817.8847 36735.0892];

% expected pseudoranges from problem 2b (in km)
rhos_exp = [41037.420 45257.860 42161.310 40959.990];

% satellite conditions for satellites 1, 2, 3, 4 respectively of QZSS
r_QZSS = [-24371.056 32026.357 -8027.077;
    -24905.318 22853.526 30095.105;
    -25385.052 33662.647 36.793;
    -34500.873 20615.646 -7899.728];

%initial guess <x_u, y_u, z_u, t_u>
guess_0 = zeros([1 4]);

% linearize to user location
receiver_pos = trilaterate_to_loc(guess_0, r_QZSS, rhos_exp, rhos_act);
% print results
disp("3)");
fprintf('The location of the receiver in ECEF coordinates are <');
fprintf('%.4f, ', receiver_pos(1:end-2));
fprintf('%.4f> km.', receiver_pos(end-1));
disp(' ');
fprintf('The receiver clock offset is ');
fprintf('%.10f seconds.', receiver_pos(end));
disp(' ');
disp(' ');

% BONUS PROBLEM
% convert to lat, lon, h
[lat, lon, h] = findLatLonH(receiver_pos(1), receiver_pos(2), receiver_pos(3));
receiver_loc = [lat lon h];
% print results, gives coordinates of Takesaki Observation Stand in Japan
disp("Bonus problem:");
disp("The latitude of the receiver is " + num2str(receiver_loc(1), 8) + " degrees.");
disp("The longitude of the receiver is " + num2str(receiver_loc(2), 9) + " degrees.");
disp("The altitude of the receiver is " + num2str(receiver_loc(3), 9) + " km.");

function [r_user_f] = trilaterate_to_loc(r_user, r_sats, pseudo_exp, pseudo_act)
    % initial values needed
    c = 299792.000; % speed of light, in km/s
    
    % measured - expected pseudoranges
    delta_rho = pseudo_act - pseudo_exp;
    
    % initialize matrix H
    H = [0 0 0 c;
        0 0 0 c;
        0 0 0 c;
        0 0 0 c];
    
    % edit H to proper form
    for i = 1:4
        r_i_hat = sqrt((r_sats(i,1) - r_user(1))^2 + (r_sats(i,2) - r_user(2))^2 + (r_sats(i,3) - r_user(3))^2);
        for j = 1:3
            a = -(r_sats(i, j) - r_user(j)) / r_i_hat;
            H(i, j) = a;
        end
    end
    
    % calculate correction term (inverse H times transposed delta_rho)
    delta_x_u = H\delta_rho.';
    delta_x_u = delta_x_u.'; % transpose delta_x_u
    
    % update r_user with correction term
    r_user = r_user + delta_x_u;
    
    % come up with pseudo_guess
    pseudo_guess = zeros([1 4]); % initialize pseudo_guess
    for i = 1:4
        pseudo_guess(i) = sqrt((r_sats(i,1) - r_user(1))^2 + (r_sats(i,2) - r_user(2))^2 + (r_sats(i,3) - r_user(3))^2) + c*r_user(4);
    end
    
    % set final answer if converges, recurse otherwise
    if isConvergence(pseudo_guess, pseudo_act) == 1
        r_user_f = r_user;
    else    
        r_user_f = trilaterate_to_loc(r_user, r_sats, pseudo_guess, pseudo_act); % recurse until condition met
    end
end

function [x] = isConvergence(a, b)
    % function checking if convergence criteria is met
    x = 0; % pseudo boolean variable
    switchy = 0; % true if not converge, false if converges
    
    % compare each element of the two vectors
    for i = 1:length(a)
        if (a(i) - b(i)) > 1*10^-4
            switchy = 1;
        end
    end
    
    % return true if convergence
    if switchy == 0
        x = 1;
    end
end

function [lat, lon, h] = findLatLonH(x, y, z)
    % converts ECEF vector to lat, lon, h on Earth surface
    % needed constants
    r_E = 6371; % radius of Earth, in km
    
    % solving for lon, lat, and h using algebra
    lon = atan2d(y, x); % calculating latitude (in degrees)
    lat = atand(z*sind(lon) / y); % calculating longitude (in degrees)
    h = z/sind(lat) - r_E; % calculating altitude (in meters)
end