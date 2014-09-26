function [time_array, path] = BrownianPath(n, power)
% This function consumes a postive integer n, and power is either 1 or 2
% This function returns an array of Brownian paths with n sub-intervals

delt_t = 1/n;    % time increment
path_Z = zeros(1, n+1);  % array for storing Brownian Motion
path_Z(1) = 0;           % Z(t = 0) = 0

if (power == 1)
    for i = 2:n+1
        % from the formula for standard Brownian Motion
        path_Z(i) = path_Z(i-1) + sqrt(delt_t)*randn;
    end
else
    for i = 2:n+1
        % from the formula for squares of the Brownian Motion
        path_Z(i) = path_Z(i-1) + delt_t*(randn^2);
    end
    
end
time_array = 0:delt_t:1;
% return the array
path = path_Z;
