function [price, delta, gamma] = AMER_PDESolve(S0,S,K,T,sig,r,Nsteps,optionType, timestepping, dnorm)
% This function consumes nine parameters
% S0 is the initial asset price
% S is the asset price grid vector
% r is risk-free rate
% sigma is the volatility
% K is strike price of the option
% T is the time to expiry
% N is the number of timesteps
% optionType is a string that's either 'call' or 'put'
% timestepping is a string that's either 'constant' or 'variable'
% dnorm is used to timestep selector
% This function returns price, which represents the Black-Scholes 
% price of the American option at S0, using Penalty Method, and returns 
% The delta and gamma of the option at t = 0
% The timestepping method used is CN-Rannacher

dtau = T/Nsteps; % timestep size in tau direction
row = length(S); % number of node points for asset prices

tol = sqrt(eps);
Large = 1/tol;
% pre-load array/matrix for option values
% row represents the tau direction
% columns presents the asset price direction

V = zeros(row, Nsteps+1); 

% determine the payoff for the option at tau = 0
if strcmpi(optionType, 'call')
    V(:, 1) = max(S'-K,0);  
    V_star = V(:, 1);
else
    V(:, 1) = max(K-S',0);
    V_star = V(:, 1);
end

% determine the coefficients for central weighting
Sf = S(3:row);
Sb = S(1:row-2);
Sm = S(2:row-1);
sig2 = sig*sig;

% central difference
aj = sig2*(Sm.^2)./((Sm-Sb).*(Sf-Sb)) - (r*Sm./(Sf-Sb));
bj = sig2*(Sm.^2)./((Sf-Sm).*(Sf-Sb)) + (r*Sm./(Sf-Sb));

% forward difference
af = sig2*(Sm.^2)./((Sm-Sb).*(Sf-Sb));
bf = sig2*(Sm.^2)./((Sf-Sm).*(Sf-Sb)) + (r*Sm./(Sf-Sm));

% upstream
aj = aj.*(aj > 0) + af.*(aj < 0);
bj = bj.*(bj > 0) + bf.*(bj < 0);

arow = length(aj);
I = eye(row);

[M, M_hat] = triMatrix(dtau, aj, bj, r, row, arow);

if strcmpi(timestepping, 'constant')
    % constant timestepping
    last_col = Nsteps+1;
    for i = 2:last_col 
        V(:, i) = V(:, i-1);
        cond = true;
        while cond
            P = Large*(V(:, i) < V_star);
            P = spdiags(P, 0, row, row);
            if i <=3 
                % fully implicit
                V_old = V(:, i);
                Y = V(:, i-1)+P*V_star;
                [L_fi, U_fi] = lu(I + M + P);
                V_new = U_fi\(L_fi\Y);
                V(:, i) = V_new;
            else
                % Crank-Nicolson
                V_old = V(:, i);
                Y = (I-M_hat)*V(:, i-1)+P*V_star;
                [L_cn, U_cn] = lu(I + M_hat + P);                
                V_new = U_cn\(L_cn\Y);
                V(:, i) = V_new;
            end
            % determine if 
            maxi = max(abs(V_new - V_old)./(max(1, abs(V_new))));
            cond = (maxi >= tol);
            
        end % end while loop
    end % end for loop
        
else
    % variable timestepping
    D = 1;
    Tau = T;
    i = 2;
    while Tau > 0
        V(:, i) = V(:, i-1);
        cond = true;
        while cond
            P = Large*(V(:, i) < V_star);
            P = spdiags(P, 0, row, row);
            if i <=3 
                % fully implicit
                V_old = V(:, i);
                Y = V(:, i-1)+P*V_star;
                [L_fi, U_fi] = lu(I + M + P);
                V_new = U_fi\(L_fi\Y);
                V(:, i) = V_new;
            else
                % Crank-Nicolson
                V_old = V(:, i);
                Y = (I-M_hat)*V(:, i-1)+P*V_star;
                [L_cn, U_cn] = lu(I + M_hat + P);                
                V_new = U_cn\(L_cn\Y);
                V(:, i) = V_new;
            end
            % determine if 
            maxi = max(abs(V_new - V_old)./(max(1, abs(V_new))));
            cond = (maxi >= tol);
            
        end % end while loop
        
        Tau = Tau - dtau;
        % update the new dtau
        Vm = max(abs(V(:, i)), abs(V(:, i-1)));
        Vm = max(D, Vm);
        MaxRelChange = max(abs(V(:, i)-V(:, i-1))./ Vm);
        dtau = dnorm/MaxRelChange*dtau;
        
        % disp(dtau)
        
        if Tau - dtau < 0
            dtau = Tau;
        end 
        % update the tridiagonal matrix
        if i <= 3
            [M, ~] = triMatrix(dtau, aj, bj, r, row, arow);
        else
            [~, M_hat] = triMatrix(dtau, aj, bj, r, row, arow);
        end
        
        i = i + 1;
    end % end Tau while loop
    
    last_col = i - 1;
end % end if


dX = (Sf - Sb)';
fV = V(:, last_col);
dV = fV(3:end) - fV(1:end-2);
delta = dV./dX;

ddV = (fV(3:end) - fV(2:end-1))./(Sf - Sm)' ...
    + (fV(1:end-2) - fV(2:end-1))./(Sm - Sb)';
ddX = dX/2;
gamma = ddV./ddX;

% price = V;
% return the price of the American option from linear interpolation
price = interp1(S', fV, S0);
delta = interp1(Sm', delta, S0);
gamma = interp1(Sm', gamma, S0);

end
