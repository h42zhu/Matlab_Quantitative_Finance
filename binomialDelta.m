function [V, L] = binomialDelta(S0, r, sigma, T, N, fpayoff)
% This function consumes six parameters
% S0 is the current stock price
% r is the risk-free interest rate
% sigma is the volatility
% T is the expiry time
% N is the number of timesteps
% fpayoff specifies a function which returns the payoff of an option

% It returns the initial option value V, and a structured array L, 
% with two fields, S and alpha, representing the underlying price 
% and the delta hedging position

delt = T/N;  % timestep size

% binomial lattice parameters
u = exp(sigma*sqrt(delt));
d = 1 ./ u;
a = exp(r*delt);
p = (a-d)/(u-d);    

% range of possible asset values at the end of the tree
W = S0*d.^([N:-1:0]').*u.^([0:N]'); 
St = W;   

% compute the payoff of the option at the final nodes
W = fpayoff(W);

% create the structure array
field1 = 'S';  value1 = {};
field2 = 'alpha';  value2 = {};
L = struct(field1,value1,field2,value2);

% Note: L(N) is values at t = t(N-1)
% backward recursion for option value
for i = N:-1:1
        
    % compute the delta hedging position
    alphas = (W(2:i+1) - W(1:i))./(St(2:i+1)-St(1:i));
    % binomial backward recursion
    W = exp(-r*delt)*(p*W(2:i+1) + (1-p)*W(1:i));
    % asset value for the previous step
    St = u*St(1:i);  
    % store the vector of values in the structure array
    L(i).S = St;
    L(i).alpha = alphas;

end    
% end for loop
% return option value
V = W(1);


end
