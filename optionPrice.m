function [price] = optionPrice(S0,K,T,r,sigma,opttype, Eu_Am, Power, Nsteps)
% function for price options with binomial tree, consumes three integers
% and produces the price of the option, which is a float
% S0 = current stock price
% K = strike price
% T = expiry time
% r = risk-free interest rate
% sigma = volatility
% opptype is an integer of either 1 or 0; 0 for call, 1 for put
% Eu_Am is an integer of either 1 or 0; 0 for European option, 1 for
% American option
% Power is for Power Options
% Nsteps is an integer greater than 0 denoting the number of time steps

delt = T/Nsteps;  % timestep

% tree parameters

u = exp(sigma*sqrt(delt));
d = 1 ./ u;
a = exp(r*delt);
p = (a-d)/(u-d);    


W = S0*d.^([Nsteps:-1:0]').*u.^([0:Nsteps]'); % range of possible asset values at the end of the tree
S = W;  % storing the values of the asset price for American Options
    
% payoff at t = T

    if (opttype == 0)  
        W = (max(W-K,0)).^Power;   % call option
    else 
        W = (max(K-W,0).^Power);   % put option
    end   
    %end if

if (Eu_Am == 0) 
% European option case

	% backward recursion
	for i = Nsteps:-1:1
        W = exp(-r*delt)*(p*W(2:i+1) + (1-p)*W(1:i));
	end     
    
else
% American option case

    % backward recursion
    for i = Nsteps:-1:1
       S = u*S(1:i);   % asset value for the previous step
       W = exp(-r*delt)*(p*W(2:i+1) + (1-p)*W(1:i));  % European option value
       if (opttype == 0) 
           payoff = S-K;   % payoff at time t for American Call if exercise
       else
           payoff = K-S;   % payoff at time t for American Put if exercise
       end
       W = max(W, payoff); 
       % choose the higher value of the payoff or the value of holding the option
    end

   
end
% end if

price = W(1);  % price of the option, output of the function

end 
% end function
