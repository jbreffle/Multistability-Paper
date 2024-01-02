% phi_logistic.m 
% simply supplies the logistic function with input x, half-max mu (=x_th),
% and steepness parameter alpha (= Delta)
function y = phi_logistic(x,mu,alpha)
y = 1./(1 + exp ( (mu-x)/alpha) );
end

