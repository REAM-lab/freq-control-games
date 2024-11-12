function output = integration_xQx(tin, tfin, Q, tsr)
% Function to compute the integral of x'Qx. It can be used for u'Ru in LQR
% or linear quadratic differential games.
% tsr : timeseries in which all the states are saved in columns.
% tin : initial time. It also the lower bound of the integral
% tfi : final time. It is also upper bound of the integral
% Q :  weighted matrix  
% x :  state vector in the form:
%        x1 x2 x3 ....
%     t1 
%     t2 
%
% Note 1: The trapz function requies that the size (time axis) of state 
% vector must coincide with integral domain. 
%
% By P.Serna-Torre (2024)

timepoints = tsr.Time ; % vector of timepoints

t_mask = (tin <= timepoints) & (timepoints < tfin); % logical matrix 

int_domain = timepoints(t_mask); % integration domain 

x = tsr.Data(t_mask,:) ; % x within the integration domain

xss = x(end,:); % end of the vector x above

xdev = x - xss; % deviation state. 

nintervals = length(int_domain);

output.integrand = zeros(nintervals,1);

for t=1:nintervals 

    x_at_t = xdev(t,:)';    % here we integrate the deviation x, but if
                            % you want to integrate x, then put x instead
                            % of xdev.

    output.integrand(t,1) = x_at_t'*Q*x_at_t;

end

output.integral = trapz(int_domain, output.integrand);

output.tps = int_domain;

return