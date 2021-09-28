function [q2] = timestepper(tspan,q1,L)

disp(['--> Calling time stepper with Dt=' num2str(tspan(end))]);

LNSEfun = @(t,q1)LNSE(t,q1,L);

options = odeset('Stats','on');
[~,q2]  = ode45(LNSEfun,tspan,q1,options);

q2      = q2(end,:);

end
