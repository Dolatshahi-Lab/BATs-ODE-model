% The ODEs for a simple predator/prey model with 2 equations and 4 parameters.
% This model is from page 6, section 3.1 of "A Methodology For Performing
% Global Uncertainty And Sensitivity Analysis In Systems Biology", Marino, et.
% al., Journal of Theoretical Biology, 2008-09-07, doi:
% 10.1016/j.jtbi.2008.04.011.

% From the above cited paper, with y(1) = Q and y(2) = 

% It has two state variables (Q, P) and several parameters. Q represents the
% density of prey, P represents the density of predators, alpha is the intrinsic
% rate of prey population increase, beta is the predation rate coefficient, sigma is
% the predator mortality rate, and delta is the reproduction rate of predators per
% prey consumed.

% This file can be used as a template for creating a new ODE model. To do so
% copy this file and edit it, replacing the predator/prey parameters and
% equations with those for the new model.
%
% When creating a new model, you will also need to copy and edit file
% lhs_ode_predator_prey_settings_new.m to create a settings file for the new
% model, to define the model parameters values, initial conditions, etc.

function dy = RB_lhs_ode_predator_prey_ode_c(t, y, params,MN2,ML2,MP2,MD2,CN2,CL2,CP2,CD2,TN2,TL2,TP2,TD2)

neq = size(y,1);
dy = zeros(neq,1);


Matp = 0.1287;%doubling time
Mmax = 1.4;%max tumor cells
Catp = 0.06;%doubling time
Cmax = 3.2;%max tumor cells
Tatp = 0.16;%doubling time
Tmax = 1.2;%max tumor cells

aKill_N = (params(1));
aKill_P = (params(2));
aKill_L = (params(3));
aKill_D = (params(4));


%%MCF7
%dTumor/dt y(1)
dy(1) = Matp*y(1)-(Matp*y(1)^2)/Mmax...
       - aKill_N*y(2)*y(1)...
       - aKill_L*y(3)*y(1)...
       - aKill_P*y(4)*y(1)...
       - aKill_D*y(5)*y(1);

%dBAT_N/dt y(5)
dy(2) = ppval(MN2,t);

%dBAT_L/dt y(8)
dy(3) = ppval(ML2,t);

%dBAT_P/dt y(11)
dy(4) = ppval(MP2,t);

%dBAT_D/dt y(11)
dy(5) = ppval(MD2,t);


%%CAMA1
%dTumor/dt y(1)
dy(6) = Catp*y(6)-(Catp*y(6)^2)/Cmax...
       - aKill_N*y(7)*y(6)...
       - aKill_L*y(8)*y(6)...
       - aKill_P*y(9)*y(6)...
       - aKill_D*y(10)*y(6);

%dBAT_N/dt y(5)
dy(7) = ppval(CN2,t);

%dBAT_L/dt y(8)
dy(8) = ppval(CL2,t);

%dBAT_P/dt y(11)
dy(9) = ppval(CP2,t);

%dBAT_D/dt y(11)
dy(10) = ppval(CD2,t);


%%MCF7
%dTumor/dt y(1)
dy(11) = Tatp*y(11)-(Tatp*y(11)^2)/Tmax...
       - aKill_N*y(12)*y(11)...
       - aKill_L*y(13)*y(11)...
       - aKill_P*y(14)*y(11)...
       - aKill_D*y(15)*y(11);

%dBAT_N/dt y(5)
dy(12) = ppval(TN2,t);

%dBAT_L/dt y(8)
dy(13) = ppval(TL2,t);

%dBAT_P/dt y(11)
dy(14) = ppval(TP2,t);

%dBAT_D/dt y(11)
dy(15) = ppval(TD2,t);

end
