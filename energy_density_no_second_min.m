function [outputArg1] = energy_density_no_second_min(lambda,J_0m,tilt_0,diaphrag_outer_R,radius,chi,kappa_m)

DJ= Delta_splay(lambda,J_0m,tilt_0,diaphrag_outer_R,radius,0);
K= Saddle_Splay(lambda,tilt_0,diaphrag_outer_R,radius,0);
T=Tilt(lambda,tilt_0,diaphrag_outer_R,radius,0);
outputArg1=0.5*kappa_m*DJ^2+chi*kappa_m*K+0.5*kappa_m/lambda^2*T^2-0.5*kappa_m*J_0m^2;

end