function [outputArg1] = Delta_splay(lambda,J_0m,tilt_0,diaphrag_outer_R,radius,chi)

n=(chi^2/4+1)^0.5;

A=-0.5/lambda*tilt_0/besseli(n,diaphrag_outer_R/lambda)*(radius/diaphrag_outer_R)^(-chi/2);
B = (besseli(n+1,radius/lambda)+besseli(n-1,radius/lambda)+(2-chi)*(lambda/radius)*besseli(n,radius/lambda));
outputArg1=A*B-J_0m;
end