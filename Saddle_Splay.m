function [outputArg1] = Saddle_Splay(lambda,tilt_0,diaphrag_outer_R,radius,chi)

n=(chi^2/4+1)^0.5;

A=0.5/(radius*lambda)*tilt_0^2*(radius/diaphrag_outer_R)^(-chi)*besseli(n,radius/lambda)/(besseli(n,diaphrag_outer_R/lambda))^2;
B = (besseli(n+1,radius/lambda)+besseli(n-1,radius/lambda)-chi*(lambda/radius)*besseli(n,radius/lambda));
outputArg1=A*B;
end