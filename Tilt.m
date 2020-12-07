function [outputArg1] = Tilt(lambda,tilt_0,diaphrag_outer_R,radius,chi)
n=(chi^2/4+1)^0.5;
outputArg1=tilt_0*(radius/diaphrag_outer_R)^(-chi/2)*besseli(n,radius/lambda)/besseli(n,diaphrag_outer_R/lambda);
end