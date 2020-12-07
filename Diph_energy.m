function [outputArg1] = Diph_energy(lambda,J_0m,tilt_0,diaphrag_outer_R,chi,kappa_m)
    Max_point_num=500;
    dr=diaphrag_outer_R/Max_point_num;
    radius=dr;
    Energy=0;
    while (radius<=diaphrag_outer_R)

        dE = energy_density(lambda,J_0m,tilt_0,diaphrag_outer_R,radius,chi,kappa_m);
        Energy=Energy+4*pi*dE*radius*dr;
        radius=radius+dr;
    end
    outputArg1=Energy;

end
