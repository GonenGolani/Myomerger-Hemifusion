function [Barrier_energy,Critical_Radius] = Pore_barrier_energy_no_second_min(lambda,J_0m,tilt_0,diaphrag_outer_R,chi,kappa_m,line_tension)
    Max_point_num=200;
    dr=diaphrag_outer_R/Max_point_num;
    radius=dr;
    Surface_Energy=0;
    int=1;
    while (radius<=diaphrag_outer_R)

        dE =energy_density_no_second_min(lambda,J_0m,tilt_0,diaphrag_outer_R,radius,chi,kappa_m);
        Surface_Energy=Surface_Energy+4*pi*dE*radius*dr;
        Line_energy=2*pi*radius*line_tension;
        Pore_energy(int)=Line_energy-Surface_Energy;
        radius_vec(int)=radius;
        radius=radius+dr;
        int=int+1;
    end
    
    [Barrier_energy,index]=max(Pore_energy);
    Critical_Radius=radius_vec(index);

end
