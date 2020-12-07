
clc
clear all
close all


kT=4.11*10^-21;
nm=10^-9;
kappa_m=10*kT;
%
l=2.5*nm;
%
kappa_t_m=kappa_m/l^2;
back_J0=-0.1*nm^-1;

pore_line_tension=30*10^-12;
Radius_max=8;
Radius_min=1.05;
Radius_res=0.01;
% folders name
main_folder_Diaphragm='D:\Evolver_calc\Hemifission\Diaphragm\l=2.0nm\R=free x=5 y=15 fix\J0_up -0.097';
J0_max=0.100;
number_of_points_for_fit=300;



cd(main_folder_Diaphragm);

Radius=Radius_min;
Radius_int=1;
while (Radius<=Radius_max)

    Radius_folder_name = sprintf('Radius %.2f',Radius);
    if exist(Radius_folder_name)==7
        cd(Radius_folder_name);
        theta_u=0;
        int_angle=1;
        while (theta_u<=90)
 
                Up_folder_name = sprintf('Up %.1f',theta_u);
                if exist(Up_folder_name)==7
                    cd(Up_folder_name);
                    if exist('Minimal_energy.txt')==2
                        energy=importdata('Minimal_energy.txt');
                        if (energy<20)
              
                            theta_u_vec(Radius_int,int_angle)=theta_u;
                            Energy_vec(Radius_int,int_angle)=energy;

                            int_angle=int_angle+1;
                        end
                    end
                  cd(main_folder_Diaphragm);
                  cd(Radius_folder_name);
                end
          theta_u=theta_u+0.1;
        end
    
    %fit energy for single Radius    
    theta_u_vec__inter(Radius_int,:)=linspace(theta_u_vec(Radius_int,1),theta_u_vec(Radius_int,length(theta_u_vec(Radius_int,:))),number_of_points_for_fit); 
    P_angle = polyfit(theta_u_vec(Radius_int,:),Energy_vec(Radius_int,:),2);
    fit_angle_energy(Radius_int,:) = P_angle(1)*theta_u_vec__inter(Radius_int,:).^2+P_angle(2)*theta_u_vec__inter(Radius_int,:).^1+P_angle(3);

    [minimal_energy_angle_min(Radius_int),angle_index]=min(fit_angle_energy(Radius_int,:));
    minimal_angle_min(Radius_int)=theta_u_vec__inter(Radius_int,angle_index);  
    

    
     [delta_angle,Theta_index]=min(abs(theta_u_vec-minimal_angle_min(Radius_int)));
    %plot data
    length_this=ones(1,length(Energy_vec(Radius_int,:)));
    figure (1);
    scatter3(Radius*length_this,theta_u_vec(Radius_int,:),Energy_vec(Radius_int,:)*kappa_m/kT,'r');
    hold on
    p1=plot3(Radius*ones(1,number_of_points_for_fit),theta_u_vec__inter(Radius_int,:),fit_angle_energy(Radius_int,:)*kappa_m/kT,'g');
    s3=scatter3(Radius,minimal_angle_min(Radius_int), minimal_energy_angle_min(Radius_int)*kappa_m/kT);
    s3.MarkerFaceColor='cyan'; s3.LineWidth=0.5; 
    % plot data end
    
    clear Energy_vec;
    clear theta_u_vec;
    Radius_vec_diaph(Radius_int)=Radius;
    Radius_int=Radius_int+1;
    cd(main_folder_Diaphragm);
    end

    Radius=Radius+Radius_res;
end

Radius_inter_vec=linspace(min(Radius_vec_diaph),max(Radius_vec_diaph),number_of_points_for_fit);
P_radius = polyfit(Radius_vec_diaph,minimal_energy_angle_min,2);
fit_radius_energy = P_radius(1)*Radius_inter_vec.^2+P_radius(2)*Radius_inter_vec.^1+P_radius(3);

[min_energy_radius,min_radius_index]=min(fit_radius_energy);
%min angle
P_radius_angle = polyfit(Radius_vec_diaph,minimal_angle_min,2);
fit_angle = P_radius_angle(1)*Radius_inter_vec.^2+P_radius_angle(2)*Radius_inter_vec+P_radius_angle(3);
        

        
min_radius_fit=Radius_inter_vec(min_radius_index);    
min_angle_fit=fit_angle(min_radius_index);
min_energy=min_energy_radius;
        
%plot fit data
figure (1);
p2=plot3(Radius_inter_vec,fit_angle,fit_radius_energy*kappa_m/kT);
p2.LineWidth=4; p2.Color='black'; p2.LineStyle='--';
s4=scatter3(min_radius_fit,min_angle_fit, min_energy*kappa_m/kT);
s4.MarkerFaceColor='cyan';  s4.LineWidth=5;
xL=xlabel('Radius[nm]');
yL=ylabel('Angle');
zL=zlabel('Energy[k_bT]');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
zL.FontSize=14;
zL.FontWeight='bold';
set(gcf,'color','white')
cd('D:\');
