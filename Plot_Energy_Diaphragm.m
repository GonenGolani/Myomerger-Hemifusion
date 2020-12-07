
clc
clear all
close all


kT=4.11*10^-21;
nm=10^-9;
kappa_m=10*kT;
l=1*nm;
kappa_t_m=kappa_m/l^2;
Radius=1.25*nm;
back_J0=-0.1*nm^-1;

pore_line_tension=30*10^-12;


% folders name
main_folder_Diaphragm='D:\Evolver_calc\Hemifission\Diaphragm\R=3 x=4 y=15';
J0_max=0.100;



cd(main_folder_Diaphragm);

J0=-0.1;
int_J0=1;
int_all=1;
while (J0<=J0_max)

    J0_folder_name = sprintf('J0_up %.3f',J0);
    if exist(J0_folder_name)==7
        cd(J0_folder_name);
        theta_u=0;
        int_angle=1;
        while (theta_u<=90)
 
                Up_folder_name = sprintf('Up %.1f',theta_u);
                if exist(Up_folder_name)==7
                    cd(Up_folder_name);
                    if exist('Minimal_energy.txt')==2
                        energy=importdata('Minimal_energy.txt');
                        %sum_energy_down=importdata('sum_energy_down.txt');
                        %sum_energy_up=importdata('sum_energy_up.txt');
                        Diaphragm_energy(int_angle)=importdata('Diaphragm_energy.txt');
                        if (energy<40)
              

                            theta_u_vec(int_angle)=theta_u;
                            Energy_vec(int_angle)=energy;
                            %sum_energy_down_vec(int_angle)=sum_energy_down;
                            %sum_energy_up_vec(int_angle)=sum_energy_up;

                            int_angle=int_angle+1;
                        end
                        int_all=int_all+1;
                    end
                  cd(main_folder_Diaphragm);
                  cd(J0_folder_name);
                end
          theta_u=theta_u+0.1;
        end
        
    theta_u_vec__inter=linspace(theta_u_vec(1,1),theta_u_vec(length(theta_u_vec)),200); 
    P = polyfit(theta_u_vec,Energy_vec,2);
    %Linar_Diaph_energy = P(1)*theta_u_vec__inter.^4+P(2)*theta_u_vec__inter.^3+P(3)*theta_u_vec__inter.^2+P(4)*theta_u_vec__inter+P(5);
    %Linar_Diaph_energy = P(1)*theta_u_vec__inter.^3+P(2)*theta_u_vec__inter.^2+P(3)*theta_u_vec__inter.^1+P(4);
    Linar_Diaph_energy = P(1)*theta_u_vec__inter.^2+P(2)*theta_u_vec__inter.^1+P(3);

    [minimal_energy(int_J0),angle_index]=min(Linar_Diaph_energy);
    minimal_energy_theta_Diaph(int_J0)=theta_u_vec__inter(angle_index);   
        

     [delta_angle,Theta_index]=min(abs(theta_u_vec-minimal_energy_theta_Diaph(int_J0)));
    
    %minimal_sum_energy_down_Diaph(int_J0)=sum_energy_down_vec(Theta_index);
    %minimal_sum_energy_up_Diaph(int_J0)=sum_energy_up_vec(Theta_index);
    %minimal_avg_stress_down_Diaph(int_J0)=avg_stress_down_vec(Theta_index);
    %minimal_avg_stress_up_Diaph(int_J0)=avg_stress_up_vec(Theta_index);
    Diaphragm_energy_min(int_J0)=Diaphragm_energy(Theta_index);

  
    clear Energy_vec;
    clear theta_u_vec;
    J0_vec_diaph(int_J0)=J0;
    int_J0=int_J0+1;
    cd(main_folder_Diaphragm);
    end

    J0=J0+0.001
end



J_0_vec_inter_daiph=linspace(J0_vec_diaph(1,1),J0_vec_diaph(length(J0_vec_diaph)),200); 

P3=polyfit(J0_vec_diaph,minimal_energy_theta_Diaph,1);
%Angle_fit=P3(1)*J_0_vec_inter_daiph.^3+P3(2)*J_0_vec_inter_daiph.^2+P3(3)*J_0_vec_inter_daiph+P3(4);
%Angle_fit=P3(1)*J_0_vec_inter_daiph.^2+P3(2)*J_0_vec_inter_daiph+P3(3);
Angle_fit=P3(1)*J_0_vec_inter_daiph+P3(2);


figure (1);
scatter((J0_vec_diaph),minimal_energy_theta_Diaph);
hold on
p1=plot((J_0_vec_inter_daiph),Angle_fit,'--');
p1.LineWidth=2; p1.Color='red';
xlabel('J0[nm^-^1]');
ylabel('Angle')
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white');


figure (4);
s1=scatter(J0_vec_diaph,(minimal_energy)*kappa_m/kT,'b');
p1.LineWidth=2; p1.Color='red';

xL=xlabel('J0[nm^-^1]');
yL=ylabel('Energy[k_bT]');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white');



tilt_angle=pi/180*(90-Angle_fit)/2;
tilt_vec=tan(tilt_angle);
delta_Stress_to_stress=(tilt_vec/tilt_vec(1)).^2-1;
Stress_in_middle=kappa_t_m*(tilt_vec(1)^2*1/(besseli(1,Radius/l)).^2-2*back_J0*l*tilt_vec(1)/besseli(1,Radius/l));
Stress_at_edge=kappa_t_m*tilt_vec(1)^2*(besseli(0,Radius/l)^2/besseli(1,Radius/l)^2+1)-2*kappa_t_m*tilt_vec(1)*back_J0*l*besseli(0,Radius/l)/besseli(1,Radius/l);


prefactor=pi*pore_line_tension^2/(Stress_in_middle*kT);



fusion_rate_exponent=prefactor*delta_Stress_to_stress;

figure(5)
p1=plot(J_0_vec_inter_daiph,exp(fusion_rate_exponent));
p1.LineWidth=2; p1.Color='red';

xL=xlabel('J0 [nm^-^1]');
yL=ylabel('Rate factor');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white');


Diaphragm_bilayer_energy=pi*tilt_vec.^2*kappa_m*(Radius/l)^2*(besseli(0,Radius/l)^2-besseli(0,Radius/l)*besseli(2,Radius/l))/besseli(1,Radius/l)^2;
figure(6)
p1=plot(J_0_vec_inter_daiph,Diaphragm_bilayer_energy/(pi*Radius^2)*1000);
p1.LineWidth=2; p1.Color='red';
xL=xlabel('J0 [nm^-^1]');
yL=ylabel('Avg Stress [mN/m]');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white');

cd('D:\');
