
clc
clear all
close all

main_folder='D:\Evolver_calc\Hemifission\Diaphragm\l=1.0nm\R=free x=5 y=15 fix';
cd(main_folder)

kT=4.11*10^-21;
nm=10^-9;
kappa_m=10*kT;
Chi=-0.5;
%
l=2*nm;
energy_limit_up=10;
energy_limit_down=-50;


%
kappa_t_m=kappa_m/l^2;
back_J0=-0.1*nm^-1;

pore_line_tension=10*10^-12;
Radius_max=3;
Radius_min=0.95;
Radius_res=0.01;
% folders name
J0_res=0.001;
J0_max=-0.065;
J0=-0.1;
J0_int=1;


number_of_points_for_fit=300;



while (J0<=J0_max)
    main_folder_Diaphragm=sprintf('J0_up %.3f',J0);
    if exist(main_folder_Diaphragm)==7
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
                                % this to correct worng line in Evolver
                                t0=tan(pi/4-theta_u*pi/180/2);
                                correct_energy=4*pi*Radius*back_J0*nm*t0;
                                %
                                if (energy+correct_energy<energy_limit_up && energy+correct_energy>energy_limit_down)

                                    theta_u_vec(Radius_int,int_angle)=theta_u;
                                    Energy_vec(Radius_int,int_angle)=energy+correct_energy;

                                    int_angle=int_angle+1;
                                end
                            end
                          cd(main_folder);
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
            figure (J0_int);
            scatter3(Radius*length_this,theta_u_vec(Radius_int,:),Energy_vec(Radius_int,:)*kappa_m/kT,'r');
            hold on
            p1=plot3(Radius*ones(1,number_of_points_for_fit),theta_u_vec__inter(Radius_int,:),fit_angle_energy(Radius_int,:)*kappa_m/kT,'g');
            s3=scatter3(Radius,minimal_angle_min(Radius_int), minimal_energy_angle_min(Radius_int)*kappa_m/kT);
            s3.MarkerFaceColor='blue'; s3.LineWidth=0.5;
            title(main_folder_Diaphragm)
            xlabel('Radius');
            ylabel('Angle');
            zlabel('Energy');
            % plot data end

            clear Energy_vec;   clear theta_u_vec; clear P_angle; 
            Radius_vec_diaph(Radius_int)=Radius;
            Radius_int=Radius_int+1;
            cd(main_folder);
            cd(main_folder_Diaphragm);
            end

            Radius=Radius+Radius_res;
        end

        %fit minima to find gobal minima with respect to radius

        %min radius
        Radius_inter_vec=linspace(min(Radius_vec_diaph),max(Radius_vec_diaph),number_of_points_for_fit);
        P_radius = polyfit(Radius_vec_diaph,minimal_energy_angle_min,2);
        fit_radius_energy = P_radius(1)*Radius_inter_vec.^2+P_radius(2)*Radius_inter_vec.^1+P_radius(3);

        [min_energy_radius,min_radius_index]=min(fit_radius_energy);
        %min angle
        P_radius_angle = polyfit(Radius_vec_diaph,minimal_angle_min,2);
        fit_angle = P_radius_angle(1)*Radius_inter_vec.^2+P_radius_angle(2)*Radius_inter_vec+P_radius_angle(3);
        

        
        min_radius_fit(J0_int)=Radius_inter_vec(min_radius_index);    
        min_angle_fit(J0_int)=fit_angle(min_radius_index);
        min_energy(J0_int)=min_energy_radius;
        
        %plot fit data
        figure (J0_int);
        p2=plot3(Radius_inter_vec,fit_angle,fit_radius_energy*kappa_m/kT);
        p2.LineWidth=4; p2.Color='black'; p2.LineStyle='--';
        s4=scatter3(min_radius_fit(J0_int),min_angle_fit(J0_int), min_energy(J0_int)*kappa_m/kT);
        s4.MarkerFaceColor='cyan';  s4.LineWidth=5;
        J0_vec(J0_int)=J0;
        J0_int=J0_int+1;
        clear Radius_inter_vec;clear Radius_vec_diaph;  clear minimal_energy_angle_min; clear minimal_angle_min; clear fit_angle; clear fit_radius_energy; clear fit_angle_energy;
        cd(main_folder)

    end
    J0=J0+J0_res; 
end


cd(main_folder)
J0_inter_vec=linspace(min(J0_vec),max(J0_vec),number_of_points_for_fit);
fit_min_radius=polyfit(J0_vec,min_radius_fit,1);
cont_radius=fit_min_radius(1)*J0_inter_vec+fit_min_radius(2);

fit_min_angle=polyfit(J0_vec,min_angle_fit,1);
cont_angle=fit_min_angle(1)*J0_inter_vec+fit_min_angle(2);

fit_min_energy=polyfit(J0_vec,min_energy,1);
conr_energy=fit_min_energy(1)*J0_inter_vec+fit_min_energy(2);

save('cont_J0.txt','J0_inter_vec','-ascii');
save('J0.txt','J0_vec','-ascii');

save('fit_radius.txt','cont_radius','-ascii');
save('Radius.txt','min_radius_fit','-ascii');

save('fit_angle.txt','cont_angle','-ascii');
save('Angle.txt','min_angle_fit','-ascii');

save('fit_Energy.txt','conr_energy','-ascii');
save('Energy.txt','min_energy','-ascii');

save('Radius_fit_limit.txt','Radius_min','Radius_max','-ascii');
save('energy_fit_limit.txt','energy_limit_down','energy_limit_up','-ascii');



figure(J0_int+1)
scatter(J0_vec,min_radius_fit)
hold on
plot(J0_inter_vec,cont_radius)
xL=xlabel('J0[nm^-^1]');
yL=ylabel('Radius[nm]');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white')

figure(J0_int+2)
scatter(J0_vec,min_angle_fit)
hold on
plot(J0_inter_vec,cont_angle)
xL=xlabel('J0_vec[nm]');
yL=ylabel('Angle');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white')

figure(J0_int+3)
scatter(J0_vec,min_energy*kappa_m/kT)
hold on
plot(J0_inter_vec,conr_energy*kappa_m/kT)
xL=xlabel('J0[nm^-^1]');
yL=ylabel('Energy');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white')

%find energy barrier change
int=1;
while (int<=number_of_points_for_fit)
    rho_vec=linspace(0,cont_radius(int)*nm,number_of_points_for_fit);
    rhol=rho_vec/l;
    t0_to_besseli=tan(pi*(45-cont_angle(int)/2)/180)/besseli(1,cont_radius(int));
    pore_energy_LT=2*pi*pore_line_tension*rho_vec;
    pore_energy_kappa_1=-2*pi*kappa_m.*rhol.*besseli(1,rhol).*besseli(0,rhol)*t0_to_besseli^2;
    Pore_energy_kappa_2=4.*pi.*kappa_m*back_J0*rho_vec.*besseli(1,rhol)*t0_to_besseli;
    pore_energy=pore_energy_LT+pore_energy_kappa_1+Pore_energy_kappa_2;
    [max_pore_energy(int),index]=max(pore_energy);
    max_pore_radius(int)=rho_vec(index);
    int=int+1;
end

figure(J0_int+4)
plot(J0_inter_vec,max_pore_energy/kT)
xL=xlabel('J0[nm^-^1]');
yL=ylabel('Barrier Energy [kT]');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white')

figure(J0_int+5)
plot(J0_inter_vec,max_pore_radius/nm)
xL=xlabel('J0[nm^-^1]');
yL=ylabel('Pore Radius [nm]');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white')

Delta_E=max_pore_energy-max_pore_energy(1);

figure(J0_int+6)
plot(J0_inter_vec,exp(-Delta_E/kT))
xL=xlabel('J0[nm^-^1]');
yL=ylabel('Acceleration factor');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white')

cd('D:\');
