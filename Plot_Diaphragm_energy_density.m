
clc
clear all
close all

kT=4.11*10^-21;
nm=10^-9;
kappa_m=10*kT;
l=2*nm;
kappa_t_m=kappa_m/l^2;
Radi0=1.25*nm;



Angle1=28.47;
Angle2=26.82;

Radius=linspace(0,Radi0,100);

Tilt1=tan(pi/180*(90-Angle1)/2);
Tilt2=tan(pi/180*(90-Angle2)/2);


Energy_density_no_tilt=kappa_t_m*((besseli(0,Radius/l)).^2+(besseli(1,Radius/l)).^2)./(besseli(1,Radi0/l)).^2;

figure(1)
plot(Radius/nm,Tilt1^2*Energy_density_no_tilt*1000)
hold on
%plot(Radius/nm,Tilt2^2*Energy_density_no_tilt*1000)
xL=xlabel('Radius[nm]');
yL=ylabel('Stress[mN/m]');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white');
legend('0%','10%')

figure(2)
plot(Radius/nm,Tilt2^2*Energy_density_no_tilt*1000-Tilt1^2*Energy_density_no_tilt*1000)
xL=xlabel('Radius[nm]');
yL=ylabel('Delta Stress[mN/m]');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white');




cd('D:\');
