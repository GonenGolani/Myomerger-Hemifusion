
clc
clear all
close all

main_folder='D:\Evolver_calc\Hemifission\Diaphragm\';
cd(main_folder)

kT=4.11*10^-21;
nm=10^-9;
kappa_m=10*kT;

J_0m=-0.1*nm^-1;

pore_line_tension=10*10^-12;
chi=-0.5;

%%
%%folder 1
folder1='D:\Evolver_calc\Hemifission\Diaphragm\l=1.5nm\R=free x=4 y=15 fix';
cd(folder1);

lambda=1.5*nm;

cont_J01=importdata('cont_J0.txt');
J01=importdata('J0.txt');
fit_radius1=importdata('fit_radius.txt');
radius1=importdata('Radius.txt');
fit_angle1=importdata('fit_angle.txt');
angle1=importdata('Angle.txt');
fit_energy1=importdata('fit_Energy.txt');
energy1=importdata('Energy.txt');
number_of_points_for_fit=length(cont_J01);
tilt_at_diaph=tan(pi*(45-fit_angle1/2)/180);

int=1;
while (int<=number_of_points_for_fit)
    tilt_0=tilt_at_diaph(int);
    diaphrag_outer_R=fit_radius1(int);
    [Barrier_energy00(int),Critical_Radius00(int)] = Pore_barrier_energy(lambda/nm,J_0m*nm,tilt_0,diaphrag_outer_R,0,kappa_m/kT,pore_line_tension/kT*nm);
    [Barrier_energy025(int),Critical_Radius025(int)] = Pore_barrier_energy(lambda/nm,J_0m*nm,tilt_0,diaphrag_outer_R,-0.25,kappa_m/kT,pore_line_tension/kT*nm);
    [Barrier_energy05(int),Critical_Radius05(int)] = Pore_barrier_energy(lambda/nm,J_0m*nm,tilt_0,diaphrag_outer_R,-0.5,kappa_m/kT,pore_line_tension/kT*nm);
    [Barrier_energy1(int),Critical_Radius1(int)] = Pore_barrier_energy(lambda/nm,J_0m*nm,tilt_0,diaphrag_outer_R,-1,kappa_m/kT,pore_line_tension/kT*nm);

    int=int+1;
end

Delta_E00=Barrier_energy00(1)-Barrier_energy00;
Delta_E025=Barrier_energy025(1)-Barrier_energy025;
Delta_E05=Barrier_energy05(1)-Barrier_energy05;
Delta_E1=Barrier_energy1(1)-Barrier_energy1;

Area_ratio1=(fit_radius1./fit_radius1(1)).^2;

%%

figure(1)
hold on
p1=plot(cont_J01,(Area_ratio1.*exp(Delta_E00)));
p2=plot(cont_J01,(Area_ratio1.*exp(Delta_E025)));
p3=plot(cont_J01,(Area_ratio1.*exp(Delta_E05)));
p4=plot(cont_J01,(Area_ratio1.*exp(Delta_E1)));

xL=xlabel('J0[nm^-^1]');
yL=ylabel('Acceleration factor');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white')
legend('chi=0','chi=-0.25','chi=-0.5','chi=-1')
p1.LineWidth=2; p1.Color='blue'; 
p2.LineWidth=2; p2.Color='red'; 
p3.LineWidth=2; p3.Color='green'; 
p4.LineWidth=2; p4.Color='magenta';

figure(2)
hold on
p1=plot(cont_J01,Barrier_energy00);
p2=plot(cont_J01,Barrier_energy025);
p3=plot(cont_J01,Barrier_energy05);
p4=plot(cont_J01,Barrier_energy1);

xL=xlabel('J0[nm^-^1]');
yL=ylabel('Barrier Energy [kb_bT]');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white')
legend('chi=0','chi=-0.25','chi=-0.5','chi=-1')
p1.LineWidth=2; p1.Color='blue'; 
p2.LineWidth=2; p2.Color='red'; 
p3.LineWidth=2; p3.Color='green'; 
p4.LineWidth=2; p4.Color='magenta';


cd('D:\');

