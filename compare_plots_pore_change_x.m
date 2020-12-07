
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
chi=0;

%%
%%folder 1
folder1='D:\Evolver_calc\Hemifission\Diaphragm\l=1.0nm\R=free x=3 y=15 fix';
cd(folder1);

lambda=1*nm;

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
    [Barrier_energy1(int),Critical_Radius1(int)] = Pore_barrier_energy(lambda/nm,J_0m*nm,tilt_0,diaphrag_outer_R,chi,kappa_m/kT,pore_line_tension/kT*nm);
    int=int+1;
end

Delta_E1=Barrier_energy1(1)-Barrier_energy1;
Area_ratio1=(fit_radius1./fit_radius1(1)).^2;


%%
%%folder 2
folder2='D:\Evolver_calc\Hemifission\Diaphragm\l=1.0nm\R=free x=4 y=15 fix';
cd(folder2);

lambda=1*nm;

cont_J02=importdata('cont_J0.txt');
J02=importdata('J0.txt');
fit_radius2=importdata('fit_radius.txt');
radius2=importdata('Radius.txt');
fit_angle2=importdata('fit_angle.txt');
angle2=importdata('Angle.txt');
fit_energy2=importdata('fit_Energy.txt');
energy2=importdata('Energy.txt');
number_of_points_for_fit=length(cont_J02);
tilt_at_diaph=tan(pi*(45-fit_angle2/2)/180);

int=1;
while (int<=number_of_points_for_fit)
    tilt_0=tilt_at_diaph(int);
    diaphrag_outer_R=fit_radius2(int);
    [Barrier_energy2(int),Critical_Radius2(int)] = Pore_barrier_energy(lambda/nm,J_0m*nm,tilt_0,diaphrag_outer_R,chi,kappa_m/kT,pore_line_tension/kT*nm);
    int=int+1;
end

Delta_E2=Barrier_energy2(1)-Barrier_energy2;
Area_ratio2=(fit_radius2./fit_radius2(1)).^2;
%%
%%
%%folder 3
folder3='D:\Evolver_calc\Hemifission\Diaphragm\l=1.0nm\R=free x=5 y=15 fix';
cd(folder3);

lambda=1*nm;

cont_J03=importdata('cont_J0.txt');
J03=importdata('J0.txt');
fit_radius3=importdata('fit_radius.txt');
radius3=importdata('Radius.txt');
fit_angle3=importdata('fit_angle.txt');
angle3=importdata('Angle.txt');
fit_energy3=importdata('fit_Energy.txt');
energy3=importdata('Energy.txt');
number_of_points_for_fit=length(cont_J03);
tilt_at_diaph=tan(pi*(45-fit_angle3/2)/180);

int=1;
while (int<=number_of_points_for_fit)
    tilt_0=tilt_at_diaph(int);
    diaphrag_outer_R=fit_radius3(int);
    [Barrier_energy3(int),Critical_Radius3(int)] = Pore_barrier_energy(lambda/nm,J_0m*nm,tilt_0,diaphrag_outer_R,chi,kappa_m/kT,pore_line_tension/kT*nm);
    int=int+1;
end

Delta_E3=Barrier_energy3(1)-Barrier_energy3;
Area_ratio3=(fit_radius3./fit_radius3(1)).^2;


%%

figure(1)
hold on
p1=plot(cont_J01,fit_radius1);
p2=plot(cont_J02,fit_radius2);
p3=plot(cont_J03,fit_radius3);

s1=scatter(J01,radius1);
s2=scatter(J02,radius2);
s3=scatter(J03,radius3);

xL=xlabel('J0[nm^-^1]');
yL=ylabel('Radius[nm]');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white')
legend('x=2nm','x=4nm','x=5nm')
p1.LineWidth=2; p1.Color='blue'; 
p2.LineWidth=2; p2.Color='red'; 
p3.LineWidth=2; p3.Color='green'; 
s1.MarkerEdgeColor='blue';  
s2.MarkerEdgeColor='red'; 
s3.MarkerEdgeColor='green';  


figure(2)
hold on
p1=plot(cont_J01,fit_angle1);
p2=plot(cont_J02,fit_angle2);
p3=plot(cont_J03,fit_angle3);
s1=scatter(J01,angle1);
s2=scatter(J02,angle2);
s3=scatter(J03,angle3);
xL=xlabel('J0[nm^-^1]');
yL=ylabel('Angle');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white')
legend('x=3nm','x=4nm','x=5nm')
p1.LineWidth=2; p1.Color='blue'; 
p2.LineWidth=2; p2.Color='red'; 
p3.LineWidth=2; p3.Color='green'; 
p4.LineWidth=2; p4.Color='magenta'; 
s1.MarkerEdgeColor='blue';  
s2.MarkerEdgeColor='red'; 
s3.MarkerEdgeColor='green';  
s4.MarkerEdgeColor='magenta'; 

figure(3)
hold on
p1=plot(cont_J01,fit_energy1*kappa_m/kT);
p2=plot(cont_J02,fit_energy2*kappa_m/kT);
p3=plot(cont_J03,fit_energy3*kappa_m/kT);
s1=scatter(J01,energy1*kappa_m/kT);
s2=scatter(J02,energy2*kappa_m/kT);
s3=scatter(J03,energy3*kappa_m/kT);

xL=xlabel('J0[nm^-^1]');
yL=ylabel('Energy');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white')
legend('x=2nm','x=4nm','x=5nm')
p1.LineWidth=2; p1.Color='blue'; 
p2.LineWidth=2; p2.Color='red'; 
p3.LineWidth=2; p3.Color='green'; 
s1.MarkerEdgeColor='blue';  
s2.MarkerEdgeColor='red'; 
s3.MarkerEdgeColor='green';  

figure(4)
hold on
p1=plot(cont_J01,Barrier_energy1);
p2=plot(cont_J02,Barrier_energy2);
p3=plot(cont_J03,Barrier_energy3);

xL=xlabel('J0[nm^-^1]');
yL=ylabel('Barrier Energy [kT]');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white')
legend('x=2nm','x=4nm','x=5nm')
p1.LineWidth=2; p1.Color='blue'; 
p2.LineWidth=2; p2.Color='red'; 
p3.LineWidth=2; p3.Color='green'; 
p4.LineWidth=2; p4.Color='magenta';


figure(5)
hold on
p1=plot(cont_J01,Critical_Radius1);
p2=plot(cont_J02,Critical_Radius2);
p3=plot(cont_J03,Critical_Radius3);

xL=xlabel('J0[nm^-^1]');
yL=ylabel('Pore Radius [nm]');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white')
legend('x=2nm','x=4nm','x=5nm')
p1.LineWidth=2; p1.Color='blue'; 
p2.LineWidth=2; p2.Color='red'; 
p3.LineWidth=2; p3.Color='green'; 
p4.LineWidth=2; p4.Color='magenta';


figure(6)
hold on
p1=plot(cont_J01,log(Area_ratio1.*exp(Delta_E1)));
p2=plot(cont_J02,log(Area_ratio2.*exp(Delta_E2)));
p3=plot(cont_J03,log(Area_ratio3.*exp(Delta_E3)));
xL=xlabel('J0[nm^-^1]');
yL=ylabel('ln(Acceleration factor)');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white')
legend('x=3nm','x=4nm','x=5nm')
p1.LineWidth=2; p1.Color='blue'; 
p2.LineWidth=2; p2.Color='red'; 
p3.LineWidth=2; p3.Color='green'; 
p4.LineWidth=2; p4.Color='magenta';


figure(7)
hold on
p1=plot(cont_J01,(Area_ratio1.*exp(Delta_E1)));
p2=plot(cont_J02,(Area_ratio2.*exp(Delta_E2)));
p3=plot(cont_J03,(Area_ratio3.*exp(Delta_E3)));


xL=xlabel('J0[nm^-^1]');
yL=ylabel('Acceleration factor');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white')
legend('x=2nm','x=4nm','x=5nm')
p1.LineWidth=2; p1.Color='blue'; 
p2.LineWidth=2; p2.Color='red'; 
p3.LineWidth=2; p3.Color='green'; 
p4.LineWidth=2; p4.Color='magenta';

cd('D:\');

