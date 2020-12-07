
clc
clear all
close all

main_folder='D:\Evolver_calc\Hemifission\Diaphragm\';
cd(main_folder)

kT=4.11*10^-21;
nm=10^-9;
kappa_m=10*kT;

back_J0=-0.1*nm^-1;

pore_line_tension=10*10^-12;
%%
%%folder 1
folder1='D:\Evolver_calc\Hemifission\Diaphragm\l=1.0nm\R=free x=4 y=10';
cd(folder1);

l=1*nm;
kappa_t_m=kappa_m/l^2;

cont_J01=importdata('cont_J0.txt');
J01=importdata('J0.txt');
fit_radius1=importdata('fit_radius.txt');
radius1=importdata('Radius.txt');
fit_angle1=importdata('fit_angle.txt');
angle1=importdata('Angle.txt');
fit_energy1=importdata('fit_Energy.txt');
energy1=importdata('Energy.txt');
number_of_points_for_fit=length(cont_J01);

int=1;
while (int<=number_of_points_for_fit)
    rho_vec=linspace(0,fit_radius1(int)*nm,number_of_points_for_fit);
    rhol=rho_vec/l;
    %
    t0_to_besseli=tan(pi*fit_angle1(int)/180)/besseli(1,fit_radius1(int));
    %
    pore_energy_LT=2*pi*pore_line_tension*rho_vec;
    pore_energy_kappa_1=-2*pi*kappa_m.*rhol.*besseli(1,rhol).*besseli(0,rhol)*t0_to_besseli^2;
    Pore_energy_kappa_2=4.*pi.*kappa_m*back_J0*rho_vec.*besseli(1,rhol)*t0_to_besseli;
    
    pore_energy=pore_energy_LT+pore_energy_kappa_1+Pore_energy_kappa_2;
    [max_pore_energy1(int),index]=max(pore_energy);
    max_pore_radius1(int)=rho_vec(index);
    int=int+1;
end

Delta_E1=max_pore_energy1(1)-max_pore_energy1;
Area_ratio1=(fit_radius1./fit_radius1(1)).^2;

%%

%%
%%folder 2
folder2='D:\Evolver_calc\Hemifission\Diaphragm\l=2.0nm\R=free x=3 y=15';
cd(folder2);

l=2*nm;
kappa_t_m=kappa_m/l^2;

cont_J02=importdata('cont_J0.txt');
J02=importdata('J0.txt');
fit_radius2=importdata('fit_radius.txt');
radius2=importdata('Radius.txt');
fit_angle2=importdata('fit_angle.txt');
angle2=importdata('Angle.txt');
fit_energy2=importdata('fit_Energy.txt');
energy2=importdata('Energy.txt');
number_of_points_for_fit=length(cont_J02);

int=1;
while (int<=number_of_points_for_fit)
    rho_vec=linspace(0,fit_radius2(int)*nm,number_of_points_for_fit);
    rhol=rho_vec/l;
    %
    t0_to_besseli=tan(pi*fit_angle2(int)/180)/besseli(1,fit_radius2(int));
    %
    pore_energy_LT=2*pi*pore_line_tension*rho_vec;
    pore_energy_kappa_1=-2*pi*kappa_m.*rhol.*besseli(1,rhol).*besseli(0,rhol)*t0_to_besseli^2;
    Pore_energy_kappa_2=4.*pi.*kappa_m*back_J0*rho_vec.*besseli(1,rhol)*t0_to_besseli;
    
    pore_energy=pore_energy_LT+pore_energy_kappa_1+Pore_energy_kappa_2;
    [max_pore_energy2(int),index]=max(pore_energy);
    max_pore_radius2(int)=rho_vec(index);
    int=int+1;
end

Delta_E2=max_pore_energy2(1)-max_pore_energy2;
Area_ratio2=(fit_radius2./fit_radius2(1)).^2;


%%
%%
%%folder 3
folder3='D:\Evolver_calc\Hemifission\Diaphragm\l=2.0nm\R=free x=4 y=15';
cd(folder3);

l=2*nm;
kappa_t_m=kappa_m/l^2;

cont_J03=importdata('cont_J0.txt');
J03=importdata('J0.txt');
fit_radius3=importdata('fit_radius.txt');
radius3=importdata('Radius.txt');
fit_angle3=importdata('fit_angle.txt');
angle3=importdata('Angle.txt');
fit_energy3=importdata('fit_Energy.txt');
energy3=importdata('Energy.txt');
number_of_points_for_fit=length(cont_J03);

int=1;
while (int<=number_of_points_for_fit)
    rho_vec=linspace(0,fit_radius3(int)*nm,number_of_points_for_fit);
    rhol=rho_vec/l;
    %
    t0_to_besseli=tan(pi*fit_angle3(int)/180)/besseli(1,fit_radius3(int));
    %
    pore_energy_LT=2*pi*pore_line_tension*rho_vec;
    pore_energy_kappa_1=-2*pi*kappa_m.*rhol.*besseli(1,rhol).*besseli(0,rhol)*t0_to_besseli^2;
    Pore_energy_kappa_2=4.*pi.*kappa_m*back_J0*rho_vec.*besseli(1,rhol)*t0_to_besseli;
    
    pore_energy=pore_energy_LT+pore_energy_kappa_1+Pore_energy_kappa_2;
    [max_pore_energy3(int),index]=max(pore_energy);
    max_pore_radius3(int)=rho_vec(index);
    int=int+1;
end

Delta_E3=max_pore_energy3(1)-max_pore_energy3;
Area_ratio3=(fit_radius3./fit_radius3(1)).^2;

%%

%%
%%folder 4
folder4='D:\Evolver_calc\Hemifission\Diaphragm\l=2.0nm\R=free x=5 y=15';
cd(folder4);

l=2*nm;
kappa_t_m=kappa_m/l^2;

cont_J04=importdata('cont_J0.txt');
J04=importdata('J0.txt');
fit_radius4=importdata('fit_radius.txt');
radius4=importdata('Radius.txt');
fit_angle4=importdata('fit_angle.txt');
angle4=importdata('Angle.txt');
fit_energy4=importdata('fit_Energy.txt');
energy4=importdata('Energy.txt');
number_of_points_for_fit=length(cont_J04);

int=1;
while (int<=number_of_points_for_fit)
    rho_vec=linspace(0,fit_radius4(int)*nm,number_of_points_for_fit);
    rhol=rho_vec/l;
    %
    t0_to_besseli=tan(pi*fit_angle4(int)/180)/besseli(1,fit_radius4(int));
    %
    pore_energy_LT=2*pi*pore_line_tension*rho_vec;
    pore_energy_kappa_1=-2*pi*kappa_m.*rhol.*besseli(1,rhol).*besseli(0,rhol)*t0_to_besseli^2;
    Pore_energy_kappa_2=4.*pi.*kappa_m*back_J0*rho_vec.*besseli(1,rhol)*t0_to_besseli;
    
    pore_energy=pore_energy_LT+pore_energy_kappa_1+Pore_energy_kappa_2;
    [max_pore_energy4(int),index]=max(pore_energy);
    max_pore_radius4(int)=rho_vec(index);
    int=int+1;
end

Delta_E4=max_pore_energy4(1)-max_pore_energy4;
Area_ratio4=(fit_radius4./fit_radius4(1)).^2;

%%



figure(1)
hold on
%p1=plot(cont_J01,fit_radius1);
p2=plot(cont_J02,fit_radius2);
p3=plot(cont_J03,fit_radius3);
p4=plot(cont_J04,fit_radius4);

%s1=scatter(J01,radius1);
s2=scatter(J02,radius2);
s3=scatter(J03,radius3);
s4=scatter(J04,radius4);

xL=xlabel('J0[nm^-^1]');
yL=ylabel('Radius[nm]');
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


figure(2)
hold on
%p1=plot(cont_J01,fit_angle1);
p2=plot(cont_J02,fit_angle2);
p3=plot(cont_J03,fit_angle3);
p4=plot(cont_J04,fit_angle4)
%s1=scatter(J01,angle1);
s2=scatter(J02,angle2);
s3=scatter(J03,angle3);
s4=scatter(J04,angle4);
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
%p1=plot(cont_J01,fit_energy1*kappa_m/kT);
p2=plot(cont_J02,fit_energy2*kappa_m/kT);
p3=plot(cont_J03,fit_energy3*kappa_m/kT);
p4=plot(cont_J04,fit_energy4*kappa_m/kT);
%s1=scatter(J01,energy1*kappa_m/kT);
s2=scatter(J02,energy2*kappa_m/kT);
s3=scatter(J03,energy3*kappa_m/kT);
s4=scatter(J04,energy4*kappa_m/kT);

xL=xlabel('J0[nm^-^1]');
yL=ylabel('Energy');
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

figure(4)
hold on
%p1=plot(cont_J01,max_pore_energy1/kT);
p2=plot(cont_J02,max_pore_energy2/kT);
p3=plot(cont_J03,max_pore_energy3/kT);
p4=plot(cont_J04,max_pore_energy4/kT);

xL=xlabel('J0[nm^-^1]');
yL=ylabel('Barrier Energy [kT]');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white')
legend('x=3nm','x=4nm','x=5nm')
p1.LineWidth=2; p1.Color='blue'; 
p2.LineWidth=2; p2.Color='red'; 
p3.LineWidth=2; p3.Color='green'; 
p4.LineWidth=2; p4.Color='magenta'


figure(5)
hold on
%p1=plot(cont_J01,max_pore_radius1/nm);
p2=plot(cont_J02,max_pore_radius2/nm);
p3=plot(cont_J03,max_pore_radius3/nm);
p4=plot(cont_J04,max_pore_radius4/nm);

xL=xlabel('J0[nm^-^1]');
yL=ylabel('Pore Radius [nm]');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white')
legend('x=3nm','x=4nm','x=5nm')
p1.LineWidth=2; p1.Color='blue'; 
p2.LineWidth=2; p2.Color='red'; 
p3.LineWidth=2; p3.Color='green'; 
p4.LineWidth=2; p4.Color='magenta'


figure(6)
hold on
%p1=plot(cont_J01,Area_ratio1.*exp(Delta_E1/kT));
p2=plot(cont_J02,Area_ratio2.*exp(Delta_E2/kT));
p3=plot(cont_J03,Area_ratio3.*exp(Delta_E3/kT));
p4=plot(cont_J04,Area_ratio4.*exp(Delta_E4/kT));
xL=xlabel('J0[nm^-^1]');
yL=ylabel('Acceleration factor');
xL.FontSize=14;
xL.FontWeight='bold';
yL.FontSize=14;
yL.FontWeight='bold';
set(gcf,'color','white')
legend('x=3nm','x=4nm','x=5nm')
p1.LineWidth=2; p1.Color='blue'; 
p2.LineWidth=2; p2.Color='red'; 
p3.LineWidth=2; p3.Color='green'; 
p4.LineWidth=2; p4.Color='magenta'

figure(7)
hold on
%p1=plot(cont_J01,Area_ratio1.*exp(Delta_E1/kT));
p2=plot(cont_J02,log(Area_ratio2.*exp(Delta_E2/kT)));
p3=plot(cont_J03,log(Area_ratio3.*exp(Delta_E3/kT)));
p4=plot(cont_J04,log(Area_ratio4.*exp(Delta_E4/kT)));
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
p4.LineWidth=2; p4.Color='magenta'
cd('D:\');
