%Testing functions

close all
clear all

kT=4.11*10^-21;
nm=10^-9;
kappa_m=10;
diaphrag_outer_R=3;
tilt_0=1;
dr=0.001;
J_0m=-0.1;

lambda=2;
pore_line_tension=10*10^-12/kT*nm;


int=1;
radius=dr;
while (radius<=diaphrag_outer_R)
	Delta_splay_vec_chi0(int) = Delta_splay(lambda,J_0m,tilt_0,diaphrag_outer_R,radius,0);
   	Delta_splay_vec_chi01(int) = Delta_splay(lambda,J_0m,tilt_0,diaphrag_outer_R,radius,-0.1);
   	Delta_splay_vec_chi05(int) = Delta_splay(lambda,J_0m,tilt_0,diaphrag_outer_R,radius,-0.5);

    Saddle_splay_vec_chi0(int) = Saddle_Splay(lambda,tilt_0,diaphrag_outer_R,radius,0);
    Saddle_splay_vec_chi01(int) = Saddle_Splay(lambda,tilt_0,diaphrag_outer_R,radius,-0.1);
    Saddle_splay_vec_chi05(int) = Saddle_Splay(lambda,tilt_0,diaphrag_outer_R,radius,-0.5);
   
    Tilt_vec_chi0(int) = Tilt(lambda,tilt_0,diaphrag_outer_R,radius,0);
    Tilt_vec_chi01(int) = Tilt(lambda,tilt_0,diaphrag_outer_R,radius,-0.1);
    Tilt_vec_chi05(int) = Tilt(lambda,tilt_0,diaphrag_outer_R,radius,-0.5);
    
    energy_density_chi0(int) = energy_density_no_second_min(lambda,J_0m,tilt_0,diaphrag_outer_R,radius,0,kappa_m);
    energy_density_chi01(int) = energy_density_no_second_min(lambda,J_0m,tilt_0,diaphrag_outer_R,radius,-0.1,kappa_m);
	energy_density_chi05(int) = energy_density_no_second_min(lambda,J_0m,tilt_0,diaphrag_outer_R,radius,-0.5,kappa_m);
    
	energy_density_chi05_no_second_min(int) = energy_density_no_second_min(lambda,J_0m,tilt_0,diaphrag_outer_R,radius,-0.5,kappa_m);
    
    radius_vec(int)=radius;    
	radius=radius+dr;
    int=int+1;
end


%%
int=2;
radius=dr;
Surface_integral_chi05_R1(1)=4*pi*radius*energy_density(lambda,J_0m,tilt_0,1,radius,-0.5,kappa_m);
radius_vec_1(1)=dr;

while (radius<=1)
    radius_vec_1(int)=radius;
    Surface_integral_chi05_R1(int)=Surface_integral_chi05_R1(int-1)+4*pi*radius*energy_density(lambda,J_0m,tilt_0,1,radius,-0.5,kappa_m);
	radius=radius+dr;
    int=int+1;
end

Barrier_energy_integral_chi05_R1=2*pi*radius_vec_1*pore_line_tension-Surface_integral_chi05_R1;



%%
int=2;
radius=dr;
Surface_integral_chi05_R2(1)=4*pi*radius*energy_density(lambda,J_0m,tilt_0,2,radius,-0.5,kappa_m);
radius_vec_2(1)=dr;

while (radius<=2)
    radius_vec_2(int)=radius;
    Surface_integral_chi05_R2(int)=Surface_integral_chi05_R2(int-1)+4*pi*radius*energy_density(lambda,J_0m,tilt_0,2,radius,-0.5,kappa_m);
	radius=radius+dr;
    int=int+1;
end

Barrier_energy_integral_chi05_R2=2*pi*radius_vec_2*pore_line_tension-Surface_integral_chi05_R2;


%%
int=2;
radius=dr;
Surface_integral_chi0_R3(1)=4*pi*radius*energy_density(lambda,J_0m,tilt_0,diaphrag_outer_R,radius,0,kappa_m);
Surface_integral_chi05_R3(1)=4*pi*radius*energy_density(lambda,J_0m,tilt_0,diaphrag_outer_R,radius,-0.5,kappa_m);
Surface_integral_chi1_R3(1)=4*pi*radius*energy_density(lambda,J_0m,tilt_0,diaphrag_outer_R,radius,-1,kappa_m);
radius_vec_3(1)=dr;

while (radius<=diaphrag_outer_R)
    radius_vec_3(int)=radius;
    Surface_integral_chi0_R3(int)=Surface_integral_chi0_R3(int-1)+4*pi*radius*energy_density(lambda,J_0m,tilt_0,diaphrag_outer_R,radius,0,kappa_m);
    Surface_integral_chi05_R3(int)=Surface_integral_chi05_R3(int-1)+4*pi*radius*energy_density(lambda,J_0m,tilt_0,diaphrag_outer_R,radius,-0.5,kappa_m);
    Surface_integral_chi1_R3(int)=Surface_integral_chi1_R3(int-1)+4*pi*radius*energy_density(lambda,J_0m,tilt_0,diaphrag_outer_R,radius,-1,kappa_m);

	radius=radius+dr;
    int=int+1;
end

Barrier_energy_integral_chi0_R3=2*pi*radius_vec_3*pore_line_tension-Surface_integral_chi0_R3;
Barrier_energy_integral_chi05_R3=2*pi*radius_vec_3*pore_line_tension-Surface_integral_chi05_R3;
Barrier_energy_integral_chi1_R3=2*pi*radius_vec_3*pore_line_tension-Surface_integral_chi1_R3;
%%

chi=0.0;
int=1;
while (chi>=-1)
	Diph_energy_vec_R1(int) = Diph_energy(lambda,J_0m,tilt_0,1,chi,kappa_m);
   	Diph_energy_vec_R2(int) = Diph_energy(lambda,J_0m,tilt_0,2,chi,kappa_m);
   	Diph_energy_vec_R3(int) = Diph_energy(lambda,J_0m,tilt_0,3,chi,kappa_m);
    
   	Diph_energy_vec(int) = Diph_energy(lambda,J_0m,tilt_0,3,chi,kappa_m);
   	Diph_energy_vec_no_second_min(int) = Diph_energy_no_second_min(lambda,J_0m,tilt_0,3,chi,kappa_m);

    Pore_barrier_energy_vec_R1(int)= Pore_barrier_energy(lambda,J_0m,tilt_0,1,chi,kappa_m,pore_line_tension);
    Pore_barrier_energy_vec_R2(int)= Pore_barrier_energy(lambda,J_0m,tilt_0,2,chi,kappa_m,pore_line_tension);
    Pore_barrier_energy_vec_R3(int)= Pore_barrier_energy(lambda,J_0m,tilt_0,3,chi,kappa_m,pore_line_tension);

    chi_vec(int)=chi;
    
	chi=chi-0.01;
    int=int+1;
end

%%
R_out=1;
int=1;
while (R_out<4)
    Pore_barrier_energy_Chi0(int)= Pore_barrier_energy_no_second_min(lambda,0,tilt_0,R_out,0,kappa_m,pore_line_tension);
	Pore_barrier_energy_Chi05(int)= Pore_barrier_energy_no_second_min(lambda,0,tilt_0,R_out,-0.5,kappa_m,pore_line_tension);
    Saddle_splay_vec_chi0_at_middle(int) = Saddle_Splay(lambda,tilt_0,R_out,0.00001,0);
    Saddle_splay_vec_chi01_at_middle(int) = Saddle_Splay(lambda,tilt_0,R_out,-0.1,0);
    Saddle_splay_vec_chi05_at_middle(int) = Saddle_Splay(lambda,tilt_0,R_out,-0.5,0);

    R_out_vec(int)=R_out;
    
	R_out=R_out+0.01;
    int=int+1;
end

%%
figure(1)
hold on
p1=plot(radius_vec,Delta_splay_vec_chi0);
p2=plot(radius_vec,Delta_splay_vec_chi01);
p3=plot(radius_vec,Delta_splay_vec_chi05);
xL=xlabel('radius[nm]');
yL=ylabel('Delta splay [nm^-^1]');
xL.FontSize=14;xL.FontWeight='bold';yL.FontSize=14;yL.FontWeight='bold';set(gcf,'color','white')
p1.LineWidth=2; p1.Color='blue';
p2.LineWidth=2; p2.Color='red';
p3.LineWidth=2; p3.Color='green';
legend('chi=0','chi=-0.1','chi=-0.5')

figure(2)
hold on
p1=plot(radius_vec,Saddle_splay_vec_chi0);
p2=plot(radius_vec,Saddle_splay_vec_chi01);
p3=plot(radius_vec,Saddle_splay_vec_chi05);
xL=xlabel('radius[nm]');
yL=ylabel('saddle splay [nm^-^2]');
xL.FontSize=14;xL.FontWeight='bold';yL.FontSize=14;yL.FontWeight='bold';set(gcf,'color','white')
p1.LineWidth=2; p1.Color='blue';
p2.LineWidth=2; p2.Color='red';
p3.LineWidth=2; p3.Color='green';
legend('chi=0','chi=-0.1','chi=-0.5')

figure(3)
hold on
p1=plot(radius_vec,Tilt_vec_chi0);
p2=plot(radius_vec,Tilt_vec_chi01);
p3=plot(radius_vec,Tilt_vec_chi05);
xL=xlabel('radius[nm]');
yL=ylabel('Tilt');
xL.FontSize=14;xL.FontWeight='bold';yL.FontSize=14;yL.FontWeight='bold';set(gcf,'color','white')
p1.LineWidth=2; p1.Color='blue';
p2.LineWidth=2; p2.Color='red';
p3.LineWidth=2; p3.Color='green';
legend('chi=0','chi=-0.1','chi=-0.5')

figure(4)
hold on
p1=plot(radius_vec,energy_density_chi0);
p2=plot(radius_vec,energy_density_chi01);
p3=plot(radius_vec,energy_density_chi05);
xL=xlabel('radius[nm]');
yL=ylabel('Energy Density [k_bT/nm^2]');
xL.FontSize=14;xL.FontWeight='bold';yL.FontSize=14;yL.FontWeight='bold';set(gcf,'color','white')
p1.LineWidth=2; p1.Color='blue';
p2.LineWidth=2; p2.Color='red';
p3.LineWidth=2; p3.Color='green';

legend('chi=0','chi=-0.1','chi=-0.5','chi=-0.5 no second')

figure(5)
hold on
p1=plot(chi_vec,Diph_energy_vec_R1);
p2=plot(chi_vec,Diph_energy_vec_R2);
p3=plot(chi_vec,Diph_energy_vec_R3);
xL=xlabel('chi');
yL=ylabel('Diphragm energy [k_bT]');
xL.FontSize=14;xL.FontWeight='bold';yL.FontSize=14;yL.FontWeight='bold';set(gcf,'color','white')
p1.LineWidth=2; p1.Color='blue';
p2.LineWidth=2; p2.Color='red';
p3.LineWidth=2; p3.Color='green';
legend('R=1','R=2','R=3')
set(gca, 'XDir','reverse')

figure(6)
hold on
p1=plot(chi_vec,Diph_energy_vec);
p2=plot(chi_vec,Diph_energy_vec_no_second_min);
xL=xlabel('chi');
yL=ylabel('Diphragm energy [k_bT]');
xL.FontSize=14;xL.FontWeight='bold';yL.FontSize=14;yL.FontWeight='bold';set(gcf,'color','white')
p1.LineWidth=2; p1.Color='blue';
p2.LineWidth=2; p2.Color='red';
p3.LineWidth=2; p3.Color='green';
legend('Second min','No second')
set(gca, 'XDir','reverse')


figure(7)
hold on
p1=plot(chi_vec,Pore_barrier_energy_vec_R1);
p2=plot(chi_vec,Pore_barrier_energy_vec_R2);
p3=plot(chi_vec,Pore_barrier_energy_vec_R3);
xL=xlabel('chi');
yL=ylabel('Pore energy [k_bT]');
xL.FontSize=14;xL.FontWeight='bold';yL.FontSize=14;yL.FontWeight='bold';set(gcf,'color','white')
p1.LineWidth=2; p1.Color='blue';
p2.LineWidth=2; p2.Color='red';
p3.LineWidth=2; p3.Color='green';
legend('R=1','R=2','R=3')
set(gca, 'XDir','reverse')

figure(8)
hold on
p1=plot(radius_vec_3,Barrier_energy_integral_chi0_R3);
p2=plot(radius_vec_3,Barrier_energy_integral_chi05_R3);
p3=plot(radius_vec_3,Barrier_energy_integral_chi1_R3);
xL=xlabel('radius[nm]');
yL=ylabel('Energy [k_bT]');
xL.FontSize=14;xL.FontWeight='bold';yL.FontSize=14;yL.FontWeight='bold';set(gcf,'color','white')
p1.LineWidth=2; p1.Color='blue';
p2.LineWidth=2; p2.Color='red';
p3.LineWidth=2; p3.Color='green';
legend('chi=0','chi=-0.5','chi=-1')



figure(9)
hold on
p1=plot(radius_vec_1,Barrier_energy_integral_chi05_R1);
p2=plot(radius_vec_2,Barrier_energy_integral_chi05_R2);
p3=plot(radius_vec_3,Barrier_energy_integral_chi05_R3);
xL=xlabel('radius[nm]');
yL=ylabel('Energy [k_bT]');
xL.FontSize=14;xL.FontWeight='bold';yL.FontSize=14;yL.FontWeight='bold';set(gcf,'color','white')
p1.LineWidth=2; p1.Color='blue';
p2.LineWidth=2; p2.Color='red';
p3.LineWidth=2; p3.Color='green';
legend('R=1','R=2','R=3')

figure(10)
hold on
p1=plot(R_out_vec,max(Pore_barrier_energy_Chi0)-Pore_barrier_energy_Chi0);
p2=plot(R_out_vec,max(Pore_barrier_energy_Chi05)-Pore_barrier_energy_Chi05);
xL=xlabel('diaphragm Radius[nm]');
yL=ylabel('delta Pore Energy Barrier[k_bT]');
xL.FontSize=14;xL.FontWeight='bold';yL.FontSize=14;yL.FontWeight='bold';set(gcf,'color','white')
p1.LineWidth=2; p1.Color='blue';
p2.LineWidth=2; p2.Color='red';
p3.LineWidth=2; p3.Color='green';
legend('chi=0','chi=-0.5')

figure(11)
hold on
p1=plot(R_out_vec,Saddle_splay_vec_chi0_at_middle);
p2=plot(R_out_vec,Saddle_splay_vec_chi01_at_middle);
p3=plot(R_out_vec,Saddle_splay_vec_chi05_at_middle);

xL=xlabel('diaphragm Radius[nm]');
yL=ylabel('Saddle splay [nm^-^2]');
xL.FontSize=14;xL.FontWeight='bold';yL.FontSize=14;yL.FontWeight='bold';set(gcf,'color','white')
p1.LineWidth=2; p1.Color='blue';
p2.LineWidth=2; p2.Color='red';
p3.LineWidth=2; p3.Color='green';
legend('chi=0','chi=-0.1','chi=-0.5')


