clear variables; close all; clc
Nx = 10;
Ny = 101;
Lx = 1e-7;
Ly = 1e-6;
electron = 1.6e-19;
kB = 1.38e-23;
Temp = 273;
x = linspace(0,Lx-Lx/Nx,Nx);
y = linspace(0,Ly,Ny);
[X,Y] = meshgrid(x,y);
charge = 0.01*ones(size(X)).';

convertFactor = 9.64e4;

charge0 = -charge * convertFactor / 6.95E-10;
voltage = -5e-3;
charge_p = charge0.*exp(-electron*voltage/kB/Temp);
charge_n = charge0.*exp(electron*voltage/kB/Temp);

phi_old=fast_Poisson(Nx,Ny,Lx, Ly, charge_p*0,voltage);
omega = 0.05;

for i = 1:300
    phi=fast_Poisson(Nx,Ny,Lx, Ly, charge,voltage);
    phi = omega*phi + (1-omega)*phi_old;
    el2 = norm(phi-phi_old);
    phi_old = phi;

    
    charge_p = charge0.*exp(-electron*phi/kB/Temp);
    charge_n = charge0.*exp(electron*phi/kB/Temp);
    charge = charge_p-charge_n;
end




% contourf(X,Y,phi.');
% shading interp
% colormap('jet');
% colorbar;

figure
plot(y,phi(5,:))
figure
plot(y,charge_p(5,:)/(-convertFactor / 6.95E-10))
figure
plot(y,charge_n(5,:)/(-convertFactor / 6.95E-10))