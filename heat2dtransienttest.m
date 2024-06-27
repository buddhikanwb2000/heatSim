%MATLAB code to solve for the transient state conduction in explicit method.
clear all
close all
clc

%defining the boundary.
Lx=200e-6;
Ly=1e-6;
x=linspace(0,Lx,20);
y=linspace(0,Ly,10);
dx = abs(x(1)-x(2));
nx = length(x);
dy = abs(y(1)-y(2));
ny = length(y);

%convection
h_conv=10;
T_ambient=300;
%%
%defining the boundary conditions.
t_top = 300;%600;
t_bottom = 300;%900;
t_left = 300;%400;
t_right = 300;%800;

%defining the initial conditions.
t = 300*ones(nx,ny);
% t(1,:) = t_top;
% t(ny,:) = t_bottom;
% t(:,1) = t_left;
% t(:,nx) = t_right;
% t(1,1) = (t_left + t_bottom)/2;
% t(1,ny) = (t_top + t_left)/2;
% t(nx,ny) = (t_top + t_right)/2;
% t(nx,1) = (t_right + t_bottom)/2;

%initialising the t_initial and t_old for future computaions
t_initial = t;
t_old = t;

%defining the thermal diffusivity.
rho=10400;
cp=235;
K = 429;%bulk silver density,C_p and thermal conductivity
alpha = K/(rho*cp);    % Thermal diffusivity

%defining the time steps
T = 1000;
dt = min(dy,dx)^2/alpha/2;
%dt=0.1;

%% 
%source
p=0.5;
R=0.69;
r=100e-6;
a=1.144e-8;
gamma=1/a;

source=zeros(length(y),length(x));
disp(size(source));
for n=1:length(y)
    for q=1:length(x)
        source(n,q)=p*(1-R)*exp(-(x(q)-Lx/2)^2/r^2)*a*exp(-gamma*y(n))/pi/r^2;
    end
end
figure(1); % Creates a new figure window
surf(x, y, source); % Plots the source distribution in 3D
xlabel('x [m]'); % Label for the x-axis
ylabel('y [m]'); % Label for the y-axis
zlabel('Source Intensity'); % Label for the z-axis
title('3D Plot of the Source Distribution'); % Title for the plot
colorbar;

%% 
%initialising the k1 and k2 for future simplifications
k1 = alpha*(dt/(dx^2));
k2 = alpha*(dt/(dy^2));

%initialising the iteration count.
iteration = 0;
%% 
%starting the time loop. The loop is predetermined to run for 500 time steps.

for k = 1:1e5
    %t = t + dx^2*dt*source'/cp/rho;
    
    
    %starting the spatial loops.
    
    for j = 2:(ny-1)
        
        for i = 2:(nx-1)
            
            %calculating the temperature at a perticular point.
            term1 = (t_old(i-1,j) - 2*t_old(i,j) + t_old(i+1,j));
            term2 = (t_old(i,j-1) - 2*t_old(i,j) + t_old(i,j+1));
            t(i,j) = t_old(i,j) + (term1*k1) + (term2*k2);
            
        end
        
    end
    
    %updating the t_old and the iteration count,
    t = t + source'/cp/rho;
    t_old = t;
    iteration = iteration + 1;
    
    %boundary conditions
     t(1,:)=t( 2,:)+0.01 * (T_ambient - t(1, :)) * dt;
     t(end,:)=t(19, :) + 0.01* (T_ambient - t(end, :)) * dt;
      %t(:,1)=t(:, 2)+ 0.01 * (T_ambient - t(:, 1)) * dt;
     %t(:,end)=t(:, 9);%+ 30 * (T_ambient - t(:, end)) * dt;
    %plotting the solution at each iteration.
%     [c,h] = contourf(x,y,t);
%     clabel(c,h);
%     colorbar
%     colormap(jet)
%     set(gca,'YDIR','reverse');
%     xlabel('X Axis');
%     ylabel('Y Axis');
%     title({['2D  Heat Conduction in Transient State.'],['Number of iterations = ',num2str(iteration)]});
%     pause(0.01);
%     
%     surf(y, x,t);%axis([0 Lx 0 Ly 280 400]);
%     title(['Time = ' num2str(k*dt)]);
%     xlabel('X-axis');
%     ylabel('Y-axis');
%     zlabel('Temperature');
%     %axis([0 Lx 0 Ly 280 max(max(u))]);  % Adjust axis limits if needed
%     shading interp;  % Interpolate colors for smoother visualization
%     drawnow;
end
%%
figure(2)
surf(y, x,t);%axis([0 Lx 0 Ly 280 400]);
title(['Temperature of the cross section at Time = ' num2str(k*dt)]);
xlabel('X-axis'); 
ylabel('Y-axis');
zlabel('Temperature');
%axis([0 Lx 0 Ly 280 max(max(u))]);  % Adjust axis limits if needed
%shading interp;  % Interpolate colors for smoother visualization

