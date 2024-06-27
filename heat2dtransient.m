%MATLAB code to solve for the transient state conduction in explicit method.
clear all
close all
clc

%defining the boundary.
x = 1:10;
dx = abs(x(1)-x(2));
nx = length(x);
y = 1:10;
dy = abs(y(1)-y(2));
ny = length(y);

%defining the time steps
T = 1000;
dt = 0.1;

%defining the boundary conditions.
t_top = 300;%600;
t_bottom = 300;%900;
t_left = 300;%400;
t_right = 300;%800;

%defining the initial conditions.
t = 300*ones(nx,ny);
t(1,:) = t_top;
t(ny,:) = t_bottom;
t(:,1) = t_left;
t(:,nx) = t_right;
t(1,1) = (t_left + t_bottom)/2;
t(1,ny) = (t_top + t_left)/2;
t(nx,ny) = (t_top + t_right)/2;
t(nx,1) = (t_right + t_bottom)/2;

%initialising the t_initial and t_old for future computaions
t_initial = t;
t_old = t;

%defining the thermal diffusivity.
alpha = 1.2;

%initialising the k1 and k2 for future simplifications
k1 = alpha*(dt/(dx^2));
k2 = alpha*(dt/(dy^2));

%initialising the iteration count.
iteration = 0;

%starting the time loop. The loop is predetermined to run for 500 time steps.
for k = 1:500
    t(1 , :) = t(ny, :) + 10;
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
    t_old = t;
    iteration = iteration + 1;
    %boundary conditions
    %t(1,:)=0;
    t(end,:)=t(end-1,:);
    %t(:,1)=t(:,2);
    t(:,end)=t(:,end-1);
    %plotting the solution at each iteration.
    [c,h] = contourf(x,y,t);
    clabel(c,h);
    colorbar
    colormap(jet)
    set(gca,'YDIR','reverse');
    xlabel('X Axis');
    ylabel('Y Axis');
    title({['2D  Heat Conduction in Transient State.'],['Number of iterations = ',num2str(iteration)]});
    pause(0.01);
    
end
