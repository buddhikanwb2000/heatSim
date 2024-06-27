%Tt=A*Txx+Tyy+Bq
clear;clf;
          % Length of the domain in y-direction
T = 1;           % Total simulation time
        % Number of time steps (adjust as needed)
rho=10490;
cp=235;
K = 429;

%domain
Lx = 200;          % Length of the domain in x-direction
Ly = 1;

Nx=1001;Ny=51;
Nt=500;
dx=Lx/(Nx-1);
dy=Ly/(Ny-1);

%CFL condition
c=1;
C=0.1;
dt=C*dx/c;

%Field Variables

Tn=zeros(Ny,Nx);
x=linspace(0,Lx,Nx);
y=linspace(0,Ly,Ny);
[X,Y]=meshgrid(x,y);
%The thermal conductivity
K = ones(Ny,Nx)*429;%Thermal conductivity
%initial condition
Tn(:,:)=300;t=0;

%loop
for n = 1:Nt
    Tc=Tn;
    t=t+dt;
    for i=2:Nx-1
        for j=2:Ny-1
        Tn(j,i)=Tc(j,i)+dt*(K(j,i)/rho/cp)*((Tc(j,i+1)+Tc(j+1,i)-4*Tc(j,i)+Tc(j,i-1)+Tc(j-1,i)))/dx/dx;
        end
    end
%Source term
    Sx=round(7*Nx/Lx);
    Sy=round(3*Ny/Ly);
    Sy=1;
    if (t<0.5)
        Tn(Sy,Sx)=Tn(Sy,Sx)+dt*100/rho/cp;
    end
    %Boundary conditions
    Tn(1,:)=Tn(2,:);
    Tn(end,:)=Tn(end-1,:);
    Tn(:,1)=Tn(:,2);
    Tn(:,end)=Tn(:,end-1);
    %Visualize
    mesh(x,y,Tn);axis([0 Lx 0 Ly 280 350]);
    xlabel("Distance Along the rod");ylabel("Temperature");
    title("Time =%f seconds",t);
    pause(0.01);
%     figure(2);
%     imagesc(x, y, Tn);
%     colorbar;
%     title('Temperature Field Heatmap');
%     xlabel('X-axis');
%     ylabel('Y-axis');
%     axis xy;  % To ensure the correct orientation of the y-axis

% Optionally, add grid lines for better visualization
grid on;
end
return;