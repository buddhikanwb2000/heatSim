p=0.5;
R=0.69;
r=100e-6;
a=1.144e-6;
gamma=1/a;
Lx=200e-6;
Ly=200e-6;
x=linspace(0,Lx,100);
y=linspace(0,Ly,10);
source=zeros(length(y),length(x));
disp(size(source));
for n=1:length(y)
    for q=1:length(x)
        source(n,q)=p*(1-R)*exp(-((x(q)-Lx/2)^2+(y(n)-Ly/2)^2)/r^2)*a*exp(-gamma*y(n))/pi/r^2;
    end
end
figure(1); % Creates a new figure window
surf(x, y, source); % Plots the source distribution in 3D
xlabel('x [m]'); % Label for the x-axis
ylabel('y [m]'); % Label for the y-axis
zlabel('Source Intensity'); % Label for the z-axis
title('3D Plot of the Source Distribution'); % Title for the plot
colorbar;