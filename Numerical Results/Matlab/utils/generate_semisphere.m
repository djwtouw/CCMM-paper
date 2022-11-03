%% Generate Semi-spherical Shell
function [X,label] = generate_semisphere(r1_inner,r2_inner,r1_outer,r2_outer,n)
%% Using polar coordinate
radius_inner = (r2_inner-r1_inner)*rand(1,n) + r1_inner;
radius_outer = (r2_outer-r1_outer)*rand(1,n) + r1_outer;
theta_inner = 2*pi*rand(1,n);
psi_inner = 0.5*pi*rand(1,n);
theta_outer = 2*pi*rand(1,n);
psi_outer = 0.5*pi*rand(1,n);
x_inner = radius_inner.*(sin(psi_inner).*cos(theta_inner));
y_inner = radius_inner.*(sin(psi_inner).*sin(theta_inner));
z_inner = radius_inner.*cos(psi_inner);
x_outer = radius_outer.*(sin(psi_outer).*cos(theta_outer));
y_outer = radius_outer.*(sin(psi_outer).*sin(theta_outer));
z_outer = radius_outer.*cos(psi_outer);
figure;
scatter3(x_inner,y_inner,z_inner,'r.');
hold on;
scatter3(x_outer,y_outer,z_outer,'b.');
hold off;
A_inner = [x_inner;y_inner;z_inner];
A_outer = [x_outer;y_outer;z_outer];
X = [A_inner,A_outer];
label = [ones(n,1);2*ones(n,1)];
end