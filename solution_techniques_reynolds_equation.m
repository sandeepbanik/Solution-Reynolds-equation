
%% Film thickness lubrication
clear all;
clc;
%% Long bearing approximation
X = 10;                                     % Length of bearing
Z = 100;                                    % Breadth of bearing
U = 20;                                     % Velocity (m/s)
hmax = 0.04;                                % Maximum wedge thickness
hmin = 0.02;                                % Minimum wedge thickness
C = (hmax + hmin)/2;                        % average wedge thickness
eta = 10e-03;
RHS = -0.03*2/(10*3);
factor = X^2/Z^2;
fun_hbar = @(x_bar,X)((hmax - (hmax - hmin)*x_bar)/C);
fun_hhalf = @(h)(-((h(2)-h(1))/2) + h(2));  % fluid film half node (i-0.5)
fun_hhalf2 = @(h)(((h(2)-h(1))/2) + h(1));  % fluid film half node (i+0.5)
del_x = 0.1;
x_bar = 0:del_x:1;
h = fun_hbar(x_bar,X);
h_halfminus = zeros(size(x_bar,2),1);
for i=2:length(h)-1
   h_halfminus(i) = fun_hhalf([h(i-1) h(i)]);
end
h_halfplus = zeros(size(x_bar,2),1);
for i=2:length(h)-1
   h_halfplus(i) = fun_hhalf2([h(i) h(i+1)]);
end

h_plus = h_halfplus.^3;
h_minus = h_halfminus.^3;
k = 20;
p = zeros(length(x_bar),1);
error = 100;
init_p = p;
p_bar = p;
while error > 0.01
    for i=2:size(p,1)-1
        p(i) = (h_plus(i)*p(i+1) + h_minus(i)*p(i-1) -RHS*del_x^2)/(h_plus(i) + h_minus(i));
    end
    error = (norm(p,1)-norm(init_p,1))/norm(p,1);
    init_p = p;
    p_bar = [p_bar p];
end
% Rescaling pressure 
P = p_bar(:,end)*6*eta*U*(X/1000)^2/(C/1000)^3;

% Comparison of p value
p_bar1 = zeros(length(x_bar),1);
for i=1:length(x_bar)
    p_bar1(i) = 0.0067*((x_bar(i)*(1-x_bar(i)))/(2-x_bar(i))^2);
end

figure;
subplot(1,2,1);
plot(x_bar*X,P);
xlabel('X(mm)');
ylabel('Pressure');
title('Pressure numerical');
subplot(1,2,2);
plot(x_bar*X,p_bar1);
xlabel('X(mm)');
ylabel('Pressure');
title('Pressure long bearing approximation');
pause;
close all;

%% Short bearing approximation

fun_p2d = @(z,x)((1e-06/(2/3*(2-x))^3)*z*(1-z));
z = 0:0.1:1;
p_compare = zeros(length(z),length(x_bar));
for i=1:length(z)
    for j=1:length(x_bar)
        p_compare(j,i) =  fun_p2d(z(i),x_bar(j));
    end
end
colormap(jet);
surf(x_bar*X,z,p_compare);
xlabel('Y(mm)');
ylabel('X(mm)');
ylabel('Pressure');
title('Pressure short bearing approximation');
pause;
close all;

%% Hybrid approach
% using the specifications as above
clear X
clear Z
X = 100;
Z = 10;
del_z = 0.1;
del_x = 0.1;
x_hyd = 0:del_x:1;
z_hyd = 0:del_z:1;
p_app1 = @(x)(0.0009/4*((3/(2-x))-(2/(2-x)^2)-1));
p_app2 = @(x,z)(2.7/80000*((z-z^2)/(2-x)^3));
p_hyd = zeros(length(x_hyd),length(z_hyd));
for i=2:size(p_hyd,1)-1
    for j=2:size(p_hyd,2)-1
        p_hyd(i,j) = (p_app1(x_hyd(i))*p_app2(x_hyd(i),z_hyd(j)))/(p_app1(x_hyd(i)) + p_app2(x_hyd(i),z_hyd(j)));
    end
end
figure;
colormap(jet);
surf(z_hyd*Z,x_hyd*X,p_hyd);
xlabel('Z-direction(mm)');
ylabel('X-direction(mm)');
zlabel('Non-dimensional pressure');

clear X
clear Z
X = 10;
Z = 100;
p_app1 = @(x)(0.009/4*((3/(2-x))-(2/(2-x)^2)-1));
p_app2 = @(x,z)(2.7/8*((z-z^2)/(2-x)^3));
p_hyd = zeros(length(x_hyd),length(z_hyd));
for i=2:size(p_hyd,1)-1
    for j=2:size(p_hyd,2)-1
        p_hyd(i,j) = (p_app1(x_hyd(i))*p_app2(x_hyd(i),z_hyd(j)))/(p_app1(x_hyd(i)) + p_app2(x_hyd(i),z_hyd(j)));
    end
end
figure;
colormap(jet);
surf(z_hyd*Z,x_hyd*X,p_hyd);
xlabel('Z-direction(mm)');
ylabel('X-direction(mm)');
zlabel('Non-dimensional pressure');
pause;
close all;

%% Finite difference method
clear X;
clear Z;
clear C;
X = 10;
Z = 25;
n = 50;                                 % No of nodes on x direction
m = 50;                                 % No of nodes on z direction
del_x = 1/n;
del_z = 1/m;
% eta_fd = 0.01;                        % Different eta
% eta_fd = 18.6e-03;                    % Low viscous oil
eta_fd = 0.1586;                        % High viscous oil
x_bar = 0:del_x:1;
z_bar = 0:del_z:1;
U_fd = 20;
hmax = 0.04;                            % Maximum wedge thickness
hmin = 0.02;                            % Minimum wedge thickness
C = (hmax + hmin)/2;                    % average wedge thickness
h_bar_fd = fun_hbar(x_bar,X);           % non dimensionalized film thickness
h_m_fd = zeros(size(x_bar,2),1);        % non dimensionalized half node i-0.5
for i=2:length(h_bar_fd)-1
   h_m_fd(i) = fun_hhalf([h_bar_fd(i-1) h_bar_fd(i)]);
end
h_p_fd = zeros(size(x_bar,2),1);        % non dimensionalized half node i + 0.5
for i=2:length(h_bar_fd)-1
   h_p_fd(i) = fun_hhalf2([h_bar_fd(i) h_bar_fd(i+1)]);
end
h_fd = h_bar_fd.^3;
h_plus_fd = h_p_fd.^3;
h_minus_fd = h_m_fd.^3;

X_Z_const = X^2/Z^2;

iter = 1000;
p_bar_fd = zeros(length(x_bar),length(z_bar));
A_fun = @(h_p,h_m,h)(h_p/(h_p + h_m + 2*X_Z_const*h));
B_fun = @(h_p,h_m,h)(h_m/(h_p + h_m + 2*X_Z_const*h));
C_fun = @(h_p,h_m,h)((X_Z_const*h)/(h_p + h_m + 2*X_Z_const*h));
E_fun = @(h_p,h_m,h_b1,h_b2,h)((-del_x*C/(2*X))*(h_b2-h_b1)/(h_p + h_m + 2*X_Z_const*h));
p_init = p_bar_fd;
error_target = 0.0001;
h = waitbar(0,'Please wait...');
for k=1:iter
    for j=2:size(p_bar_fd,2)-1
        for i=2:size(p_bar_fd,1)-1
            p_bar_fd(i,j) = A_fun(h_plus_fd(i),h_minus_fd(i),h_fd(i))*p_bar_fd(i+1,j) +...
                B_fun(h_plus_fd(i),h_minus_fd(i),h_fd(i))*p_bar_fd(i-1,j) +...
                C_fun(h_plus_fd(i),h_minus_fd(i),h_fd(i))*p_bar_fd(i,j+1) +...
                C_fun(h_plus_fd(i),h_minus_fd(i),h_fd(i))*p_bar_fd(i,j-1) +...
                E_fun(h_plus_fd(i),h_minus_fd(i),h_bar_fd(i-1),h_bar_fd(i+1),h_fd(i));
        end
    end
    error = (sum(p_bar_fd(:)) - sum(p_init(:)))/sum(p_bar_fd(:));
    p_init = p_bar_fd;
    if error < error_target
        break;
    end
    waitbar(k/iter);
end
close(h);
P_scale = 1000*p_bar_fd*6*eta_fd*U_fd*X^2/(C^3);
x_scale = x_bar.*X;
z_scale = z_bar.*Z;
colormap(jet);
surf(z_scale,x_scale,P_scale);
title('Pressure profile');
xlabel('Y(mm)');
ylabel('X(mm)');
zlabel('Pressure');
pause;
close all;

%% Hydrodynamic lubrication with viscosity variation
al = 0.020730e-06;                          % Pressure coefficient
% al = 20e-09;
eta_fun = @(p)(eta_fd*exp(al*p));
for i=1:size(P_scale,1)
    for j=1:size(P_scale,2)
        p_vis(i,j) = (1/al)*log(1/(1-P_scale(i,j)*al));
    end
end
colormap(jet);
surf(z_scale,x_scale,p_vis);
title('Pressure viscosity profile');
xlabel('Y(m)');
ylabel('X(m)');
zlabel('Pressure');
pause;
close all;

%% Estimating elastic deformation
clear all;

X = 10;
Z = 25;
n = 50;                                                         % No of nodes on x direction
m = 50;                                                         % No of nodes on z direction
del_x = 1/n;
del_z = 1/m;
eta_fd = 0.1586;
x_bar = 0:del_x:1;
z_bar = 0:del_z:1;
U_fd = 20;
hmax = 0.04;                                                    % Maximum wedge thickness
hmin = 0.02;                                                    % Minimum wedge thickness
C = (hmax + hmin)/2;                                            % average wedge thickness
fun_hbar = @(x_bar,X)((hmax - (hmax - hmin)*x_bar)/C);

h_bar_x = fun_hbar(x_bar,X);                                   % non dimensionalized film thickness
h_bar = zeros(length(x_bar),length(z_bar));

for i=1:size(h_bar,2)
    for j=1:size(h_bar,1)
        h_bar(:,i)  = h_bar_x(i);
    end
end
h_bar = h_bar';
X_Z_const = X^2/Z^2;
nu = 0.305;                                                     % Poisson's ratio
E_s = 180*10^9;                                                 % Elastic modulus
iter = 1000;
p_bar_fd = zeros(length(x_bar),length(z_bar));
A_fun = @(h,h_p,h_m)((0.75*(h_p - h_m)+h)/(2*h*(1+X_Z_const)));
B_fun = @(h,h_p,h_m)((h-0.75*(h_p - h_m))/(2*h*(1+X_Z_const)));
C_fun = @(h,h_p,h_m)(X_Z_const*(0.75*(h_p - h_m)+h)/(2*h*(1+X_Z_const)));
D_fun = @(h,h_p,h_m)(X_Z_const*(h-0.75*(h_p - h_m))/(2*h*(1+X_Z_const)));
E_fun = @(h,h_p,h_m)(C*del_x/(X*4)*(h_m-h_p)/(h^3*(1 + X_Z_const)));
p_init = p_bar_fd;
error_target = 0.0001;
delta_flag = 0;                                                 % Condition for elastic deformation
h = waitbar(0,'Please wait...');
for k=1:iter
    for i=2:size(p_bar_fd,2)-1
        for j=2:size(p_bar_fd,1)-1
            p_bar_fd(i,j) = A_fun(h_bar(i,j),h_bar(i+1,j),h_bar(i-1,j))*p_bar_fd(i+1,j) +...
                B_fun(h_bar(i,j),h_bar(i+1,j),h_bar(i-1,j))*p_bar_fd(i-1,j) +...
                C_fun(h_bar(i,j),h_bar(i,j+1),h_bar(i,j-1))*p_bar_fd(i,j+1) +...
                D_fun(h_bar(i,j),h_bar(i,j+1),h_bar(i,j-1))*p_bar_fd(i,j-1) +...
                E_fun(h_bar(i,j),h_bar(i+1,j),h_bar(i-1,j));
        end
    end
    if delta_flag == 1
        delta = 0;
    else
        delta = zeros(length(x_bar),length(z_bar));
        for m=2:size(delta,1)-1
            for l = 2:size(delta,2)-1
                temp = 0;
                del = 0;
                for i=2:size(p_bar_fd,2)-1
                    for j=2:size(p_bar_fd,1)-1
                        if i == l || j==m
                            temp = 0;
                        else
                            temp = (1-nu^2)/(pi*E_s*C)*(6*eta_fd*U_fd*X^3/C^3)*(p_bar_fd(i,j)*del_x/sqrt((i-l)^2*X_Z_const + (j-m)^2));
                        end
                        del = del + temp;
                    end
                end
                delta(m,l) = del;
            end
        end
    end
    h_bar = h_bar + delta;
    error = (sum(p_bar_fd(:)) - sum(p_init(:)))/sum(p_bar_fd(:));
    p_init = p_bar_fd;
    if error < error_target
        break;
    end
    waitbar(k / iter)
end
close(h);
pause;
P_scale = 1000*p_bar_fd*6*eta_fd*U_fd*X^2/(C^3);
x_scale = x_bar.*X;
z_scale = z_bar.*Z;
colormap(jet);
surf(z_scale,x_scale,P_scale);
title('Pressure profile');
xlabel('Y(m)');
ylabel('X(m)');
zlabel('Pressure');

% Only deformation after pressure 

if delta_flag == 1
    delta = zeros(length(x_bar),length(z_bar));
    for k=2:size(delta,1)-1
        for l = 2:size(delta,2)-1
            temp = 0;
            del = 0;
            for i=2:size(p_bar_fd,2)-1
                for j=2:size(p_bar_fd,1)-1
                    if i == l || j==k
                        temp = 0;
                    else
                        temp = (1-nu^2)/(pi*E_s*C)*(6*eta_fd*U_fd*X^3/C^3)*(p_bar_fd(i,j)*del_x/sqrt((i-l)^2*X_Z_const + (j-k)^2));
                    end
                    del = del + temp;
                end
            end
            delta(l,k) = del;
        end
    end
end
figure;
surf(z_scale,x_scale,delta);
title('Deformation profile');
xlabel('Y(mm)');
ylabel('X(mm)');
zlabel('Deformation');


al = 0.020730e-06;
% p_vis = 1/al*ln(1/(1 - P_scale.*al));
for i=1:size(P_scale,1)
    for j=1:size(P_scale,2)
        p_vis(i,j) = (1/al)*log(1/(1-P_scale(i,j)*al));
    end
end

colormap(jet);
figure;
surf(z_scale,x_scale,p_vis);
title('Pressure profile');
xlabel('Y(m)');
ylabel('X(m)');
zlabel('Pressure');
pause;
close all;

%% Hydrodynamic lubrication journal bearing 
clear all;

R = 100e-03;                                                    % Radial diameter of bearing
X = 2*pi;                                                       % Circumference of bearing
Z = 160e-03;                                                    % Length of the bearing
n = 50;                                                         % No of nodes on x direction
m = 50;                                                         % No of nodes on z direction
del_x = 1/n;
del_z = 1/m;
eta_fd = 0.1586;
x_bar = 0:del_x:1;
z_bar = 0:del_z:1;
rps = 1200/60;
U_fd = R*rps;                                                   % average wedge thickness
C = 0.18e-03;
e_p = 0.7;
e = e_p*C;
fun_hbar = @(theta)(e_p*cos(theta) + 1);
theta = x_bar*2*pi;
del_theta = theta(end)-theta(end-1);
h_bar_x = fun_hbar(theta)';                                   % non dimensionalized film thickness
dh = (1+e_p*cos(theta));


fun_hhalf = @(h)(-((h(2)-h(1))/2) + h(2));
fun_hhalf2 = @(h)(((h(2)-h(1))/2) + h(1));
h_halfminus = zeros(size(x_bar,2),1);
for i=2:length(h_bar_x)-1
   h_halfminus(i) = fun_hhalf([h_bar_x(i-1) h_bar_x(i)]);
end
h_halfplus = zeros(size(x_bar,2),1);
for i=2:length(h_bar_x)-1
   h_halfplus(i) = fun_hhalf2([h_bar_x(i) h_bar_x(i+1)]);
end

h_plus = h_halfplus.^3;
h_minus = h_halfminus.^3;

h_bar = h_bar_x;
% h_bar = h_bar';
X_Z_const = R^2/Z^2;
nu = 0.305;
E_s = 180*10^9;
iter = 5000;
p_bar_fd = zeros(length(x_bar),length(z_bar));
A_fun = @(h,h_p,h_m)(h_p*del_z^2/((h_p + h_m)*del_z^2 + 2*X_Z_const*h^3*del_theta^2));
B_fun = @(h,h_p,h_m)(h_m*del_z^2/((h_p + h_m)*del_z^2 + 2*X_Z_const*h^3*del_theta^2));
C_fun = @(h,h_p,h_m)(X_Z_const*h^3*del_theta^2/((h_p + h_m)*del_z^2 + 2*X_Z_const*h^3*del_theta^2));
E_fun = @(theta,h,h_p,h_m)(e_p*sin(theta)*del_z^2*del_theta^2/((h_p + h_m)*del_z^2 + 2*X_Z_const*h^3*del_theta^2));
p_init = p_bar_fd;
error_target = 0.0001;
h = waitbar(0,'Please wait...');
for k=1:iter
    for j=2:size(p_bar_fd,2)-1
        for i=2:size(p_bar_fd,1)-1
            p_bar_fd(i,j) = A_fun(h_bar(i),h_plus(i),h_minus(i))*p_bar_fd(i+1,j) +...
                B_fun(h_bar(i),h_plus(i),h_minus(i))*p_bar_fd(i-1,j) +...
                C_fun(h_bar(i),h_plus(i),h_minus(i))*p_bar_fd(i,j+1) +...
                C_fun(h_bar(i),h_plus(i),h_minus(i))*p_bar_fd(i,j-1) +...
                E_fun(theta(i),h_bar(i),h_plus(i),h_minus(i));
            if p_bar_fd(i,j) < -0.017
                p_bar_fd(i,j) =-0.017;
            end
        end
    end
    error = (sum(p_bar_fd(:)) - sum(p_init(:)))/sum(p_bar_fd(:));
    p_init = p_bar_fd;
    if error < error_target
        break;
    end
    waitbar(k / iter)
    
end
close(h);
x_scale = x_bar;
z_scale = z_bar.*Z;
colormap(jet);
surf(z_scale,theta,p_init);
title('Pressure profile');
xlabel('Y(m)');
ylabel('\theta');
zlabel('Pressure');




