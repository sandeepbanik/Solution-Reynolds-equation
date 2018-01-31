%% Contour maps for rectangular contact
% Sandeep Banik - 30/04/2017
% gh = 2.45*gv^0*ge^0 - isoviscous rigid regim
% gh = 1.05*gv^(2/3)*ge^0 - viscous rigid
% gh = 1.654*gv^(0.54)*ge^(0.06) - viscous elastic regim
% gh = 2.45*gv^0*ge^0.8 - isoviscous elastic regim
% This is a trail map and may be incorrect

clear all;
clc;

fun_gv_1 = @(gh)[(gh/1.05)^(3/2)];
fun_ge_1 = @(gh,gv)[(gh/(1.654*gv^(0.54)))^(1/0.06)];
fun_ge_2 = @(gh)[(gh/2.45)^(1/0.8)];
fun_gv_2 = @(gh,ge)[(gh/(1.654*ge^(0.06)))^(1/0.54)];

gh_xp = 20; % Film thickness (Unstable after 40)
n_xp = length(gh_xp);
x = zeros(n_xp,4);
y = zeros(n_xp,4);
x_line1 = zeros(n_xp,1);
x_line2 = zeros(n_xp,1);
y_line1 = zeros(n_xp,1);
y_line2 = zeros(n_xp,1);

for i=1:length(gh_xp)
    gv_high = fun_gv_1(gh_xp(i));
    ge_low = fun_ge_1(gh_xp(i),gv_high);
    ge_high = fun_ge_2(gh_xp(i));
    gv_low = fun_gv_2(gh_xp(i),ge_high);
   
    x(i,2) = ge_low;
    x_line1(i) = x(i,2);
    x(i,3) = ge_high;
    x_line2(i) = x(i,3);
    x(i,4) = ge_high;
    
    y(i,1) = gv_high;
    y(i,2) = gv_high;
    y_line1(i) = y(i,2);
    y(i,3) = gv_low;
    y_line2(i) = y(i,3);
    plot(x(i,:),y(i,:),'b');hold on;
end
plot(x_line1,y_line1,':','LineWidth',2,'Color','g');
plot(x_line2,y_line2,':','LineWidth',2,'Color','g');
axis([0.1 20 0.1 100]);
title('EHL regimes for rectangular contact');
xlabel('Elasticity parameter ge');
ylabel('Viscocity parameter gv');

ge_xp = ge_high:-0.1:ge_low;
for i=1:length(ge_xp)
    z(2*i) = fun_gv_2(gh_xp,ge_xp(i));
    z(2*i-1) = fun_ge_1(gh_xp,z(2*i));
end
