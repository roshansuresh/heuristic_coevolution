%% feasible designs fitness weight vs number of feasible designs
clear
close all
clc

%% Test plot
x = linspace(0,30,31);
y = 10*exp(-x);
figure
plot(x,y,'-*')
xlabel('Number of Feasible Designs')
ylabel('$\alpha$','Interpreter','Latex')
