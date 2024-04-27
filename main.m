clear;
clc;

savePlots = true;

const.m = 1;
const.I = const.m/3;
const.cMac = 0.1;
const.xac = 0.25;
const.aw = 6;
const.ckphi = 2;
const.cky = 10;
const.cmuy = 0.02;
const.cmuphi = 0.02;

%% Question 1

a1 = const.cMac/const.I;
b1 = (const.xac*const.aw) / const.I;
c1 = -const.ckphi/const.I;
d1 = -const.cmuphi/const.I;

b2 = const.aw / const.m;
c2 = -const.cky/const.m;
d2 = -const.cmuy/const.m;

B = [0; 0; 0; a1];

A = [0, 1, 0, 0;
    c2, d2-b2, b2, 0;
    0, 0, 0, 1;
    0, -b1, b1+c1, d1];

%% Question 2

eigenvaluesA = eig(A);
absEigenvaluesA = sort(abs(eigenvaluesA), 'descend');

fprintf('Eigenvalues of A:\n');

for j = 1:4
    fprintf('%f + i*(%f)\n',real(eigenvaluesA(j)),imag(eigenvaluesA(j)));
end

fprintf('\n\n');

dt0 = 1/absEigenvaluesA(1);
dt1 = dt0 * (1/2);
dt = dt0 * (1/2)^2;

fprintf('Initial time step: %f\n\n', dt0);
fprintf('Stable time step: %f\n\n', dt);

figures.f1 = figure(1);
plot(real(eigenvaluesA).*dt0, imag(eigenvaluesA).*dt0, 'x', 'Color', "blue")
hold on
plot(real(eigenvaluesA).*dt1, imag(eigenvaluesA).*dt1, 'x', 'Color', "green")
plot(real(eigenvaluesA).*dt, imag(eigenvaluesA).*dt, 'x', 'Color', "red")
drawCircle(-1, 0, 1);
yline(0)
xline(0)
xlabel("Re", "Interpreter", "latex", FontSize = 15)
ylabel("Im", "Interpreter", "latex",'rotation',0, FontSize = 15)
title("Stability Region", "Interpreter", "LaTeX", FontSize = 20)
legend(sprintf("$\\Delta t^*= %f$", dt0), sprintf("$\\Delta t^*= %f$", dt1), sprintf("$\\Delta t^*= %f$", dt), "Interpreter", "latex")
xlim([-1.5, 0.2])
ylim([-1.1 1.1])
hold off

if savePlots
    saveas(gcf, "stability.png");
end

%% Question 3

x0 = [0; 0; 0; 0];

tExpEuler = [0, dt];
xExpEuler = [x0, eulerExplict_oneStep(@(tExpEuler,xExpEuler) diffEq(tExpEuler,xExpEuler,A,B), tExpEuler(end), x0, dt)];

condition = 1e-8;

while (abs(xExpEuler(2, end)) > condition) | (abs(xExpEuler(4, end)) > condition)

    xExpEuler = [xExpEuler, eulerExplict_oneStep(@(tExpEuler,xExpEuler) diffEq(tExpEuler,xExpEuler,A,B), tExpEuler(end), xExpEuler(:, end), dt)];
    tExpEuler = [tExpEuler, tExpEuler(end) + dt];

end

figures.f2 = figure(2);
plot(tExpEuler, xExpEuler(1,:))
hold on
plot(tExpEuler, xExpEuler(2,:))
plot(tExpEuler, xExpEuler(3,:))
plot(tExpEuler, xExpEuler(4,:))
xlabel("Time (s)", "Interpreter", "latex", FontSize = 15)
title("y, $\dot y$,$\phi$ and $\dot \phi$ V.S. time (Euler Explicit Method)", "Interpreter", "LaTeX", FontSize = 15)
legend("y","$\dot y$","$\phi$","$\dot \phi$", "Interpreter", "LaTeX", FontSize = 15)

if savePlots
    saveas(gcf, "eulerExplicit.png");
end

%% Question 4

tTrap = [0, dt];
xTrap = [x0, trapezoidal_oneStep(x0, dt, A, B)];

condition = 1e-8;

while (abs(xTrap(2, end)) > condition) | (abs(xTrap(4, end)) > condition)

    xTrap = [xTrap, trapezoidal_oneStep(xTrap(:,end), dt, A, B)];
    tTrap = [tTrap, tTrap(end) + dt];

end

figures.f3 = figure(3);
plot(tTrap, xTrap(1,:))
hold on
plot(tTrap, xTrap(2,:))
plot(tTrap, xTrap(3,:))
plot(tTrap, xTrap(4,:))
xlabel("Time (s)", "Interpreter", "latex", FontSize = 15)
title("y, $\dot y$,$\phi$ and $\dot \phi$ V.S. time (Trapezoidal Method)", "Interpreter", "LaTeX", FontSize = 15)
legend("y","$\dot y$","$\phi$","$\dot \phi$", "Interpreter", "LaTeX", FontSize = 15)

if savePlots
    saveas(gcf, "trapezoidal.png");
end

%% Question 5

fprintf("Computation of question 5 might take a while... (takes around 5min)\n\n")

tend = 7;

case1.dt = 1e-2;
case2.dt = 5e-3;
case3.dt = 1e-3;
case4.dt = 5e-4;
case5.dt = 1e-4;

cases = [case1, case2, case3, case4, case5];

for j = 1:5

    

    cases(j).tExpEuler = [0, cases(j).dt];
    cases(j).xExpEuler = [x0, eulerExplict_oneStep(@(t,x) diffEq(t,x,A,B), cases(j).tExpEuler(end), x0, cases(j).dt)];

    cases(j).tTrap = [0, cases(j).dt];
    cases(j).xTrap = [x0, trapezoidal_oneStep(x0, cases(j).dt, A, B)];

    % Euler Explicit

    while cases(j).tExpEuler < tend

        cases(j).xExpEuler = [cases(j).xExpEuler, eulerExplict_oneStep(@(t,x) diffEq(t,x,A,B), cases(j).tExpEuler(end), cases(j).xExpEuler(:, end), cases(j).dt)];
        cases(j).tExpEuler = [cases(j).tExpEuler, cases(j).tExpEuler(end) + cases(j).dt];

    end

    cases(j).index2Euler = find(cases(j).tExpEuler >= 2, 1);
    cases(j).index4Euler = find(cases(j).tExpEuler >= 4, 1);
    cases(j).index6Euler = find(cases(j).tExpEuler >= 6, 1);

    % Trapezoidal

    while cases(j).tTrap < tend

        cases(j).xTrap = [cases(j).xTrap, trapezoidal_oneStep(cases(j).xTrap(:,end), cases(j).dt, A, B)];
        cases(j).tTrap = [cases(j).tTrap, cases(j).tTrap(end) + cases(j).dt];

    end

    cases(j).index2Trap = find(cases(j).tTrap >= 2, 1);
    cases(j).index4Trap = find(cases(j).tTrap >= 4, 1);
    cases(j).index6Trap = find(cases(j).tTrap >= 6, 1);

    switch cases(j).dt
        case 1e-2
            fprintf("Case 1 computed.\n")
        case 5e-3
            fprintf("Case 2 computed.\n")
        case 1e-3
            fprintf("Case 3 computed. Be patient. \n")
        case 5e-4
            fprintf("Case 4 computed. No, the program didn't crash.\n")
        case 1e-4
            fprintf("Case 5 computed. Finally !!\n")
    end
    
end

cases(5).errorEuler = 0;
cases(5).errorTrap = 0;

for j = 1:4

    cases(j).error2Euler = norm(cases(j).xExpEuler(:, cases(j).index2Euler) - cases(5).xExpEuler(:, cases(5).index2Euler));
    cases(j).error4Euler = norm(cases(j).xExpEuler(:, cases(j).index4Euler) - cases(5).xExpEuler(:, cases(5).index4Euler));
    cases(j).error6Euler = norm(cases(j).xExpEuler(:, cases(j).index6Euler) - cases(5).xExpEuler(:, cases(5).index6Euler));

    cases(j).error2Trap = norm(cases(j).xTrap(:, cases(j).index2Trap) - cases(5).xTrap(:, cases(5).index2Trap));
    cases(j).error4Trap = norm(cases(j).xTrap(:, cases(j).index4Trap) - cases(5).xTrap(:, cases(5).index4Trap));
    cases(j).error6Trap = norm(cases(j).xTrap(:, cases(j).index6Trap) - cases(5).xTrap(:, cases(5).index6Trap));

end
%%
%load("question5.mat")

dtCases = [cases(1).dt, cases(2).dt, cases(3).dt, cases(4).dt, cases(5).dt];
t2errorEuler = [cases(1).error2Euler, cases(2).error2Euler, cases(3).error2Euler, cases(4).error2Euler, cases(5).errorEuler];
t4errorEuler = [cases(1).error4Euler, cases(2).error4Euler, cases(3).error4Euler, cases(4).error4Euler, cases(5).errorEuler];
t6errorEuler = [cases(1).error6Euler, cases(2).error6Euler, cases(3).error6Euler, cases(4).error6Euler, cases(5).errorEuler];

t2errorTrap = [cases(1).error2Trap, cases(2).error2Trap, cases(3).error2Trap, cases(4).error2Trap, cases(5).errorTrap];
t4errorTrap = [cases(1).error4Trap, cases(2).error4Trap, cases(3).error4Trap, cases(4).error4Trap, cases(5).errorTrap];
t6errorTrap = [cases(1).error6Trap, cases(2).error6Trap, cases(3).error6Trap, cases(4).error6Trap, cases(5).errorTrap];

figures.f4 = figure(4);
tiledlayout(2,1)

nexttile;
loglog(dtCases, t2errorEuler, 'x-')
hold on
loglog(dtCases, t4errorEuler, 'x-')
loglog(dtCases, t6errorEuler, 'x-')
xlabel("$\Delta t$", "Interpreter", "latex", FontSize = 15)
ylabel("Error", "Interpreter", "latex", FontSize = 15)
title("Error of Euler Explicit Method", "Interpreter", "LaTeX", FontSize = 15)
legend("t = 2s", "t = 4s", "t = 6s", "Interpreter", "latex", Location="southeast")
grid on

nexttile;
loglog(dtCases, t2errorTrap, 'x-')
hold on
loglog(dtCases, t4errorTrap, 'x-')
loglog(dtCases, t6errorTrap, 'x-')
xlabel("$\Delta t$", "Interpreter", "latex", FontSize = 15)
ylabel("Error", "Interpreter", "latex", FontSize = 15)
title("Error of Trapezoidal Method", "Interpreter", "LaTeX", FontSize = 15)
legend("t = 2s", "t = 4s", "t = 6s", "Interpreter", "latex", Location="southeast")
grid on

if savePlots
    saveas(gcf, "error.png");
end


function xdot = diffEq (t, x, A, B)    
%   Differential equation to be sovled
%
%   Args:
%       t: time
%       x: state vector
%       A: state matrix
%       B: independent values vector
%
%   Returns:
%       xdot: derivative of the state vector

    xdot = A*x + B;

end



function p = drawCircle(x,y,r)
%   Draws a circle in the plot given the position of the centre and the radius.
%
%   Args:
%       x: x position of the centre
%       y: y position of the centre
%       r: radius of the circle
%
%   Returns:
%       p: plot object

    theta = 0:pi/50:2*pi;

    xPos = r * cos(theta) + x;
    yPos = r * sin(theta) + y;

    p = plot(xPos, yPos);

end