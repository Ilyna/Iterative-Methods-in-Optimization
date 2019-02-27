clc
clear all
close all
format long

% Function Definition:
syms X Y;
f = 1000*X^(1/4)*Y^(1/4)-20*X-10*Y;

% Initial Guess:
x(1) = 1;
y(1) = 1;
x(2) = 2.86;
y(2) = 2.9;
e = 10^(-8);
i = 2;

% Gradient and Hessian Computation:
df_dx = diff(f, X);
df_dy = diff(f, Y);
J = [subs(df_dx,[X,Y], [x(1),y(1)]) subs(df_dy, [X,Y], [x(1),y(1)])]; % Gradient
J2 = [subs(df_dx,[X,Y], [x(2),y(2)]) subs(df_dy, [X,Y], [x(2),y(2)])];

% Optimization Loop:
while norm(J) > e
    I = [x(i-1),y(i-1)]';
    I2 = [x(i),y(i)]';
    x(i+1) = I2(1)-J2(1)*(I2(1)-I(1))/(J2(1)-J(1));
    y(i+1) = I2(2)-J2(2)*(I2(2)-I(2))/(J2(2)-J(2));
    i = i+1;
    J = J2; %Update Jacobian
    J2 = [subs(df_dx,[X,Y], [x(i),y(i)]) subs(df_dy, [X,Y], [x(i),y(i)])];
end

% Result Table:`
Iter = 1:i;
X_coordinate = x';
Y_coordinate = y';
Iterations = Iter';
T = table(Iterations,X_coordinate,Y_coordinate);

% Plots:
fcontour(f, [0 500 0 900], 'Fill', 'On');
hold on;
plot(x,y,'*-r');
grid on;

% Error Analysis:
Y_actual = (250*10^(-.75)*20^-.25)^2; 
Y_error = (y(57)-Y_actual)/Y_actual;
X_actual = (250*20^(-.75)*10^-.25)^2;
X_error = (x(57)-X_actual)/X_actual;

% Output:
fprintf('Initial Function Value: %d\n\n',subs(f,[X,Y], [x(1),y(1)]));
if (norm(J) < e)
    fprintf('Maximum succesfully obtained...\n\n');
end
fprintf('Number of Iterations for Convergence: %d\n\n', i);
fprintf('Point of Maxima: [%d,%d]\n\n', x(i), y(i));
fprintf('Function Maximum Value after Optimization: %f\n\n', subs(f,[X,Y], [x(i),y(i)]));
disp(T)
fprintf('The percent error in X is: %d\n', X_error);
fprintf('The percent error in Y is: %d\n', Y_error);
