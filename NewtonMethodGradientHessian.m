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
e = 10^(-8); 
i = 1; 

% Gradient and Hessian Computation:
df_dx = diff(f, X);
df_dy = diff(f, Y);
J = [subs(df_dx,[X,Y], [x(1),y(1)]) subs(df_dy, [X,Y], [x(1),y(1)])]; % Gradient
ddf_ddx = diff(df_dx,X);
ddf_ddy = diff(df_dy,Y);
ddf_dxdy = diff(df_dx,Y);
ddf_ddx_1 = subs(ddf_ddx, [X,Y], [x(1),y(1)]);
ddf_ddy_1 = subs(ddf_ddy, [X,Y], [x(1),y(1)]);
ddf_dxdy_1 = subs(ddf_dxdy, [X,Y], [x(1),y(1)]);
H = [ddf_ddx_1, ddf_dxdy_1; ddf_dxdy_1, ddf_ddy_1]; % Hessian
S = inv(H); 

% Optimization Loop:
while norm(J) > e
    I = [x(i),y(i)]';
    x(i+1) = I(1)-S(1,:)*J';
    y(i+1) = I(2)-S(2,:)*J';
    i = i+1;
    J = [subs(df_dx,[X,Y], [x(i),y(i)]) subs(df_dy, [X,Y], [x(i),y(i)])]; % Update Jacobian
    ddf_ddx_1 = subs(ddf_ddx, [X,Y], [x(i),y(i)]);
    ddf_ddy_1 = subs(ddf_ddy, [X,Y], [x(i),y(i)]);
    ddf_dxdy_1 = subs(ddf_dxdy, [X,Y], [x(i),y(i)]);
    H = [ddf_ddx_1, ddf_dxdy_1; ddf_dxdy_1, ddf_ddy_1]; % Update Hessian
    S = inv(H); % New Search Direction
end

% Result Table:`
Iter = 1:i;
X_coordinate = x';
Y_coordinate = y';
Iterations = Iter';
T = table(Iterations,X_coordinate,Y_coordinate);

% Plots:
fcontour(f,[0 500 0 900],'Fill', 'On');
hold on;
plot(x,y,'*-r');
grid on;

% Error Analysis:
Y_actual = (250*10^(-.75)*20^-.25)^2 
Y_error = (y(12)-Y_actual)/Y_actual;
X_actual = (250*20^(-.75)*10^-.25)^2
X_error = (x(12)-X_actual)/X_actual;

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

