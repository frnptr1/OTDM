clear;

% canonical form of the function, matlab read this form 
%f_rosen = @(x) 100*(x(2) - x(1).^2).^2 + (1 - x(1)).^2;
f_rosen = @(x) (1-x(1)).^2 + 100*(x(2)-x(1).^2).^2;

% first derivative/gradient method of the function
% g_rosen = @(x) [400*x(1).^3-400*x(1)*x(2)+2*x(1)-2 ; 200*x(2)-200*x(1).^2];
grad = @(x)[-400*(x(2) - x(1)^2)*x(1) - 2*(1 - x(1)), 200*(x(2) - x(1)^2)];

% second derivative/second order gradient/ hessian matrix

% h_rosen = @(x) [1200*x(1).^2-400*x(2)+2 -400*x(1);-400*x(1) 200];
hess = @(x)[1200*x(1)^2 - 400*x(2) + 2, -400*x(1);
            -400*x(1), 200];

% initialize the parameters
options=zeros(1,11);
% set the flag to 1 for 8th param
options( 2) = 15000;
options( 8) = 0;
options(10) = 1;
% options(20) = 1;

% starting point
x = [0.5;0.5];

[x, xk, alk, dk, gk, fk, iout] = otdm_uo_students(f_rosen, grad, hess, x, options);
xylim = [0,0,0,0];
otdm_uo_plot(f_rosen, xk, gk, xylim);

% value of the objective function
% euclidean norm of the gradient in the current point
% min eigenvalue of the Hessian matrix ( min lamba value)
% alpha at the current point, steplength

%algorithm stop because the objective function reached a desireable value
%under the tolerance


%h1 = hp12(x1) hessian matrix in point x1
%[L,D,P] = ldl(h1)  P is a permutation matrix
% ludenberger modification of the D matrix in a way that the minimum value
% of the diagonal matrix to be > delta = user parameter usually between
% 0.1-0.8, if this value increase we are asking for a deep modification of
% the hessian matrix
% mu= max(0, delta-min(diag(D))
%DL = D+mu*eye(2 - or n)
%H1L = P*L*DL*L'*P' the new updated hessian
%d = -inv(H1L)*gradient(x1)'
