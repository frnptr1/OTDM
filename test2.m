clear;

% canonical form of the function, matlab read this form 
f_x = @(x) x(1)^2+(1/4)*x(2)^3+0.5*x(1)*x(2)+(15/4)*x(1)-(9/16)*x(2);
% first derivative/gradient method of the function
g_x = @(x) [2*x(1)+0.5*x(2)+(15/4), (3/4)*x(2)^2+0.5*x(1)-(9/16)];
% second derivative/second order gradient/ hessian matrix
h_x = @(x)[2, 0.5;
            0.5, (3/2)*x(2)];

% initialize the parameters
options=zeros(1,11);
% set the flag to 1 for 8th param
options( 8) = 0;
options(10) = 7;
options(20) = 1;

% starting point
x = [-3;-1];

[x, xk, alk, dk, gk, fk, iout] = otdm_uo_students(f_x, g_x, h_x, x, options);
xylim = [0,0,0,0];
%otdm_uo_plot(f_rosen, xk, gk, xylim);

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
