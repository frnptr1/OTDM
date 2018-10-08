% Procedure otdm_uo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, xk, alk, dk, gk, fk, iout] = otdm_uo_students(f_fun, g_fun, h_fun, x, options)
%
% iout = 0 : optimal solution found within tolerances.
% iout = 1 : too many iterations
% iout = 2 : ascent direction
% iout = 3 : linesearch failure.
%
% options( 1) = eps : zero tolerance
% options( 2) = maxiter
% options( 3) = max steplength coef.
% options( 4) = First strong Wolfe condition parameter c1
% options( 5) = Second strong Wolfe condition parameter c2
% options( 6) = max. iterations in BTLS.
% options( 7) = min diff. |al^i-al^(i-1)| in BTLS.
% options( 8) = flag fpr quadratic o.f..
% options(10) = 1 -> Gradient.
%               2 -> Conjugate Gradient, Fletcher-Reeves update.
%               3 -> Conjugate Gradient, Polak-Ribière update.
%               4 -> Quasi-Newton, DFP.
%               5 -> Quasi-Newton, BFGS.
%               6 -> Newton.
%               7 -> Modified Newton, Luenberger.
% options(11) = min eigval of the modified reduced Hessian, MNM-Luenberger.
%
if options( 1) == 0, options( 1) = 10^-6; end; eps        = options( 1);
if options( 2) == 0, options( 2) = 100.;  end; maxiter    = options( 2);
if options( 3) == 0, options( 3) =  1.01; end; almaxcoef  = options( 3);
if options( 4) == 0, options( 4) = 1e-04; end; c1         = options( 4);
if options( 5) == 0, options( 5) = .5;    end; c2         = options( 5);
if options( 6) == 0, options( 6) = 30;    end; maxiter_ls = options( 6);
if options( 7) == 0, options( 7) = 10^-3; end; eps_ls     = options( 7);
if options( 8) ~= 1, options( 8) = 0;     end; qof        = options( 8);
if options(10) == 0, options(10) = 1;     end; sdm        = options(10);
if options(11) == 0, options(11) = 0.1;   end; deltaLuen  = options(11);
%
fprintf('[otdm_uo_optimize] ::::::::::::::::::::::::::::::::::::::::::::::::::::::\n')
fprintf('[otdm_uo_optimize] Unconstrained optimization (OTDM/MIRI).\n')
fprintf('[otdm_uo_optimize] ::::::::::::::::::::::::::::::::::::::::::::::::::::::\n')
fprintf('[otdm_uo_optimize] Stopping conditions:\n')
fprintf('[otdm_uo_optimize]    Optimality tolerance  : options( 1) = %3.1e\n', options(1))
fprintf('[otdm_uo_optimize]    Max. iterations       : options( 2) = %i\n',    options(2))
fprintf('[otdm_uo_optimize] Baktracking linesearch:\n')
fprintf('[otdm_uo_optimize]    Max. stepl. coef.     : options( 3) = %i\n',    options(3))
fprintf('[otdm_uo_optimize]    SW1 condition, c1     : options( 4) = %3.1e\n', options(4))
fprintf('[otdm_uo_optimize]    SW2 condition, c2     : options( 5) = %3.1e\n', options(5))
fprintf('[otdm_uo_optimize]    Max. iterations       : options( 6) = %i\n',    options(6))
fprintf('[otdm_uo_optimize]    Min. progress         : options( 7) = %3.1e\n', options(7))
fprintf('[otdm_uo_optimize] Quadratic o.f.           : options( 8) = %i\n',  options(8))
fprintf('[otdm_uo_optimize] Optimization method:\n')
if     (options(10)==1)
    fprintf('[otdm_uo_optimize]    Gradient Method       : options(10) = 1\n')
elseif (options(10)==2)
    fprintf('[otdm_uo_optimize]    CGM/Fletcher-Reeves   : options(10) = 2\n')
elseif (options(10)==3)
    fprintf('[otdm_uo_optimize]    CGM/Polak-Ribière     : options(10) = 3\n')
elseif (options(10)==4)
    fprintf('[otdm_uo_optimize]    QNM/BFGS              : options(10) = 4\n')
elseif (options(10)==5)
    fprintf('[otdm_uo_optimize]    QNM/DFP               : options(10) = 5\n')
elseif (options(10)==6)
    fprintf('[otdm_uo_optimize]    Newton''s Method       : options(10) = 6\n')
elseif (options(10)==7)
    fprintf('[otdm_uo_optimize]    MNM-Luenberger        : options(10) = 6\n')
    fprintf('[otdm_uo_optimize]       delta              : options(11) = %3.1f\n', deltaLuen')
end
fprintf('[otdm_uo_optimize] Initial solution:\n');
fprintf('[otdm_uo_optimize]    x1=[ %7.3f, %7.3f]\n', x(1), x(2));
fprintf('[otdm_uo_optimize] Optimizaton:\n')
%
xk  = [];
alk = [];
dk  = [];
gk  = [];
fk  = [];
bk  = [];
Bbfgs = [];
sk = [];
yk = [];
bbfgs = [];
%
n = size(x,1);
i = [1:n];
k = 0;
la = eig(h_fun(x));
g   = g_fun(x);
f   = f_fun(x);

%if f is quadratic
if qof == 1
    Q=h_fun(x);
end
iout_bls = 0;
%
while norm(g) > eps & k < maxiter
    if mod(k,100) == 0
        fprintf('[otdm_uo_optimize] k    f(xk)            ||gk||        min(lak)  alk\n');
    end
    k = k + 1;
    xk = [xk,x];
    gk = [gk,g'];
    fk = [fk,f];
    la = eig(h_fun(xk(:,k)));
    d=zeros(n,1);
    if sdm == 1  % GM
        d = -g';
    elseif sdm == 2   % CGM/FR TODO
        if k == 1
            d = -g';
        else
            bk = (g*g')/(norm(gk(:,k-1))^2);
            d = -g' + bk*dk(:,k-1);
        end
    elseif sdm == 3  % CGM/PR
        d = -g';
    elseif sdm == 4  % QNM/BFGS TODO
        if k == 1
            Bbfgs = [Bbfgs, eye(n)];
            display(Bbfgs);
            d = -g';
        else
            sk = xk(:,k) - xk(:,k-1);
            display(xk(:,k));           
            yk = gk(:,k)' - gk(:,k-1)';
            bbfgs = Bbfgs(:,[k-1,k]) - (Bbfgs(:,[k-1,k])'*sk*transpose(sk)*Bbfgs(:,[k-1,k]))/(transpose(sk)*Bbfgs(:,[k-1,k])*sk) + (yk*transpose(yk))/(yk*sk);
            display(bbfgs);
            Bbfgs = [Bbfgs, bbfgs];
            
            d = -Bbfgs(:,[k-1,k])^(-1)*gk(:,k);
        end
    elseif sdm == 5  % QNM/D
        d = -g';
    elseif sdm == 6  % MN TODO
        d = -g';
    elseif sdm == 7  % MNM-Luenberger TODO
        d = -g';
    end
    % Descent direction condition:
    if g*d >=0
        dk=[dk,d];
        break
    end
    % Linesearch
    if k == 1
        almax = 1;
    elseif sdm < 4 % GM, CGM and QNM
        almax = almaxcoef*2*(fk(k)-fk(k-1))/(gk(:,k)'*d);        %NW (3.44) pag. 54
    else           % NM and MNM
        almax = min(1,almaxcoef*2*(fk(k)-fk(k-1))/(gk(:,k)'*d)); %NW (3.44) pag. 54
    end
    if qof == 1      % Quadratic o.f., exact linesearch
        al = -g*d/(d'*Q*d);
    else             % non-quadratic o.f., backtraking line-search with strong Wolfe cond.
        [al, iout_bls] = otdm_uo_backlinesearch(f_fun,g_fun,d,x,almax,c1,c2,maxiter_ls,eps_ls);
    end
    % Step
    alk = [alk,al];
    dk  = [dk,d];
    x   = x + al*d;
    g   = g_fun(x);
    f   = f_fun(x);
    %
    fprintf('[otdm_uo_optimize] %3d  %+10.8e  %8.6e  %+4.1e  %+4.1e\n', k, f_fun(xk(:,k)), norm(gk(:,k)), la(1), al);
    if iout_bls == 1
        fprintf('[otdm_uo_optimize] Warning linesearch: too much iterations, strong Wolfe conditions not guarenteed, convergence compromised.\n');
    elseif iout_bls == 2
        fprintf('[otdm_uo_optimize] Warning linesearch: stacked, strong Wolfe conditions not guarenteed, convergence compromised.\n');      
    end
    if (f-fk(k)) >= eps
        iout = 3;
        break
    end
end
if norm(g) <= eps
    iout = 0;
    k = k+1;
    xk = [xk,x];
    fk = [fk,f];
    gk = [gk,g'];
    la = eig(h_fun(xk(:,k)));
    fprintf('[otdm_uo_optimize] %3d  %+10.8e  %8.6e  %+4.1e \n', k, f_fun(xk(:,k)), norm(gk(:,k)), la(1));
    fprintf('[otdm_uo_optimize] Optimal solution within tolerance:\n');
    fprintf('[otdm_uo_optimize]    x*=[ %7.3f, %7.3f].\n', x(1), x(2));
elseif gk(:,k)'*dk(:,k) >= 0
    iout = 2;
    fprintf('[otdm_uo_optimize] %3d  %+10.8e  %8.6e\n', k, f_fun(xk(:,k)), norm(gk(:,k)));
    fprintf('[otdm_uo_optimize] FAILURE: dk is not a descent direction: gk·dk = %+4.1e, convergence lost.\n', gk(:,k)'*dk(:,k) )
elseif k == maxiter 
    iout = 1;
    k = k+1;
    xk = [xk,x];
    fk = [fk,f];
    gk = [gk,g'];
    la = eig(h_fun(xk(:,k)));
    fprintf('[otdm_uo_optimize] %3d  %+10.8e  %8.6e  %+4.1e\n', k, f_fun(xk(:,k)), norm(gk(:,k)), la(1));
    fprintf('[otdm_uo_optimize] maxiter = %3d exceeded.\n', maxiter)
elseif iout == 3
    k = k+1;
    xk = [xk,x];
    fk = [fk,f];
    gk = [gk,g'];
    fprintf('[otdm_uo_optimize] %3d  %+10.8e  %8.6e\n', k+1, f_fun(xk(:,k)), norm(gk(:,k)));
    fprintf('[otdm_uo_optimize] Linesearch failure: dk descent direction but fk+1 > fk.\n')
end
% End procedure otdm_uo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%