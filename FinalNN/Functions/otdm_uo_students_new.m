% Procedure otdm_uo_solve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, xk, alk, dk, gk, fk, iout, options] = otdm_uo_solve(f_fun, g_fun, h_fun, x, options)
%
% iout = 0 : optimal solution found within tolerances.
% iout = 1 : too many iterations
% iout = 2 : ascent direction
% iout = 3 : linesearch failure.
%
% options( 1) = eps : zero tolerance
% options( 2) = maxiter
% options(20) = 1  -> otdm_uo_backlinesearch
%             = 2  -> fminbnd
% options( 3) = max steplength coef.
% options( 4) = First strong Wolfe condition parameter c1
% options( 5) = Second strong Wolfe condition parameter c2
% options( 6) = max. iterations in BTLS.
% options( 7) = min diff. |al^i-al^(i-1)| in BTLS.
% options( 8) = flag fpr quadratic o.f.
% options( 9) = output frequancy..
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
if options(20) == 0, options(20) =  1;    end; lsproc     = options(20);
if options( 3) == 0, options( 3) =  1.01; end; almaxcoef  = options( 3);
if options( 4) == 0, options( 4) = 1e-04; end; c1         = options( 4);
if options( 5) == 0, options( 5) = .5;    end; c2         = options( 5);
if options( 6) == 0, options( 6) = 30;    end; maxiter_ls = options( 6);
if options( 7) == 0, options( 7) = 10^-3; end; eps_ls     = options( 7);
if options( 8) ~= 1, options( 8) = 0;     end; qof        = options( 8);
if options( 9) == 0, options( 9) = 1;     end; outfreq    = options( 9);
if options(10) == 0, options(10) = 1;     end; sdm        = options(10);
if options(11) == 0, options(11) = 0.1;   end; deltaLuen  = options(11);
%
fprintf('[otdm_uo_solve] ::::::::::::::::::::::::::::::::::::::::::::::::::::::\n')
fprintf('[otdm_uo_solve] Unconstrained optimization (OTDM/MIRI).\n')
fprintf('[otdm_uo_solve] ::::::::::::::::::::::::::::::::::::::::::::::::::::::\n')
fprintf('[otdm_uo_solve] Stopping conditions:\n')
fprintf('[otdm_uo_solve]    Optimality tolerance  : options( 1) = %3.1e\n', options(1))
fprintf('[otdm_uo_solve]    Max. iterations       : options( 2) = %i\n',    options(2))
fprintf('[otdm_uo_solve] Baktracking linesearch:\n')
if     (options(20)==1)
    fprintf('[otdm_uo_solve]    bckls with SW cond.   : options(20) = 1\n')
elseif   (options(20)==2)
    fprintf('[otdm_uo_solve]    fminbnd               : options(20) = 2\n')
end
fprintf('[otdm_uo_solve]    linesearch procedure  : options( 3) = %5.3f\n', options(3))
fprintf('[otdm_uo_solve]    Max. steplength coef. : options( 3) = %5.3f\n', options(3))
fprintf('[otdm_uo_solve]    SW1 condition, c1     : options( 4) = %3.1e\n', options(4))
fprintf('[otdm_uo_solve]    SW2 condition, c2     : options( 5) = %3.1e\n', options(5))
fprintf('[otdm_uo_solve]    Max. iterations       : options( 6) = %i\n',    options(6))
fprintf('[otdm_uo_solve]    Min. progress         : options( 7) = %3.1e\n', options(7))
fprintf('[otdm_uo_solve] Quadratic o.f.           : options( 8) = %i\n',    options(8))
fprintf('[otdm_uo_solve] Output frequency         : options( 9) = %i\n',    options(9))
fprintf('[otdm_uo_solve] Optimization method:\n')
if     (options(10)==1)
    fprintf('[otdm_uo_solve]    Gradient Method       : options(10) = 1\n')
elseif (options(10)==2)
    fprintf('[otdm_uo_solve]    CGM/Fletcher-Reeves   : options(10) = 2\n')
elseif (options(10)==3)
    fprintf('[otdm_uo_solve]    CGM/Polak-Ribière     : options(10) = 3\n')
elseif (options(10)==4)
    fprintf('[otdm_uo_solve]    QNM/BFGS              : options(10) = 4\n')
elseif (options(10)==5)
    fprintf('[otdm_uo_solve]    QNM/DFP               : options(10) = 5\n')
elseif (options(10)==6)
    fprintf('[otdm_uo_solve]    Newton''s Method       : options(10) = 6\n')
elseif (options(10)==7)
    fprintf('[otdm_uo_solve]    MNM-Luenberger        : options(10) = 6\n')
    fprintf('[otdm_uo_solve]       delta              : options(11) = %3.1f\n', deltaLuen')
end
fprintf('[otdm_uo_solve] Initial solution:\n');
fprintf('[otdm_uo_solve]    x1=[ %7.3f, %7.3f]\n', x(1), x(2));
fprintf('[otdm_uo_solve] Optimizaton:\n')
%
xk  = [];
alk = [];
dk  = [];
gk  = [];
fk  = [];
Bbfgs = {};
%
n = size(x,1);
i = [1:n];
k = 0;
la = eig(h_fun(x));
f  = f_fun(x);
g  = g_fun(x);
restart = n;
if qof == 1
    Q=h_fun(x);
end
iout = 3;
iout_bls = 0;
%
fprintf('[otdm_uo_solve]     k  f(xk)            ||gk||        min(lak)  alk\n');
while norm(g) > eps && k < maxiter
    if mod(k,25*outfreq) == 0 && k~=0
        fprintf('[otdm_uo_solve]     k  f(xk)            ||gk||        min(lak)  alk\n');
    end
    k = k + 1;
    xk = [xk,x];
    gk = [gk,g'];
    fk = [fk,f];
    la = eig(h_fun(xk(:,k)));
    d=zeros(n,1);
    if sdm == 1  % GM
        d = -g';
    
    elseif sdm ==2   % CGM/FR
               
        % first mandatory iteration
        if k == 1
            d = -g';     
         %zoutendijk condition, each 2*k iteration the direction is
         %restared with the steepest descent method
         elseif mod(k,n)==0 
            %restart = 2*k;
            d = -g';    

        else 
            % evaluation of beta(Fletcher-Reeves)    
            %beta_k = (gk(:,k)'*gk(:,k)) / norm(gk(:,k-1)')^2;
            beta_k = (g*g') / norm(gk(:,k-1)')^2;
            %d = -gk(:,k)+(beta_k*dk(:,k-1));
            d = -g'+(beta_k*dk(:,k-1));
        end
         
    elseif sdm == 3  % CGM/PR
                
        % first mandatory iteration
        if k == 1
            d = -g';     
         %zoutendijk condition, each 2*k iteration the direction is
         %restared with the steepest descent method
         elseif mod(k,n)==0
            %restart = 2*k;
            d = -g';    

        else 
            % evaluation of beta(Fletcher-Reeves)    
            beta_k = (gk(:,k)'* (gk(:,k)'-gk(:,k-1)')') /(norm(gk(:,k-1)'))^2;
            d = -gk(:,k)+(beta_k*dk(:,k-1));            
        end
         
    elseif sdm == 4  % QNM/BFGS
        
        if k == 1
            % bfgs matrix is just the identity matrix
            Bbfgs{k} = eye(n);
            % search direction is computed using steepest descent method
            d = -g';
            
            
        else
            % sk and yk evaluation
            sk = xk(:,k) - xk(:,k-1);      
            yk = gk(:,k) - gk(:,k-1);
            
            % second term of bfgs matrix
            B = ((Bbfgs{k-1}*(sk*sk')*Bbfgs{k-1})/(sk'*Bbfgs{k-1}*sk));
            % third term of bfgs matrix
            C = (yk*(yk)')/(yk'*sk);
            % evaluation of bfgs matrix
            bbfgs =  Bbfgs{k-1} - B+ C;
            % adding the matrix just created in the cell array
            Bbfgs{k} = bbfgs;
            % evaluating the search direction
            d = -inv(bbfgs)*gk(:,k);
        end

    elseif sdm == 5  % QNM/DFP
        
        if k == 1
            H{k} = eye(n);
            d = -g';            
        else        
            sk = xk(:,k) - xk(:,k-1);              
            yk = gk(:,k) - gk(:,k-1);           

            C = (sk*sk')/(yk'*sk);
            B = (H{k-1}*(yk*yk')*H{k-1})/(yk'*H{k-1}*yk);
            h = H{k-1}-B+C;

            H{k} = h;

            d = -h*gk(:,k);
        end
    elseif sdm == 6  % NM
        
        hessian_term = inv(h_fun(x));          
        d = -hessian_term*gk(:,k);
        
    elseif sdm == 7  % MNM-Luenberger
       
       hessian_matrix = h_fun(x); 
       [L,D] = ldl(hessian_matrix);
       [~,n] = size(D);
       DL=D+(deltaLuen-min(min(D)))*eye(n);
       HL=L*DL*L';
       d=-inv(HL)*g';
       
    end
    % Descent direction condition:
    if g*d >=0
        dk=[dk,d];
        break
    end
    % Linesearch
    if k == 1
        almax = 1;
    elseif sdm < 6 % GM, CGM and QNM
        almax = 2*(fk(k)-fk(k-1))/(gk(:,k)'*d);        %NW (3.44) pag. 54
    else           % NM and MNM
        almax = min(1,almaxcoef*2*(fk(k)-fk(k-1))/(gk(:,k)'*d)); %NW (3.44) pag. 54
    end
    if qof == 1      % Quadratic o.f., exact linesearch
        al = -g*d/(d'*Q*d);
    else             % non-quadratic o.f., backtraking line-search with strong Wolfe cond.
        if lsproc == 1
            [al, iout_bls] = otdm_uo_backlinesearch(f_fun,g_fun,d,x,almax,c1,c2,maxiter_ls,eps_ls);
        elseif lsproc == 2
            phi=@(al) f_fun(x+al*d);
            al = fminbnd(phi,0,almax);
        end
    end
    % Step
    alk = [alk,al];
    dk  = [dk,d];
    x   = x + al*d;
    g   = g_fun(x);
    f   = f_fun(x);
    %
    if mod(k,outfreq) == 0
        fprintf('[otdm_uo_solve] %5d  %+10.8e  %8.6e  %+4.1e  %+4.1e\n', k, f_fun(xk(:,k)), norm(gk(:,k)), la(1), al);
    end
    if iout_bls == 1
        fprintf('[otdm_uo_solve] WARNING linesearch: too much iterations, strong Wolfe conditions not guaranteed, convergence compromised.\n');
    elseif iout_bls == 2
        fprintf('[otdm_uo_solve] WARNING linesearch: stacked, strong Wolfe conditions not guaranteed, convergence compromised.\n');      
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
    fprintf('[otdm_uo_solve] %5i  %+10.8e  %8.6e  %+4.1e \n', k, f_fun(xk(:,k)), norm(gk(:,k)), la(1));    
    if n==2
        fprintf('[otdm_uo_solve] Optimal solution within tolerance:\n');
        fprintf('[otdm_uo_solve]    x*=[ %7.3f, %7.3f].\n', x(1), x(2));
    elseif n==3
        fprintf('[otdm_uo_solve] Optimal solution within tolerance:\n');
        fprintf('[otdm_uo_solve]    x*=[ %7.3f, %7.3f, %7.3f].\n', x(1), x(2), x(3));
    end


elseif gk(:,k)'*dk(:,k) >= 0
    iout = 2;
    fprintf('[otdm_uo_solve] %5i  %+10.8e  %8.6e\n', k, f_fun(xk(:,k)), norm(gk(:,k)));
    fprintf('[otdm_uo_solve] FAILURE: dk is not a descent direction: gk·dk = %+4.1e, convergence lost.\n', gk(:,k)'*dk(:,k) )
    if n==2
        fprintf('[otdm_uo_solve] Last iterate:\n');
        fprintf('[otdm_uo_solve]    x*=[ %7.3f, %7.3f].\n', x(1), x(2));
    elseif n==3
        fprintf('[otdm_uo_solve] Last iterate:\n');
        fprintf('[otdm_uo_solve]    x*=[ %7.3f, %7.3f, %7.3f].\n', x(1), x(2), x(3));
    end
    
 
elseif k == maxiter 
    iout = 1;
    k = k+1;
    xk = [xk,x];
    fk = [fk,f];
    gk = [gk,g'];
    la = eig(h_fun(xk(:,k)));
    fprintf('[otdm_uo_solve] %5i  %+10.8e  %8.6e  %+4.1e\n', k, f_fun(xk(:,k)), norm(gk(:,k)), la(1));
    fprintf('[otdm_uo_solve] maxiter = %3d exceeded.\n', maxiter)
    if n==2
        fprintf('[otdm_uo_solve] Last iterate:\n');
        fprintf('[otdm_uo_solve]    x*=[ %7.3f, %7.3f].\n', x(1), x(2));
    elseif n==3
        fprintf('[otdm_uo_solve] Last iterate:\n');
        fprintf('[otdm_uo_solve]    x*=[ %7.3f, %7.3f, %7.3f].\n', x(1), x(2), x(3));
    end
    
   
elseif iout == 3
    k = k+1;
    xk = [xk,x];
    fk = [fk,f];
    gk = [gk,g'];
    fprintf('[otdm_uo_solve] %5i  %+10.8e  %8.6e\n', k+1, f_fun(xk(:,k)), norm(gk(:,k)));
    fprintf('[otdm_uo_solve] Linesearch failure: dk descent direction but fk+1 > fk.\n')
    if n==2
        fprintf('[otdm_uo_solve] Last iterate:\n');
        fprintf('[otdm_uo_solve]    x*=[ %7.3f, %7.3f].\n', x(1), x(2));
    elseif n==3
        fprintf('[otdm_uo_solve] Last iterate:\n');
        fprintf('[otdm_uo_solve]    x*=[ %7.3f, %7.3f, %7.3f].\n', x(1), x(2), x(3));
    end
end
% End procedure otdm_uo_solve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%