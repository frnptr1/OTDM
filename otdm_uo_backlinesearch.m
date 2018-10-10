function [alphas,iout] = otdm_uo_backlinesearch(f,g,d,x0,alpham,c1,c2,maxiter,eps)
% function alphas = strongwolfe(f,d,x0,alpham)
% Line search algorithm satisfying strong Wolfe conditions
% Algorithms 3.5 on pages 60-61 in Nocedal and Wright
% MATLAB code by Kartik Sivaramakrishnan
% Last modified: January 27, 2008
%
% F.-Javier Heredia, September 2018 <fjh...>
% g,c1,c2,maxiter, eps
% iout = 1 : too much iterations
% iout = 2 : stacked, alpha_[i]=alpha^[i-1]

alpha0 = 0;
alphap = alpha0;
%<fjh
iout = 0;
if c1 == 0
    c1 = 1e-4;
end
if c2==0
    c2 = 0.5;
end
%alphax = alpham*rand(1);
alphax = alpham;
%[fx0,gx0] = feval(f,x0,d);
fx0 = f(x0);
gx0 = g(x0)*d;
%>
fxp = fx0;
gxp = gx0;
i=1;
% alphap is alpha_{i-1}
% alphax is alpha_i
while (1 ~= 2 && i < maxiter)
  %<fjh
    if abs(alphap-alphax) < eps
      iout = 2;
      alphas = alphax;
      return
    end
  %>
    xx = x0 + alphax*d;
  %<fjh
  %[fxx,gxx] = feval(f,xx,d);
    fxx = f(xx);
    gxx = g(xx)*d;
  %>
    if (fxx > fx0 + c1*alphax*gx0) || ((i > 1) && (fxx >= fxp))
        [alphas,iout_zoom] = zoom(f,g,x0,d,alphap,alphax,c1,c2,eps);
        %<fjh
        if iout_zoom == 2
            iout = 2;
        end
        %>
        return;
    end
  if abs(gxx) <= -c2*gx0,
    alphas = alphax;
    return;
  end
  if gxx >= 0,
    [alphas,iout_zoom] = zoom(f,g,x0,d,alphax,alphap,c1,c2,eps);
    %<fjh
    if iout_zoom == 2
        iout = 2;
    end
    %>
    return;
  end
  alphap = alphax;
  fxp = fxx;
  gxp = gxx;
  alphax = alphax + (alpham-alphax)*rand(1);
  i = i+1;
end
if i==iter_btls
    iout = 1;
    alphas = alphax;
end

function [alphas,iout] = zoom(f,g,x0,d,alphal,alphah,c1,c2,eps)
% function alphas = zoom(f,g,x0,d,alphal,alphah)
% Algorithm 3.6 on page 61 in Nocedal and Wright
% MATLAB code by Kartik Sivaramakrishnan
% Last modified: January 27, 2008
% F.-Javier Heredia, September 2018 <fjh...>

%<fjh
%[fx0,gx0] = feval(f,x0,d);
fx0 = f(x0);
gx0 = g(x0)*d;
iout = 0;
%>
while (1~=2),
    %<fjh
    if abs(alphal-alphah) < eps
      iout = 2;
      alphas = alphax;
      return
    end
    %>

   alphax = 1/2*(alphal+alphah);
   xx = x0 + alphax*d;
   %<fjh
   %[fxx,gxx] = feval(f,xx,d);
   fxx = f(xx);
   gxx = g(xx)*d;
   %>
   xl = x0 + alphal*d;
   %<fjh
   %fxl = feval(f,xl,d);
   fxl = f(xl);
   %>
   if ((fxx > fx0 + c1*alphax*gx0) || (fxx >= fxl)),
      alphah = alphax;
   else
      if abs(gxx) <= -c2*gx0,
        alphas = alphax;
        return;
      end
      if gxx*(alphah-alphal) >= 0,
        alphah = alphal;
      end
      alphal = alphax;
   end
end 