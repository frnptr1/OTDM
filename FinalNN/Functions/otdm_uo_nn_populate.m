function [xtr, ytr] = otdm_uo_nn_populate(p, rule)

if rule ==1
    xtr= round(rand(2,p));
    ytr= zeros(p,1)';   
    for j=1:p
        if xtr(2,j) > xtr(1,j)
        ytr(j) = 1;
        end
    end
end


if rule ==2
    xtr= round(rand(2,p));
    ytr= zeros(p,1)';
    for j=1:p
        if xtr(2,j) == xtr(1,j)
            ytr(j) = 1;
        end
    end
    xtr(xtr==0) = xtr(xtr==0)-1;
end


if rule ==21
    xtr= round(rand(3,p));
    ytr= zeros(p,1)';
    xtr(xtr==0) = xtr(xtr==0)-1;
    for j=1:p
        if xtr(2,j) == xtr(1,j)
            ytr(j) = 1;            
        end
        xtr(3,j) = xtr(2,j)*xtr(1,j);
    end
end

if rule == 3
    a = -1;
    b = 1;
    xtr = round((b-a).*rand(2,p)+a);
    ytr = zeros(p,1);
    for j=1:p
        if xtr(2,j)+xtr(1,j) >0
            ytr(j) = 1;
        end
    end
end

end