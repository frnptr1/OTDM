function [xtr, ytr] = otdm_uo_nn_populate(p, rule)

training_seed= 123456;
rng(training_seed);
%p= 10;
xtr= round(rand(2,p));
ytr= zeros(p,1)';
for j=1:p
    if xtr(2,j) == xtr(1,j)
    ytr(j) = 1;
    end
end

xtr(xtr==0) = xtr(xtr==0)-1;

end
