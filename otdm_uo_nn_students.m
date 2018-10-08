%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OTDM / MIRI / FIB  F.-1javier Heredia https://gnom.upc.edu/heredia
% Procedure otdm_uo_nn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
global xtr ytr la n p;
fprintf('[otdm_uo_nn] ::::::::::::::::::::::::::::::::::::::::::::::::\n')
fprintf('[otdm_uo_nn] Classification with neural networks (OTDM/MIRI).\n')
fprintf('[otdm_uo_nn] ::::::::::::::::::::::::::::::::::::::::::::::::\n')

%
% Training data set
%
rule=21;
training_seed = 123456;
rng(training_seed);
p_training=500;
[xtr,ytr] = otdm_uo_nn_populate(p_training, rule);
[n,p] = size(xtr);
fprintf('[otdm_uo_nn]       p_training  = %i\n', p_training)
fprintf('[otdm_uo_nn]    classif. rule  = %i\n', rule)
fprintf('[otdm_uo_nn]    training_seed  = %i\n', training_seed)

%
% Optimization
%


la = .1;
fprintf('[otdm_uo_nn]    L2 reg. lambda = %f\n', la)

% ACTIVATION FUNCTION
activation = @(x) (1+exp(-x)).^(-1);

% OUTPUT LAYER
%y_func = @(x,w)  (1+exp(-sum(activation(x).*w))).^(-1);
y_func = @(x,w)  activation(activation(x)'*w)';

% LOSS FUNCTION REGULARIZED WITH LAMBDA
f = @(w) sum((y_func(xtr, w) - ytr).^2) + (la/2)*sum(w.^2);
    

% GRADIENT OF THE OB1jECTIVE FUNCTION
if n==2
    g = @(w) [ sum(2*(y_func(xtr,w) - ytr) * (y_func(xtr,w) - y_func(xtr,w).^2)' * activation( xtr(1,:) ) + la*w(1)), sum(2*(y_func(xtr,w) - ytr) * (y_func(xtr,w) - y_func(xtr,w).^2)' * activation( xtr(2,:) ) + la*w(2))];
elseif n==3
    g = @(w) [ sum(2*(y_func(xtr,w) - ytr) * (y_func(xtr,w) - y_func(xtr,w).^2)' * activation( xtr(1,:) ) + la*w(1)), sum(2*(y_func(xtr,w) - ytr) * (y_func(xtr,w) - y_func(xtr,w).^2)' * activation( xtr(2,:) ) + la*w(2)), sum(2*(y_func(xtr,w) - ytr) * (y_func(xtr,w) - y_func(xtr,w).^2)' * activation( xtr(3,:) ) + la*w(3))];
end

% HESSIAN OF THE OB1jECTIVE FUNCTION
h = @(w) eye(n);

% PARAMS
options=zeros(1,20);
w1 = zeros(n,1);
%w1 = w1';
options( 2) = 1000;
options(10) = 1;
options(20) = 1;

if n==2
    fprintf('[otdm_uo_nn]    w1=[ %7.3f, %7.3f]\n', w1(1), w1(2));
elseif n==3
    fprintf('[otdm_uo_nn]    w1=[ %7.3f, %7.3f, %7.3f]\n', w1(1), w1(2), w1(3));
end

[wo, xk, alk, dk, gk, fk, iout, options] = otdm_uo_students_new(f, g, h, w1, options);





%
% Plot
%
if n==2
    xylim = [-50,50,-50,50];
    otdm_uo_plot(f, xk, gk, xylim);
end

%
% Training accuracy
%
training_accuracy= otdm_uo_accuracy(xtr,ytr,wo, y_func);
fprintf('[otdm_uo_nn] training_accuracy = %4.1f\n', training_accuracy)

%
% Test dataset
%
p_test = 100*p_training;
test_seed = 789101;
rng(test_seed);
[xte,yte] = otdm_uo_nn_populate(p_test, rule);
%
% Test accuracy
%
test_accuracy= otdm_uo_accuracy(xte,yte,wo, y_func);
fprintf('[otdm_uo_nn]            p_test = %i\n', p_test)
fprintf('[otdm_uo_nn]         test_seed = %i\n', test_seed)
fprintf('[otdm_uo_nn]     test_accuracy = %4.1f\n', test_accuracy)
%
% Log file.
%

fileID = fopen('otdm_uo_nn.out','a');
sdm = options(10);
if sdm ==1
    csdm = 'GM';
elseif sdm==2
    csdm = 'CGM-FR';
elseif sdm==3
    csdm = 'CGM-PR';
elseif sdm==4
    csdm = 'BFGS';
elseif sdm==5
    csdm = 'DFP';
end
[kk,iter] = size(gk);
if n==2
    fprintf(fileID,' r tr_p tr_seed    la     w1(1)     w1(2) ls  c2    sdm iout iter     w*(1)     w*(2)        L* tr_acc te_acc   te_q te_seed\n');
    fprintf(fileID,'%2i %4i %i6 %5.3f %+5.2e %+5.2e  %1i %3.1f %6s   %2i %4i %+5.2e %+5.2e %5.3e  %5.1f  %5.1f %6i %i6 \n', rule, p_training, training_seed, la, w1(1), w1(2), options(20), options(5), csdm,iout, iter, wo(1),wo(2), f(wo), training_accuracy, test_accuracy, p_test, test_seed);
elseif n==3
    fprintf(fileID,' r tr_p tr_seed    la     w1(1)     w1(2)     w1(3) ls  c2    sdm iout iter     w*(1)     w*(2)     w*(3)        L* tr_acc te_acc   te_q te_seed\n');
    fprintf(fileID,'%2i %4i %i6 %5.3f %+5.2e %+5.2e %+5.2e  %1i %3.1f %6s   %2i %4i %+5.2e %+5.2e %+5.2e %5.3e  %5.1f  %5.1f %6i %i6 \n', rule, p_training, training_seed, la, w1(1), w1(2), w1(3), options(20), options(5), csdm,iout, iter, wo(1),wo(2),wo(3), f(wo), training_accuracy, test_accuracy, p_test, test_seed);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End Procedure otdm_uo_nn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

