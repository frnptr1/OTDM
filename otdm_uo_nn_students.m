function otdm_uo_nn_students(prule,pmethod,plinesearch)    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OTDM / MIRI / FIB  F.-1javier Heredia https://gnom.upc.edu/heredia
    % Procedure otdm_uo_nn
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global xtr ytr la n p;
    fprintf('[otdm_uo_nn] ::::::::::::::::::::::::::::::::::::::::::::::::\n')
    fprintf('[otdm_uo_nn] Classification with neural networks (OTDM/MIRI).\n')
    fprintf('[otdm_uo_nn] ::::::::::::::::::::::::::::::::::::::::::::::::\n')

    %
    % TRAINING DATASET
    %
    rule=prule;
    training_seed = 123456;
    rng(training_seed);
    p_training=500;
    [xtr,ytr] = otdm_uo_nn_populate(p_training, rule);
    [n,p] = size(xtr);

    %
    % PARAMS
    %
    options=zeros(1,20);
    options(2) = 15000;
    options(9) = 1000;
    options(10) = pmethod;
    options(20) = plinesearch;

    % STARTING POINT
    w1 = zeros(n,1);

    % LAMBDA REGULARIZATION PARAM
    la = 0.0;


    fprintf('[otdm_uo_nn]    p_training  = %i\n', p_training)
    fprintf('[otdm_uo_nn]    classif. rule  = %i\n', rule)
    fprintf('[otdm_uo_nn]    training_seed  = %i\n', training_seed)
    fprintf('[otdm_uo_nn]    L2 reg. lambda = %f\n', la)
    if n==2
        fprintf('[otdm_uo_nn]    w1=[ %7.3f, %7.3f]\n', w1(1), w1(2));
    elseif n==3
        fprintf('[otdm_uo_nn]    w1=[ %7.3f, %7.3f, %7.3f]\n', w1(1), w1(2), w1(3));
    end

    % OUTPUT LAYER
    f = @(w) L(w);
    g = @(w) gL(w);

    % HESSIAN OF THE OB1jECTIVE FUNCTION
    h = @(w) eye(n);

    if n==2
        % CHECKING OBJECTIVE FUNCTION AND GRADIENT
        w_test  = [2,3];
        w_test2 = [0.5, 0.5];

        fprintf('[otdm_uo_nn] ::::::::::::::::::::::::::::::::::::::::::::::::\n')
        fprintf('[otdm_uo_nn] OBJECTIVE FUNCTION AND GRADIENT CHECKING\n')
        fprintf('[otdm_uo_nn] ::::::::::::::::::::::::::::::::::::::::::::::::\n')
        fprintf("[otdm_uo_nn]    w1    : [2,3]' \n");
        fprintf("[otdm_uo_nn]    L(w1) : %7.4f \n", f(w_test));
        fprintf("[otdm_uo_nn]    g(w1) : [%7.4f, %7.4f]  \n", g(w_test));
        fprintf("[otdm_uo_nn]    w2    : [0.5,0.5]' \n");
        fprintf("[otdm_uo_nn]    L(w2) : %7.4f \n", f(w_test2));
        fprintf("[otdm_uo_nn]    g(w2) : [%7.4f, %7.4f]  \n", g(w_test2));
    end

    %
    % Optimization
    %
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
    training_accuracy= otdm_uo_accuracy(wo, xtr, ytr);
    fprintf('[otdm_uo_nn] ::::::::::::::::::::::::::::::::::::::::::::::::\n')
    fprintf('[otdm_uo_nn]    training_accuracy = %4.1f\n', training_accuracy)

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
    test_accuracy= otdm_uo_accuracy(wo, xte, yte);
    fprintf('[otdm_uo_nn]    p_test            = %i\n', p_test)
    fprintf('[otdm_uo_nn]    test_seed         = %i\n', test_seed)
    fprintf('[otdm_uo_nn]    test_accuracy     = %4.1f\n', test_accuracy)
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End Procedure otdm_uo_nn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



