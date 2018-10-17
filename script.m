rules = [1,2,21,3,32];
lambdas = [0,0.1];
methods = [1,2,3,4,5];
linesearches = [1,2];

for irule = 1:size(rules,2)
    for ilambda = 1:size(lambdas,2)
        for imethod = 1:size(methods,2)
            for ils = 1:size(linesearch,2)
                otdm_uo_nn_students(rules(irule),methods(imethod),linesearches(ils));                
            end            
        end        
    end    
end
