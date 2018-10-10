% Procedure otdm_uo_accuracy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function accuracy = otdm_uo_accuracy(wo) 
global xtr ytr p;
accuracy = 0;
for i=1:p
    if ytr(i)==round(y_func(xtr(:,i),wo))
        accuracy = accuracy+1;
    end
end
accuracy = accuracy*(100/p);
end
        

