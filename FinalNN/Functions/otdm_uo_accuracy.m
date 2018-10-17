% Procedure otdm_uo_accuracy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function accuracy = otdm_uo_accuracy(wo, xte, yte) 
p = length(yte);
accuracy = 0;
for i=1:p
    if yte(i)==round(y_func(xte(:,i),wo))
        accuracy = accuracy+1;
    end
end
accuracy = accuracy*(100/p);
end




        

