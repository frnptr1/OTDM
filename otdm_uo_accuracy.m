% Procedure otdm_uo_accuracy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function accuracy = otdm_uo_accuracy(xtr,ytr,wo, y_func) 
p = size(ytr,2);
model_response = round(y_func(xtr, wo));
accuracy = 0;
for i=1:p
    if ytr(i)==model_response(i)
        accuracy = accuracy+1;
    end
end
accuracy = accuracy*(100/p);
end
        

