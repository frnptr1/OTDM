function result = gL(w)
    global xtr ytr n p la; 
    result = zeros(3,1);
    temp = 0;
    
    for i = 1:n
        for j = 1:p
            temp = temp  + 2 * (y_func(xtr(:,j),w) - ytr(j))*(y_func(xtr(:,j),w) - y_func(xtr(:,j),w)^2)*activation(xtr(i,j)) + la*w(i);
        end
        result(i) = temp;
    end

end