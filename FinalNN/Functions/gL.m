function result = gL(w)
    global xtr ytr n p la; 
    result = zeros(n,1)';
    
    
    for i = 1:n
        temp = 0;
        for j = 1:p
            temp = temp  + 2 * (y_func(xtr(:,j),w) - ytr(j)) * (y_func(xtr(:,j),w) - y_func(xtr(:,j),w)^2) * activation(xtr(i,j));
        end
        temp = temp + la*w(i);
        result(i) = temp;
    end

end



    %g = @(w) [ sum(2*(y_func(xtr,w) - ytr) * (y_func(xtr,w) - y_func(xtr,w).^2)' * activation( xtr(1,:) ) + la*w(1)), sum(2*(y_func(xtr,w) - ytr) * (y_func(xtr,w) - y_func(xtr,w).^2)' * activation( xtr(2,:) ) + la*w(2))];
