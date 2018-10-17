function result = L(w)
    global xtr ytr la n p;
    result = 0;
    for j = 1:p
        result = result + (y_func(xtr(:,j),w) - ytr(j))^2;
        
        if (j <= n)
            result = result + (la/2)*w(j)^2;
        end
    end
    
end