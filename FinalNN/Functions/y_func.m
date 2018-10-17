function sum = y_func(x,w)
    global n;
    sum = 0;
    for i = 1:n
        sum = sum  +w(i)*activation(x(i));        
    end
    
    sum = activation(sum);
end