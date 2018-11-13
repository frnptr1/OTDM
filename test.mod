
param z;
param y_test{i in 1..z};
param x_test{i in 1..z, j in 1..n};
var results {i in 1..z} = if (sum {j in 1..n} x_test[i,j] * w[j] >= gamma) then 1 else -1;
var misc = ((sum{k in 1..z} abs(results[k] - y_test[k])/2)/z);