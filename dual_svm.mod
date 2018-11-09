# PARAMETERS
# dataset dimension
param m;
# num variable for each datapoint
param n;
param nu;
# components of dataset
param x{i in 1..m, j in 1..n};
# class coming from data file
param y{i in 1..m};
# dataset
set dataset;

# VARIABLE
var lambda {i in 1..m};
var w{i in 1..n};

# OBJECTIVE FUNCTION
# we have omitted d_i and d_j from the multiplication because
# in any case the resault would have been 1
maximize obj_func : sum{i in 1..m} lambda[i] -
0.5*sum{i in 1..m, j in 1..n}lambda[i]*lambda[j]*
(sum{k in 1..n}x[i,k]*x[j,k]);

# CONSTRAINT
subject to constraint: sum{i in 1..m}lambda[i]*y[i]=0;
subject to lam1{i in 1..m}: lambda[i]>=0;
subject to lam2{i in 1..m}: lambda[i]<=nu ;
subject to retrieve_w{i in 1..n}: w[i]>=sum{j in 1..m}lambda[i]*y[i]*x[j,i];
