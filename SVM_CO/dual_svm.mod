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


# VARIABLE
var lambda {i in 1..m};
var w{i in 1..n};


# OBJECTIVE FUNCTION
# we have omitted d_i and d_j from the multiplication because
# in any case the resault would have been 1
maximize obj_func : sum{i in 1..m} lambda[i] -
0.5*sum{i in 1..m, j in 1..m}lambda[i]*lambda[j]*y[i]*y[j]*
(sum{k in 1..n}x[i,k]*x[j,k]);

# CONSTRAINT
subject to constraint: sum{i in 1..m}lambda[i]*y[i]=0;
subject to lam1{i in 1..m}: lambda[i]>=0;
subject to lam2{i in 1..m}: lambda[i]<=nu;
subject to retrieve_w{i in 1..n}: w[i]=sum{j in 1..m}lambda[j]*y[j]*x[j,i];

#subject to retrieve_gamma{i in 1..m}: lambda[i]>0 ==> gamma = (sum{j in 1..n} w[j]*x[i,j])-(1/y[i]);
#subjecto to sp{i in 1..m} if (lambda[i]>0) support[i]=1 else 0;

