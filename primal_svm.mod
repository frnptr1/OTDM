# PARAMETERS
# dataset dimension
param m;
# num variable for each datapoint
param n;
# nu
param nu;
# components of dataset
param x{i in 1..m, j in 1..n};
# class coming from data file
param y{i in 1..m};
# dataset
set dataset;

# VARIABLES
# variables for optimization task
var w{i in 1..n};
# intercept of hyperplane
var gamma;
# slack variables
var s{i in 1..m};

# OBJECTIVE FUNCTION
# objective function
minimize obj_func: 0.5*sum{i in 1..n} (w[i]^2) + nu*sum{i in 1..m}s[i];

# CONSTRAINTS
# fix the datapoint, so for each point..
subject to constr{i in 1..m}:
y[i] * ( sum{j in 1..n} ( w[j] * x[i,j] )-gamma) >= 1-s[i];

subject to positiveness{i in 1..m}:
s[i]>=0;
















