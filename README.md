# APG-stochatic-optimal control problem

This repository contains the matlab implementaion of the APG algorithm to solve the 
stochatic optimal control problem for a system with both multiplicative and 
additive disturbances.


SYSTEM GENERATION:

This algorithm is tested on spring-mass systems. In this case, the 
dynamics can have both additive and multiplicative disturbances. The 
function "tree_generation_multiple.m" generates the system matrices and 
tree structure. There is an option include multiplicative disturbance in the  
system dynamics or not.

PRECONDITIONING:
Preconditioning the system inproves the convergence of the algorithm. The precondioning 
performed here is based on huristics. There are a couple of functions that do this. 
This section can be inproved. 

ALGORITHMS:

Before the main algorithm, offline matrices that are used in  
calculation of the dual gradient are computed as factor step. The 
function "GPAD_factor_step_smpc" caclulate the factor step for a system 
with multiplicative disturbances. There is another version for regular systems with 
additive disturbances ("GPAD_factor_regular").

Now the main algorithm is "GPAD_solve_smpc.m". This function performs the dual gradient 
projection algorithm. 

TEST FUNCTIONS:

The test file"test_APG.m" tests the APG algorithm. Here, a random of 100 points are 
generated and sampling points. It plots the number of iterations it require in solving the
algorithm.

In other test functions, the perfomace of APG algorithm in CPU, in GPU and other 
commercial solvers. To solve it in GPUs, we create a optimization problem in headers 
and solve it. For more information check the TB-GPAD repository. These test files 
have complicated options to check the perfomace over a 100 inital points with different 
scenarios and different solvers. 

