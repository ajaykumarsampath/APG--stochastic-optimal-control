# APG-stochatic-optimal control problem

This repository contains the matlab implementaion of the APG algorithm to solve the 
stochatic optimal control problem for a system with both multiplicative and 
additive disturbances. 

SYSTEM GENERATION:

To test this algorithms are tested on spring-mass systems. In this case, the 
dynamics can have both additive and multiplicative disturbances. The 
function "tree_generation_multiple.m" which generate the system matrices and 
tree structure. There is a option to add multiplicative disturbance in the  
system dynamics or not. 


ALGORITHMS:

Before the main algorithm, matrices are calculated offline that will be used in 
caculation of the dual gradient. This step is given as factor step. The 
function "GPAD_factor_step_smpc" caclulate the factor step for a system 
with multiplicative disturbances. There is another version for regular systems with 
additive disturbances ("GPAD_factor_regular").

Now the main algorithm is "GPAD_solve_smpc.m". This fucntion performs the dual gradient 
projection algorithm. 

TEST FUNCTIONS:

The test-function compare the perfomace of APG algorithm in CPU, in GPU and other 
commercial solvers. To solve it in GPUs, we create a optimization problem in headers 
and solve it. For more information check the TB-GPAD repository. These test files 
have complicated options to check the perfomace over a 100 inital points with different 
scenarios and different solvers.  

