In Step 2 of the tutorial, we will do the following:

1. Add a Dirichlet boundary condition by subclassing the BCStrategy. Add this to the BCFactory and apply a simple (constant) dirichlet condition. Turn off the Helmholtz operator to use a Laplacian

2. Add a nonlinear source term (say u*u). Modify the solve loop to do a Newton loop using the model evaluator.
