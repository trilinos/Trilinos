// Simple Power method algorithm
double power_method(const Epetra_CrsMatrix& A) {  
  // variable needed for iteration
  double lambda = 0.0;
  int niters = A.RowMap().NumGlobalElements()*10;
  double tolerance = 1.0e-10;
  // Create vectors
  Epetra_Vector q(A.RowMap());
  Epetra_Vector z(A.RowMap());
  Epetra_Vector resid(A.RowMap());
  // Fill z with random Numbers
  z.Random();
  // variable needed for iteration
  double normz;
  double residual = 1.0 + tolerance;
  int iter = 0;
  while (iter < niters && residual > tolerance) {
    z.Norm2(&normz); // Compute 2-norm of z
    q.Scale(1.0/normz, z);
    A.Multiply(false, q, z); // Compute z = A*q
    q.Dot(z, &lambda); // Approximate maximum eigenvalue
    if (iter%10==0 || iter+1==niters) {
      resid.Update(1.0, z, -lambda, q, 0.0); // Compute A*q - lambda*q
      resid.Norm2(&residual);
      cout << "Iter = " << iter << "  Lambda = " << lambda 
	   << "  Two-norm of A*q - lambda*q = " 
	   << residual << endl;
    } 
    iter++;
  }
  return(lambda);
}
