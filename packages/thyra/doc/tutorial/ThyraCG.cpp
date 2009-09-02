#include "Vector.hpp"
#include "LinearOperator.hpp"

using namespace Thyra;

enum convStat = ( CONVERGENCE, MAX_ITERATIONS_EXCEEDED);


/* Definition of the class  */
class ThryaCG
{
public: 

  /* empty constructor */
  ThyraCG(){;}


  /* Main run method
  convStat = run(const LinearOperator<double>& A, 
		 const Vector<double>& b, 
		 Vector<double>& x) const;
};



/* The implementation of the run method for the class 
   ThyraCG  */

  convStat  ThyraCG::run(const LinearOperator<double>& A,  
		      const Vector<double>& b, 
		      const double& tol, 
		      const int* maxit, 
		      Vector<double>& x) const
{
  /* Initialize the vector x to zero and set up the vectors 
     r and p*/

  x.zero();

  Vector<double> r = - b; 
  Vector<double> p = r.copy();

  /* compute the initial value of A * p  */
  Vector Ap = A * p;

  /* compute the norm of b to use as a scaling factor in the
     convergence tolerance test  */
  double bNorm = b.norm2();

  /* relErr will hold the value of the relative error in the 
     current estimate of the solution.  That is, it is ||Ax - b|| /
     ||b||.  By defining it here, it will be available outside the
     loop if no convergence is obtained in maxit iterations.  */
  double relErr; 


  /* Main loop for the CG method; it will do a maximum of maxit
     iterations */
  for (int iter = 1; iter <= maxit; iter++)
    {
      double denom = p * Ap;
      double alpha = - (r * p) / denom;
      x = x + alpha * p;
      r = A * x - b;
      /* Apply relative error test  */
      relErr = r.norm2() / bNorm;
      if (relErr <= tol) 
	{
	  cout << "Convergence achieved in " << iter 
	       << " iterations" << endl;
	  return CONVERGENCE;
	}
      double beta = (r * Ap) / denom;
      p = -r + beta * p;
      Ap = A * p;
    }

  /* Loop exit because maximum number of iterations exceeded  */
  cout << "No convergence: Maximum iterations exceeded. "
       << "Relative error = " << relError << endl;
      return MAX_ITERATIONS_EXCEEDED;
}
