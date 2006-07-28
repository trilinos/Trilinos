#ifndef THYRASOLVERCG_HPP
#define THYRASOLVERCG_HPP

#include "TSFConfigDefs.hpp"
#include "TSFKrylovSolver.hpp"
#include "TSFHandleable.hpp"
#include "TSFPrintable.hpp"
#include "TSFDescribable.hpp"
#include "TSFLinearCombination.hpp"

namespace Thyra
{
  using namespace Teuchos;
  /**
   *  A simple conjugate gradient (CG) solver to illustrate features
   *  of Thyra.  It derives from IterativeSolver, Handleable,
   *  Printable, and Describable.  The virtual methods that have to be
   *  implemented are noted.
   */
  template <class Scalar>
  class ThyraSolverCG :  public IterativeSolver<Scalar>,
                         public Handleable<LinearSolverBase<Scalar> >,
                         public Printable,
                         public Describable
  {
  public:
    /** Constructor */
    ThyraSolverCG(const ParameterList& params)     
      : IterativeSolver<Scalar>(params)
    {;}


    /** Virtual destructor */
    virtual ~ThyraSolverCG() {;}

    /** \name IterativeSolver interface */
    //@{
    /**  The main solve method that takes a linear operator (A) and
	 right hand side (b) and returns a solution (x). */
    SolverState<Scalar> solve(const LinearOperator<Scalar>& A, 
			      const Vector<Scalar>& b, 
			      Vector<Scalar>& x) const ; 
 
    /** \name Printable interface */
    //@{
    /** Write to a stream.  This just prints out the description of
	the method and th parameters that are in effect at this
	call  */
    void print(ostream& os) const 
    {
      os << description() << "[" << endl;
      os << this->parameters() << endl;
      os << "]" << endl;
    }
    //@}
    
    /** \name Describable interface */
    //@{
    /** Return a brief description */
    string description() const {return "ThyraSolverCG";}
    //@}

    /** \name Handleable interface */
    //@{
    /** Return a ref count pointer to a newly created object */
    virtual RefCountPtr<LinearSolverBase<Scalar> > getRcp() 
    {return rcp(this);}
    //@}
    
  protected:


    
  };

  template <class Scalar> inline
  SolverState<Scalar> ThyraSolverCG<Scalar>
  ::solve(const LinearOperator<Scalar>& A,
                const Vector<Scalar>& b,
                Vector<Scalar>& x) const
  {
    int maxit = this->getMaxiters();
    Scalar tol = this->getTol();
    int verbosity = this->getVerbosity();
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

#endif
