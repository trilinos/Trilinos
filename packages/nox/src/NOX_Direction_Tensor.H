//@HEADER
//@HEADER

#ifndef NOX_DIRECTION_TENSOR_H
#define NOX_DIRECTION_TENSOR_H

#include "NOX_Direction_Generic.H" // base class
#include "NOX_Parameter_List.H"	   // class data element 

namespace NOX {
namespace Direction {

//! %Tensor direction computation
/*! 
  Computes the Tensor direction by finding the root or smallest magnitude 
  minimizer of the local model
  \f[ M_T(x_c+d) = F_c + Jd + a_c(s^Td)^2. \f]
  
  To use this direction, specify that the "Method" is "Tensor" in the
  "Direction" sublist of the parameters that are passed to the solver
  (see NOX::Direction::Manager for more information on choosing the
  search direction).
*/  

class Tensor : public Generic {

public:

  //! Constructor
  Tensor(Parameter::List& params);

  //! Destructor
  virtual ~Tensor();

  // derived
  virtual bool reset(Parameter::List& params);

  // derived
  virtual bool compute(Abstract::Vector& dir, Abstract::Group& grp, 
			  const Solver::Generic& solver);


protected:

private:

  // Various functions that are taken from C to at least get a tensor
  // method implemented.  At some point these functions (as well as the
  // non-abstract vector/matrix data structures) should be replaced
  // with a class or abstraction layer.  I realize this is not elegant, but
  // it is quick and dirty and will work okay for the moment.
  void** allocate_matrix(int m, int n, double**& a);
  void* allocate_vector(int n, int*& x);
  void* allocate_vector(int n, double*& x);
  void delete_matrix(double** A);
  void print_matrix(int m, int n, double** A);
  void print_vector(int n, double* x);
  void print_vector(int n, int* x);
  double inner_product(int n, double* x, double* y);
  double norm(int n, double* x);
  void givapp(double* c, double* s, double* v, int k);
  double* backsolve(double** U, double* b, int* p, int n, int dim=0);
  
private:

  //! Types of preconditioning 
  enum preconditioningType {None, Left, Right};

  //! Types of preconditioning 
  enum reorthogonalizeType {Never, AsNeeded, Always};

  //! Flag for type of preconditioning
  preconditioningType precondition;

  //! Flag for reorthogonalizing basis after first pass of Gram-Schmidt
  reorthogonalizeType reorth;

  //! "Direction" sublist with parameters for the direction vector
  Parameter::List* paramsPtr;

  //! Vector containing the predicted RHS
  Abstract::Vector* predf;

  //! Array of orthonormal basis vectors
  Abstract::Vector** basisVptr;

  //! Flag for "continuing" with reset() inside compute() where n is known
  bool isFreshlyReset;

  //! Problem size
  int n;

  //! Maximum number of Krylov iterations (kmax <= n)
  int kmax;

  //! Number of starting vectors in initial basis
  int p;

  //! Maximum dimension of data structures (maxDim = kmax + p)
  int maxDim;

  //! Threshold tolerance for tensor solve
  double tol;

  //! Band Hessenberg matrix from Arnoldi process
  double** hess;

  //! Matrix of coefficients for Givens rotations
  double** givensC;
  double** givensS;

  //! Working vector for the new column of Hessenberg matrix
  double* newHessCol;

  //! Projected right hand sides
  double* vecg;
  double* vecq;

  //! Vector for calculating the residual error
  double* terrvec; 
    
};
} // namespace Direction
} // namespace NOX
#endif
