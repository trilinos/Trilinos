//@HEADER
//@HEADER

#define DEBUG_LEVEL 2

#include "NOX_Direction_Tensor.H" // class definition
#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"
#include "stdio.h"  // for printf() -- I know, I know, this is C++...

#include "../example/NOX_Example_Vector.H"  // bwb


using namespace NOX;
using namespace NOX::Direction;

Tensor::Tensor(Parameter::List& p) 
{
  predf = NULL;
  stepdir = NULL;
  reset(p);
}

Tensor::~Tensor()
{
  cout << "\n\nDestroying Tensor object\n";  // bwb
  delete predf;
  
  cout << "1/4 Finished with Tensor object\n";
  for (int i=0; i<maxDim; i++)
    delete basisVptr[i];
  if (maxDim > 0)
    delete  basisVptr;
  cout << "1/2 Finished with Tensor object\n";
  free_matrix(hess);
  free_matrix(givensC);
  free_matrix(givensS);
  cout << "3/4 Finished with Tensor object\n";
  free(newHessCol);
  free(terrvec);
  free(x);
  free(vecg);
  free(vecq);

  cout << "Finished with Tensor object\n";
}

bool Tensor::reset(Parameter::List& params)
{
  paramsptr = &params;
  if (!paramsptr->sublist("Linear Solver").isParameter("Tolerance"))
    paramsptr->sublist("Linear Solver").setParameter("Tolerance", 1.0e-10);

  Parameter::List& localParams = paramsptr->sublist("Linear Solver");
  tol = localParams.getParameter("Tolerance", 1.0e-6);
  reorth = localParams.getParameter("Reorthogonalize", 1);
  kmax = localParams.getParameter("kmax", 20);
  p = 3;

  if (basisVptr  &&  maxDim > 0) {
    "Non null pointer, must free previously allocated space. \n";
    delete basisVptr;  // bwb - need to fix this
  }
  maxDim = 0;   //  updated in compute() once n is known

  terrvec = allocate_vector_double(p);
  
  isFreshlyReset = true;
  
  return true;
}

bool Tensor::compute(Abstract::Vector& dir, 
		     Abstract::Group& soln, 
		     const Solver::Generic& solver)
{
  if (isFreshlyReset) {
    // Set more parameters
    n = soln.getX().length(); // bwb - okay or better way? Always correct?

    // Manipulate some parameters further
    if (kmax > n) 
      kmax = n;
    maxDim = kmax + p;

    /* may need to delete these if reset more than once - bwb */
    basisVptr = new Abstract::Vector* [maxDim];
    for (int i=0; i<maxDim; i++)
      basisVptr[i] = soln.getX().clone(ShapeCopy);

    hess = allocate_matrix(maxDim,kmax);
    print_matrix(maxDim, kmax, hess);
    for (int i=0; i<maxDim; i++)
      for (int j=0; j<kmax; j++)
	hess[i][j] = 0;
    x = allocate_vector_double(kmax);
    givensC = allocate_matrix(p+1, maxDim);
    givensS = allocate_matrix(p+1, maxDim);

    vecg = allocate_vector_double(maxDim);
    vecq = allocate_vector_double(maxDim);
    newHessCol = allocate_vector_double(maxDim);
    
    isFreshlyReset = false;
  }

  double** hess2 = allocate_matrix(maxDim,kmax);
  for (int i=0; i<maxDim; i++)
    for (int j=0; j<kmax; j++)
      hess2[i][j] = 0;

  // Iteration-specific local parameters
  double errTol = soln.getNormF() * tol;
  bool breakdown = false;
  Abstract::Vector* vecw = soln.getX().clone(ShapeCopy);
  double w1;
  double w2;
  double y1;
  double qval;
  
  cout << "Local tolerance = " << tol << "  " << errTol << "  " << n;
  cout << "  " << kmax << endl;

  // Allocate storage for vectors and matrices
  //  Abstract::Vector& basisV = *basisVptr[0];
  //  cout << "Basis vector: " << basisV.length() <<
  //    dynamic_cast<const NOX::Example::Vector&>(basisV) << endl;


  // Clear storage
  for (int i=0; i<maxDim; i++) {
    vecg[i] = 0;
    vecq[i] = 0;
  }


  // Compute F at current solution
  bool ok = soln.computeF();
  if (!ok) {
    if (Utils::doPrint(Utils::Warning))
      cout << "NOX::Direction::Tensor::compute - " <<
	"Unable to compute F." << endl;
    return false;
  }

  // Compute Jacobian at current solution.
  ok = soln.computeJacobian();
  if (!ok) {
    if (Utils::doPrint(Utils::Warning))
      cout << "NOX::Direction::Tensor::compute - " <<
	"Unable to compute Jacobian." << endl;
    return false;
  }

  cout << "Current Point: " << soln.getX().length() <<
    "  " << soln.getX().norm() << "  " <<
    dynamic_cast<const NOX::Example::Vector&>(soln.getX()) << endl;

  
  // Compute the tensor direction


  // Compute the previous step direction, sc
  *basisVptr[0] = soln.getX();
  basisVptr[0]->update(1.0, solver.getPreviousSolutionGroup().getX(), -1.0);
  double normS = basisVptr[0]->norm();
  cout << "Previous direction: " << basisVptr[0]->length() <<
    "  " << normS << "  " <<
    dynamic_cast<const NOX::Example::Vector&>(*basisVptr[0]) << endl;
  if (normS == 0) {
    ok = soln.computeNewton(paramsptr->sublist("Linear Solver"));
    dir = soln.getNewton();
    return true;
  }
  
  // Compute the tensor term, ac
  Abstract::Vector* acPtr = soln.getF().clone(ShapeCopy);
  Abstract::Vector* jcxscPtr = soln.getF().clone(ShapeCopy);
  soln.applyJacobian(*basisVptr[0], *jcxscPtr);
  *acPtr = *jcxscPtr;
  acPtr->update(1.0, solver.getPreviousSolutionGroup().getF(),
		-1.0, soln.getF(), -1.0);    
  cout << "Tensor term: " << acPtr->length() <<
    "  " << acPtr->norm() << "  " <<
    dynamic_cast<const NOX::Example::Vector&>(*acPtr) << endl;

  
  // Perform partial QR factorization on n x 3 matrix
  basisVptr[0]->scale(1.0/normS);

  *basisVptr[1] = solver.getPreviousSolutionGroup().getF();
  vecq[0] = basisVptr[1]->dot(*basisVptr[0]);
  basisVptr[1]->update(-vecq[0], *basisVptr[0], 1.0);
  vecq[1] = basisVptr[1]->norm();
  basisVptr[1]->scale(1.0/vecq[1]);

  *basisVptr[2] = soln.getF();
  vecg[0] = basisVptr[2]->dot(*basisVptr[0]);
  basisVptr[2]->update(-vecg[0], *basisVptr[0], 1.0);
  vecg[1] = basisVptr[2]->dot(*basisVptr[1]);
  basisVptr[2]->update(-vecg[1], *basisVptr[1], 1.0);
  vecg[2] = basisVptr[2]->norm();
  basisVptr[2]->scale(1.0/vecg[2]);

#if DEBUG_LEVEL > 0
  cout << "vecg ";
  print_vector(maxDim, vecg);
#endif
  
#if DEBUG_LEVEL > 1
  cout << "prev F: " <<
    dynamic_cast<const NOX::Example::Vector&>(solver.getPreviousSolutionGroup().getF()) << endl;
  cout << "Basis 0: " <<
    dynamic_cast<const NOX::Example::Vector&>(*basisVptr[0]) << endl;
  cout << "Basis 1: " <<
    dynamic_cast<const NOX::Example::Vector&>(*basisVptr[1]) << endl;
  cout << "Basis 2: " <<
    dynamic_cast<const NOX::Example::Vector&>(*basisVptr[2]) << endl;
  cout << "old norm: " << normS << " new: " << basisVptr[0]->norm() << "\n"; 
#endif
  
  
  // Begin iterating
  int k=-1;
  double terr = norm(3,vecg);
  cout << "terr = " << terr << endl;
  while((terr > errTol)  &&  (k+1 < kmax)) {
    k = k+1;
    *vecw = *basisVptr[k];
    ok = soln.applyJacobian(*vecw, *basisVptr[k+p]);  // bwb - use error check
    double normav = basisVptr[k+p]->norm();
  
    // Modified Gram-Schmidt
    for (int j=0; j<k+p; j++) {
      newHessCol[j] = basisVptr[k+p]->dot(*basisVptr[j]);
      basisVptr[k+p]->update(-newHessCol[j], *basisVptr[j], 1.0);
    }
    newHessCol[k+p] = basisVptr[k+p]->norm();
    double normav2 = newHessCol[k+p];

    // Reorthogonalize, if necessary
    if ((reorth == 3) || ((reorth ==1) &&
			  (normav + 0.001 * normav2 == normav))) {
      cout << "Reorthogonalize...\n";
      for (int j=0; j<k+p; j++) {
	double hr = basisVptr[k+p]->dot(*basisVptr[j]);
	basisVptr[k+p]->update(-hr, *basisVptr[j], 1.0);
      }
      newHessCol[k+p] = basisVptr[k+p]->norm();
    }

    // Watch out for happy breakdown
    if (newHessCol[k+p] != 0)
      basisVptr[k+p]->scale(1.0/newHessCol[k+p]);
    else {
      cout << "Breakdown\n";
      breakdown = true;
    }

#if DEBUG_LEVEL > 1
    cout << "Iteration " << k << ":\n";
    print_vector(k+p+1, newHessCol);
#endif 

    for (int i=0; i<k+p+1; i++)
      hess2[i][k] = newHessCol[i];

    // Givens Rotations
    if (k>0) {
      for (int i=p-1; i>=0; i--) {
	givapp(givensC[i], givensS[i], &newHessCol[i], k);
      }
    }
    

    for (int j=k+p-1; j>=k; j--) {
      double nu = sqrt(newHessCol[j]*newHessCol[j] +
		       newHessCol[j+1]*newHessCol[j+1]);
      if (nu != 0) {
	double c1 = newHessCol[j] / nu;
	double s1 = -newHessCol[j+1] / nu;
	newHessCol[j] = c1 * newHessCol[j] - s1 * newHessCol[j+1];
	newHessCol[j+1] = 0;
	givapp(&c1, &s1, &vecg[j], 1);
	givapp(&c1, &s1, &vecq[j], 1);
	givensC[j-k][k] = c1;
	givensS[j-k][k] = s1;
      }
    }

    // Eliminate elements in first row of Hessenberg
    if (k>1) {
      for (int i=1; i<k-1; i++) {
	w1 = givensC[p][i] * newHessCol[i] - givensS[p][i] * newHessCol[0];
	w2 = givensS[p][i] * newHessCol[i] + givensC[p][i] * newHessCol[0];
	newHessCol[i] = w1;
	newHessCol[0] = w2;  // should be 0
      }
    }
    if (k>0) {
      double nu = sqrt(newHessCol[k]*newHessCol[k] +
		       newHessCol[0]*newHessCol[0]);
      if (nu != 0) {
	double c1 = newHessCol[k] / nu;
	double s1 = -newHessCol[0] / nu;
	newHessCol[k] = c1 * newHessCol[k] - s1 * newHessCol[0];
	newHessCol[0] = 0;
	
	w1 = c1 * hess[k][0] - s1 * hess[0][0];
	w2 = s1 * hess[k][0] + c1 * hess[0][0];
	hess[k][0] = w1;
	hess[0][0] = w2;  // should be 0

	w1 = c1 * vecg[k] - s1 * vecg[0];
	w2 = s1 * vecg[k] + c1 * vecg[0];
	vecg[k] = w1;
	vecg[0] = w2;

	w1 = c1 * vecq[k] - s1 * vecq[0];
	w2 = s1 * vecq[k] + c1 * vecq[0];
	vecq[k] = w1;
	vecq[0] = w2;

	givensC[p][k] = c1;
	givensS[p][k] = s1;
      }
    }
	
    // Put everything back in place in the Hessenberg matrix
    for (int i=0; i<=k+p; i++)
      hess[i][k] = newHessCol[i];
    
    // Find root of q(beta) equation
    double qa = vecq[0] * normS;
    double qb = hess[0][0];
    double qc = vecg[0];
    double discriminant = qb*qb - 4*qa*qc;
    if (discriminant < 0) {
      y1 = -qb/qa/2;
      qval = (qa*y1*y1 + qb*y1 + qc) * hess[0][0];
      cout << "Discriminant is negative! \n";
    }
    else {
      qval = 0;
      if (abs(qa/qb) < 1e-8)
	y1 = -qc/qb;
      else {
	double tmp1 = (-qb + sqrt(discriminant)) / (2*qa);
	double tmp2 = (-qb - sqrt(discriminant)) / (2*qa);
	y1 = (abs(tmp1) < abs(tmp2)) ? tmp1 : tmp2;
      }
    }	

#if DEBUG_LEVEL > 1
    cout << "normS and y1: " << normS << "  " << y1 << endl;
    print_vector(maxDim, vecg);
    print_vector(maxDim, vecq);
#endif
    
    // Update the residual norm
    for (int i=0; i<p; i++)
      terrvec[i] = -vecg[k+i+1] - vecq[k+i+1] * (normS*normS * y1*y1);
    terr = norm(p,terrvec);

    // print_vector(p,terrvec);
    
    cout << "Iteration: " << k << "  residual error: " << terr << endl;
  }  // end while loop


#if DEBUG_LEVEL > 0  
  cout << "Original Hessenberg: \n";
  print_matrix(maxDim,kmax,hess2);
  cout << "\n\nModified Hessenberg: \n";
  print_matrix(maxDim,kmax,hess);
  cout << "modified vecg ";
  print_vector(maxDim, vecg);
#endif
  
  int iterations = k+1;
  
  // Solve linear system for Newton direction
  int* pindex;
  pindex = allocate_vector_int(iterations);
  pindex[iterations-1] = 0;
  for (int i=0; i<iterations-1; i++)
    pindex[i] = i+1;
  double* yn = backsolve(hess, vecg, pindex, iterations);

#if DEBUG_LEVEL > 0 
  cout << "yn ";
  print_vector(iterations, yn);
#endif

#ifdef ACTUAL_VALUES
  yn[0] =   0.21756434933327;
  yn[1] =  -0.00166296904318;
  yn[2] =   0.01756165993543;
  yn[3] =   0.00335222653519;
  yn[4] =   0.00102653837575;
  yn[5] =   0.00083610648822;
#endif
  
  dir.init(0);
  for (int i=0; i<iterations; i++)
    dir.update(-yn[i], *basisVptr[i], 1);
  cout << "Newton Direction 1: " <<
    dynamic_cast<const NOX::Example::Vector&>(dir) << endl;
  /* ...but we can't update Newton vector in group... */
  
  // Compute the Newton direction directly
  ok = soln.computeNewton(paramsptr->sublist("Linear Solver"));
  cout << "Newton Direction 2: " <<
    dynamic_cast<const NOX::Example::Vector&>(soln.getNewton()) << endl;

  // It didn't work, but maybe it's ok anyway...
  if (!ok) {

    if (predf == NULL) {
      predf = soln.getF().clone(ShapeCopy);
    }

    soln.applyJacobian(soln.getNewton(), *predf);    
    predf->update(-1.0, soln.getF(), 1.0);
    double accuracy = predf->norm();
    if (accuracy < 1) {
      ok = true;
      if (Utils::doPrint(Utils::Warning)) 
	cout << "WARNING: NOX::Direction::Tensor::compute - Tensor solve failure.\n" 
	     << "Desired accuracy is " 
	     << Utils::sci(paramsptr->sublist("Linear Solver").getParameter("Tolerance", 1.0e-10)) << ".\n"
	     << "Using solution with accuracy of " << Utils::sci(accuracy) << "." << endl;
    }
  }

  if (!ok) {
    if (Utils::doPrint(Utils::Warning))
      cout << "NOX::Direction::Tensor::compute - Unable to compute Tensor direction." << endl;
    return false;
  }


  // Solve system for tensor direction
  double* vecz = allocate_vector_double(maxDim);
  for (int i=0; i<maxDim; i++) {
    vecz[i] = -vecg[i] - vecq[i] * normS*normS*y1*y1 - hess[i][0]*y1;
  }
  double* yt = backsolve(hess, vecz, &pindex[1], iterations-1);

#if DEBUG_LEVEL > 1
  cout << "yt ";
  print_vector(iterations, yt);
#endif
  
  dir.update(y1,*basisVptr[0]);
  for (int i=0; i<iterations-1; i++)
    dir.update(yt[i], *basisVptr[i], 1);
  cout << "Tensor Direction: " <<
    dynamic_cast<const NOX::Example::Vector&>(dir) << endl;
  
  
  // Clean up
  delete vecw;   // maybe get rid of this one?
  delete acPtr;
  delete jcxscPtr;
  free(pindex);
  free(yn);
  free(yt);
  free(vecz);  // maybe try overwriting some other vector - bwb

  free_matrix(hess2);

  // Set search direction.
  // dir = soln.getNewton();
  // dir.scale(-1.0);   //bwb

  return ok;
}



// private

double** Tensor::allocate_matrix(int m, int n)
{
  double **a;
  int i;

  /* allocate memory for storing a rectangular m x n matrix */
  a = (double **) malloc( (size_t)(m*sizeof(double)) );
  if (!a) printf("Allocation error with malloc()\n");
  a[0] = (double *) malloc( (size_t)(m*n*sizeof(double)) );
  if (!a) printf("Allocation error with malloc()\n");

  for (i=1; i<m; i++) a[i] = a[i-1] + n;

  return a;
}

double** Tensor::copy_matrix(int m, int n, double **A)
{
  int i,j;
  double **C;

  C = allocate_matrix(m,n);
  for(i=0; i<m; i++)
    for (j=0; j<n; j++)
      C[i][j] = A[i][j];
  return C;
}


void Tensor::free_matrix(double **A)
   /* Free the memory previously allocated in initialize_matrix()  */
{
  free((double *) (A[0]));
  free((double **) (A));
}


void Tensor::print_matrix(int m, int n, double **A)
   /*  Print a matrix in conventional form.  */
{
  int i,j;

  for (i=0; i<m; i++) 
    {
      for (j=0; j<n;j++)
          printf("%16.8e ", A[i][j]);
      printf("\n");
    }
}


double* Tensor::allocate_vector_double(int n)
{
  double *x;

  x = (double *) malloc( (size_t)(n*sizeof(double)) );
  if (!x) printf("Allocation error with malloc()\n");
  return x;
}


int* Tensor::allocate_vector_int(int n)
{
  int *x;

  x = (int *) malloc( (size_t)(n*sizeof(int)) );
  if (!x) printf("Allocation error with malloc()\n");
  return x;
}


void Tensor::print_vector(int n, double *x)
{
  int i;

  printf("Vector = \n");
  for (i=0; i<n; i++)
    printf("   %16.8e\n", x[i]);
  printf("\n");
}

void Tensor::print_vector(int n, int *x)
{
  int i;

  printf("Vector = \n");
  for (i=0; i<n; i++)
    printf("   %d\n", x[i]);
  printf("\n");
}


double Tensor::inner_product(int n, double *x, double *y)
{
  int i;
  double sum = 0.0;

  for (i=0; i<n; i++)
    {
      sum += x[i]*y[i];
    }

  return sum;
}


double* Tensor::multiply(int m, int n, double** A, double* x)
{
  int i,j;
  double* y;
  
  y = allocate_vector_double(m);
  
  for (i=0; i<m; i++) {
    y[i] = 0;
    for (j=0; j<n; j++) 
      y[i] += A[i][j]*x[j];
  }
  return y;
}


double Tensor::norm(int n, double* x)
{
  double temp = inner_product(n,x,x);
  return sqrt(temp);
}


void Tensor::givapp(double* c, double* s, double* v, int k)
{
  double w1;
  double w2;
  
  for (int i=0; i<k; i++) {
    w1 = c[i] * v[i] - s[i] * v[i+1];
    w2 = s[i] * v[i] + c[i] * v[i+1];
    v[i]   = w1;
    v[i+1] = w2;
  }

}


double* Tensor::backsolve(double **U, double *b, int *perm, int n)
     /* This function solves the system Ux=b when provided an upper
      * triangular matrix U. The array perm is a permutation array.
      * The pointer returned is the solution vector x.
      */
{
  double *x;
  double temp;
  int i,j;
  int pi, pj;

  /*  Allocate memory for solution vector */
  x = allocate_vector_double(n);

  /*  Solve Ux=y with backward substitution  */
  for (i=n-1; i>=0; i--) {
    pi = perm[i];
    temp = b[pi];
    for (j=i+1; j<n; j++) {
      pj = perm[j];
      temp -= U[pi][pj] * x[pj];
    }
    x[pi] = temp / U[pi][pi];
  }  

  return x;
}


  /*
  // Test matrix stuff
  x[0] = 1;
  x[1] = 2;
  for (int i=0; i<maxDim; i++)
    for (int j=0; j<kmax; j++)
      hess[i][j] = i*kmax + j;
  cout << "Hessenberg: \n";
  print_matrix(maxDim,kmax,hess);
  double* b = multiply(maxDim,kmax,hess,x);
  print_vector(kmax,x);
  print_vector(maxDim,b);
  free(b);
  */
  
