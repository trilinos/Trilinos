//@HEADER
//@HEADER

#define DEBUG_LEVEL 1
#undef STORE_HESSENBERG
#define CHECK_RESIDUALS


#include "NOX_Direction_Tensor.H" // class definition
#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"
#include "stdio.h"  // for printf() -- I know, I know, this is C++...

#include "../example/NOX_Example_Vector.H"  // for printing vectors (bwb)


using namespace NOX;
using namespace NOX::Direction;

Tensor::Tensor(Parameter::List& p) 
{
  predf = NULL;
  hess = NULL;
  givensC = NULL;
  givensS = NULL;
  newHessCol = NULL;
  terrvec = NULL; 
  vecg = NULL;
  vecq = NULL;
  
  reset(p);
}


Tensor::~Tensor()
{
  delete predf;
  
  for (int i=0; i<maxDim; i++)
    delete basisVptr[i];
  if (maxDim > 0)
    delete  basisVptr;
  delete_matrix(hess);
  delete_matrix(givensC);
  delete_matrix(givensS);
  delete newHessCol;
  delete terrvec;
  delete vecg;
  delete vecq;
}


bool Tensor::reset(Parameter::List& params)
{
  // Clean up if this is not a fresh reset....
  if (basisVptr  &&  maxDim > 0) {
    cout << "Non null pointer, must free previously allocated space. \n";
    for (int i=0; i<maxDim; i++)
      delete basisVptr[i];
    delete basisVptr;
  }

  if (hess)       delete_matrix(hess);
  if (givensC)    delete_matrix(givensC);
  if (givensS)    delete_matrix(givensS);
  if (newHessCol) delete newHessCol;
  if (terrvec)    delete terrvec;
  if (vecg)       delete vecg;
  if (vecq)       delete vecq;


  paramsptr = &params;
  // if (!paramsptr->sublist("Linear Solver").isParameter("Tolerance"))
  //   paramsptr->sublist("Linear Solver").setParameter("Tolerance", 1.0e-10);

  Parameter::List& localParams = paramsptr->sublist("Linear Solver");
  tol = localParams.getParameter("Tolerance", 1.0e-4);
  reorth = localParams.getParameter("Reorthogonalize", 1);
  kmax = localParams.getParameter("kmax", 100);

  p = 3;
  allocate_vector(p, terrvec);

  maxDim = 0;   //  updated in compute() once n is known

  isFreshlyReset = true;
  
  return true;
}


bool Tensor::compute(Abstract::Vector& dir, 
		     Abstract::Group& soln, 
		     const Solver::Generic& solver)
{
  if (isFreshlyReset) {

    // Set more parameters that need info from the Group....
    n = soln.getX().length(); // bwb - Okay or better way? Always correct?

    // Manipulate some parameters further....
    if (kmax > n) 
      kmax = n;
    maxDim = kmax + p;

    /* may need to delete these pointers if reset more than once - bwb */
    basisVptr = new Abstract::Vector* [maxDim];
    for (int i=0; i<maxDim; i++)
      basisVptr[i] = soln.getX().clone(ShapeCopy);

    allocate_matrix(maxDim, kmax, hess);
    for (int i=0; i<maxDim; i++)
      for (int j=0; j<kmax; j++)
	hess[i][j] = 0;
    allocate_matrix(p+1, maxDim, givensC);
    allocate_matrix(p+1, maxDim, givensS);

    allocate_vector(maxDim, vecg);
    allocate_vector(maxDim, vecq);
    allocate_vector(maxDim, newHessCol);
    
    isFreshlyReset = false;
  }

#ifdef STORE_HESSENBERG
  double** hess2 = NULL;
  allocate_matrix(maxDim, kmax, hess2);
  for (int i=0; i<maxDim; i++)
    for (int j=0; j<kmax; j++)
      hess2[i][j] = 0;
#endif

  // Set iteration-specific local parameters....
  bool breakdown = false;
  double errTol = soln.getNormF() * tol;
  double w1 = 0;
  double w2 = 0;
  double y1 = 0;
  double qval = 0;
  Abstract::Vector* vecw = soln.getX().clone(ShapeCopy);

  // Initialize storage (this may not be needed, but it is good practice)...
  for (int i=0; i<maxDim; i++) {
    vecg[i] = 0;
    vecq[i] = 0;
    newHessCol[i] = 0;
    for (int j=0; j<p+1; j++) {
      givensC[j][i] = givensS[j][i] = 0;
    }
  }

  // Compute F at current solution....
  bool ok = soln.computeF();
  if (!ok) {
    if (Utils::doPrint(Utils::Warning))
      cout << "NOX::Direction::Tensor::compute - " 
	   << "Unable to compute F." << endl;
    return false;
  }

  // Compute Jacobian at current solution....
  ok = soln.computeJacobian();
  if (!ok) {
    if (Utils::doPrint(Utils::Warning))
      cout << "NOX::Direction::Tensor::compute - "
	   << "Unable to compute Jacobian." << endl;
    return false;
  }


  /*** Compute the tensor direction ***/

  // Compute the previous step direction, sc, and use as first basis vector...
  *basisVptr[0] = soln.getX();
  basisVptr[0]->update(1.0, solver.getPreviousSolutionGroup().getX(), -1.0);
  double normS = basisVptr[0]->norm();
  if (normS == 0) {
    ok = soln.computeNewton(paramsptr->sublist("Linear Solver"));
    dir = soln.getNewton();
    return true;
  }
  
#ifdef CHECK_RESIDUALS
  // Compute the tensor term ac....
  Abstract::Vector* acPtr = soln.getF().clone(ShapeCopy);
  Abstract::Vector* jcxscPtr = soln.getF().clone(ShapeCopy);
  soln.applyJacobian(*basisVptr[0], *jcxscPtr);
  *acPtr = *jcxscPtr;
  acPtr->update(1.0, solver.getPreviousSolutionGroup().getF(),
		-1.0, soln.getF(), -1.0);    
  acPtr->scale(1/(normS*normS*normS*normS));
#endif


#if DEBUG_LEVEL > 1
  cout << "Current Point: " << soln.getX().length() 
       << "  " << soln.getX().norm() << "  "; 
  soln.getX().print();
  cout << "Previous direction: " << basisVptr[0]->length() 
       << "  " << normS << "  "; 
  basisVptr[0]->print();
  cout << "Current ";
  soln.getF().print();
  cout << "Previous F: "; 
  solver.getPreviousSolutionGroup().getF().print();
#ifdef CHECK_RESIDUALS
  cout << "Jcxsc: " << jcxscPtr->length() << "  " 
       << jcxscPtr->norm() << "  "; 
  jcxscPtr->print();
  cout << "Tensor term: " << acPtr->length() 
       << "  " << acPtr->norm() << "  ";
  acPtr->print();
#endif
#endif

  
  // Perform partial QR factorization of n x 3 matrix...
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

#if DEBUG_LEVEL > 1
  cout << "vecg ";
  print_vector(maxDim, vecg);
#endif
  
#if DEBUG_LEVEL > 1
  cout << "prev F: "; 
  solver.getPreviousSolutionGroup().getF().print();
  cout << "Basis 0: ";
  basisVptr[0]->print();
  cout << "Basis 1: ";
  basisVptr[1]->print();
  cout << "Basis 2: ";
  basisVptr[2]->print();
  cout << "old norm: " << normS << " new: " << basisVptr[0]->norm() << "\n"; 
#endif
  
  
  // Begin iterating...
  int k = -1;
  double terr = norm(p,vecg);  // not quite right if x0 != 0  - bwb
  printf("Iteration: %2d  Tensor residual: %8e  Newton residual: %8e\n", 
	 k+1, terr, terr);
  //  cout << "Iteration: " << k+1 << "  Tensor residual error: " << terr
  //     << "  Newton residual error: " << terr << endl;
  while((terr > errTol)  &&  (k+1 < kmax)) {
    k++;
    *vecw = *basisVptr[k];
    ok = soln.applyJacobian(*vecw, *basisVptr[k+p]);  // bwb - use error check
    double normav = basisVptr[k+p]->norm();
  
    // Perform modified Gram-Schmidt....
    for (int j=0; j<k+p; j++) {
      newHessCol[j] = basisVptr[k+p]->dot(*basisVptr[j]);
      basisVptr[k+p]->update(-newHessCol[j], *basisVptr[j], 1.0);
    }
    newHessCol[k+p] = basisVptr[k+p]->norm();
    double normav2 = newHessCol[k+p];

    // Reorthogonalize, if necessary....
    if ((reorth == 3) ||
	((reorth ==1) &&
	 (normav2 == 0 || log10(normav/normav2) > 10))) {
      cout << "Reorthogonalize...\n";
      for (int j=0; j<k+p; j++) {
	double hr = basisVptr[k+p]->dot(*basisVptr[j]);
	basisVptr[k+p]->update(-hr, *basisVptr[j], 1.0);
      }
      newHessCol[k+p] = basisVptr[k+p]->norm();
    }

    // Watch out for breakdown of Arnoldi process...
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

#ifdef STORE_HESSENBERG
    for (int i=0; i<maxDim; i++)
      hess2[i][k] = newHessCol[i];
#endif


    // Calculate projected tensor term using a shortcut method...
    if (k==0) {
      double tmp = normS*normS*normS*normS;
      for (int i=0; i<p+1; i++) {
	vecq[i] = (vecq[i] - vecg[i] - newHessCol[i]*normS) / tmp;
      }
    }

    // Perform Givens Rotations to put Hessenberg in special form....
    // ... by doing previous rotations on new Hessenberg column
    if (k>0) {
      for (int i=p-1; i>=0; i--) {
	givapp(givensC[i], givensS[i], &newHessCol[i], k);
      }
    }
    // ... and then performing current rotations 
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

    // Eliminate elements in first row of Hessenberg...
    // ... by first performing rotations from previous iterations
    if (k>1) {
      for (int i=1; i<k; i++) {
	w1 = givensC[p][i] * newHessCol[i] - givensS[p][i] * newHessCol[0];
	w2 = givensS[p][i] * newHessCol[i] + givensC[p][i] * newHessCol[0];
	newHessCol[i] = w1;
	newHessCol[0] = w2;
      }
    }
    // ... and then by performing rotations from the current iteration
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
	hess[0][0] = w2;

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
	
    // Put the working vector into the Hessenberg matrix...
    for (int i=0; i<maxDim ; i++)
      hess[i][k] = newHessCol[i];
    
    // Find smallest magnitude root of the quadratic equation (first row)...
    double qa = vecq[0] * normS * normS;
    double qb = hess[0][0];
    double qc = vecg[0];
    double discriminant = qb*qb - 4*qa*qc;
    if (discriminant < 0) {
      y1 = -qb/qa/2;
      qval = (qa*y1*y1 + qb*y1 + qc) * hess[0][0];
#if DEBUG_LEVEL > 0
      cout << "Discriminant is negative! \n";
#endif
    }
    else {
      qval = 0;
      if (abs(qa/qb) < 1e-8) {
#if DEBUG_LEVEL > 0
	cout << "qa is relatively small\n";
#endif 
	y1 = -qc/qb;
      }
      else {
	double tmp1 = (-qb + sqrt(discriminant)) / (2*qa);
	double tmp2 = (-qb - sqrt(discriminant)) / (2*qa);
	y1 = (abs(tmp1) < abs(tmp2)) ? tmp1 : tmp2;
      }
    }	

#if DEBUG_LEVEL > 1
    cout << "normS and y1: " << normS << "  " << y1 << endl;
    cout << "NormS * y1 = " << normS*y1 << endl;
    print_vector(maxDim, vecg);
    print_vector(maxDim, vecq);
#endif

    // Update the tensor residual norm...
    double beta2 = normS * normS * y1 * y1;
    for (int i=0; i<p; i++) {
      terrvec[i] = -vecg[k+i+1] - vecq[k+i+1] * beta2;
    }
    terr = norm(p,terrvec);

    // Update the Newton residual norm...
    for (int i=0; i<p; i++)
      terrvec[i] = -vecg[k+i+1];
    double nerr = norm(p,terrvec);

    // Print residual errors...
    printf("Iteration: %2d  Tensor residual: %8e  Newton residual: %8e\n", 
	   k+1, terr, nerr);
    //cout << "Iteration: " << k+1 << "  Tensor residual error: " << terr << 
    //  "  Newton residual error: " << nerr << endl;
  }  // end while loop


#if DEBUG_LEVEL > 1
#ifdef STORE_HESSENBERG
  cout << "Original Hessenberg: \n";
  print_matrix(maxDim,kmax,hess2);
#endif
  cout << "\n\nModified Hessenberg: \n";
  print_matrix(maxDim,kmax,hess);
  cout << "modified vecg ";
  print_vector(maxDim, vecg);
#endif
  
  int iterations = k+1;
  int* pindex = NULL;
  
  // Solve linear system for Newton direction...
  /* ...but we can't update Newton vector in group... */
  allocate_vector(iterations, pindex);
  //pindex = allocate_vector_int(iterations);
  pindex[iterations-1] = 0;
  for (int i=0; i<iterations-1; i++)
    pindex[i] = i+1;
  double* yn = backsolve(hess, vecg, pindex, iterations);
#if DEBUG_LEVEL > 1
  cout << "yn ";
  print_vector(iterations, yn);
#endif
  dir.init(0);
  for (int i=0; i<iterations; i++)
    dir.update(-yn[i], *basisVptr[i], 1);
#if DEBUG_LEVEL > 1
  cout << "Newton Direction 1: ";
  dir.print();
#endif


#if DEBUG_LEVEL > 1
  if (predf == NULL) {
    predf = soln.getF().clone(ShapeCopy);
  }
  soln.applyJacobian(dir, *predf);    
  predf->update(1.0, soln.getF(), 1.0);
  double residual = predf->norm();
  printf("Actual Newton1 residual = %8e\n", residual);
#endif

  
  // Compute the Newton direction directly (temporary code)
  ok = soln.computeNewton(paramsptr->sublist("Linear Solver"));

#if DEBUG_LEVEL > 1
  cout << "Newton Direction 2: ";
  soln.getNewton().print();
  soln.applyJacobian(soln.getNewton(), *predf);    
  predf->update(1.0, soln.getF(), 1.0);
  double residual2 = predf->norm();
  printf("Actual Newton2 residual = %8e\n", residual2);
#endif

  // Computing the Newton step failed, but maybe the step is okay to use...
  if (!ok) {

    if (predf == NULL) {
      predf = soln.getF().clone(ShapeCopy);
    }

    soln.applyJacobian(soln.getNewton(), *predf);    
    predf->update(1.0, soln.getF(), 1.0);
    double accuracy = predf->norm();
    if (accuracy < soln.getNormF() ) {
      ok = true;
      if (Utils::doPrint(Utils::Warning)) 
	cout << "WARNING: NOX::Direction::Tensor::compute - "
	     << "Tensor solve failure.\n" 
	     << "Desired accuracy is " 
	     << Utils::sci(paramsptr->sublist("Linear Solver").getParameter("Tolerance", 1.0e-10)) << ".\n"
	     << "Using solution with accuracy of " << Utils::sci(accuracy) 
	     << "." << endl;
    }
  }

  if (!ok) {
    if (Utils::doPrint(Utils::Warning))
      cout << "NOX::Direction::Tensor::compute - "
	   << "Unable to compute Tensor direction." << endl;
    return false;
  }


  // Solve system for tensor direction...
  double beta2 = normS * normS * y1 * y1;
  double* vecz = NULL;
  allocate_vector(maxDim, vecz);
  for (int i=0; i<maxDim; i++) {
    vecz[i] = -vecg[i] - vecq[i] * beta2 - hess[i][0] * y1;
  }
  double* yt = backsolve(hess, vecz, pindex, iterations-1, iterations);
  yt[0] = y1;
  dir.init(0);
  for (int i=0; i<iterations; i++)
    dir.update(yt[i], *basisVptr[i], 1);

#if DEBUG_LEVEL > 1
  cout << "Tensor Direction: "; 
  dir.print();
#endif
#if DEBUG_LEVEL > 0
  double beta = basisVptr[0]->dot(dir) * normS;
  printf("Beta %e, ||sc||*y1 %e\n", beta, normS*y1);
#ifdef CHECK_RESIDUALS
  if (predf == NULL) {
    predf = soln.getF().clone(ShapeCopy);
  }
  soln.applyJacobian(dir, *predf);    
  predf->update(1.0, soln.getF(), beta*beta, *acPtr, 1.0);
  double residual3 = predf->norm();
  printf("Actual tensor residual = %8e\n", residual3);
#endif
#endif

#if DEBUG_LEVEL > 1
  print_matrix(maxDim, kmax, hess);
  print_vector(maxDim, vecg);
  print_vector(maxDim, vecq);
  printf("normS = %e, y1 = %e\n", normS, y1);
  print_vector(maxDim, vecz);
  print_vector(iterations, pindex);
  print_vector(iterations, yt);
#endif
  
  
  // Clean up
  delete vecw;   // maybe get rid of this one?
  delete pindex;
  delete yn;
  delete yt;
  delete vecz;  // maybe try overwriting some other vector - bwb
#ifdef CHECK_RESIDUALS
  delete acPtr;
  delete jcxscPtr;
#endif
#ifdef STORE_HESSENBERG
  delete_matrix(hess2);
#endif

  // Set search direction...
  // dir = soln.getNewton();
  // dir.scale(-1.0);   //bwb - to test linesearch

  return ok;
}



// private

void** Tensor::allocate_matrix(int m, int n, double**& a)
{
  if (a) {
    // delete_matrix(a);
    cout << "Warning: Possibly a previously allocated matrix\n";
  }

  // allocate memory for storing a rectangular m x n matrix
  // a = (double **) malloc( (size_t)(m*sizeof(double)) );
  a = new double* [m];
  if (!a) printf("Allocation error in allocate_matrix()\n");
  // a[0] = (double *) malloc( (size_t)(m*n*sizeof(double)) );
  a[0] = new double [m*n];
  if (!a) printf("Allocation error in allocate_matrix()\n");

  for (int i=1; i<m; i++) 
    a[i] = a[i-1] + n;

  return (void**) a;
}


void* Tensor::allocate_vector(int n, int*& x)
{
  if (x) {
    // delete x;
    cout << "Warning: Possibly a previously allocated vector\n";
  }

  x = new int [n];
  if (!x) printf("Allocation error with malloc()\n");
  return (void*) x;
}


void* Tensor::allocate_vector(int n, double*& x)
{
  if (x) {
    // delete x;
    cout << "Warning: Possibly a previously allocated vector\n";
  }

  x = new double [n];
  if (!x) printf("Memory allocation error in allocate_vector().\n");
  return (void*) x;
}


void Tensor::delete_matrix(double** A)
   /* Free the memory previously allocated in allocate_matrix()  */
{
  delete A[0];
  delete A;
}


void Tensor::print_matrix(int m, int n, double** A)
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


void Tensor::print_vector(int n, double* x)
{
  int i;

  printf("Vector = \n");
  for (i=0; i<n; i++)
    printf("   %16.8e\n", x[i]);
  printf("\n");
}


void Tensor::print_vector(int n, int* x)
{
  int i;

  printf("Vector = \n");
  for (i=0; i<n; i++)
    printf("   %d\n", x[i]);
  printf("\n");
}


double Tensor::inner_product(int n, double* x, double* y)
{
  int i;
  double sum = 0.0;

  for (i=0; i<n; i++)
    {
      sum += x[i]*y[i];
    }

  return sum;
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


double* Tensor::backsolve(double** U, double* b, int* perm, int n, int dim=0)
     /* This function solves the triangular system Ux=b when provided an 
      * upper triangular matrix U. The array perm is a permutation array.
      * The pointer returned is the newly created solution vector x.
      * The variable dim is the maximum dimension of the solution vector, 
      * which may be different due to the indices in the permutation array.
      */
{
  double* x = NULL;
  double temp;
  int i, j;
  int pi, pj;

  dim = (dim>n) ? dim : n;

  //  Allocate memory for solution vector
  allocate_vector(dim, x);

  //  Solve Ux=y with backward substitution
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
 
