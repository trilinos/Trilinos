//@HEADER
//@HEADER

#define DEBUG_LEVEL 0
#undef STORE_HESSENBERG
#undef CHECK_ORTHOGONALITY
#define CHECK_RESIDUALS
#undef USE_NEWTON_DIR
#define BUILD_UP_BASIS

#include "NOX_Direction_Tensor.H" // class definition
#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"
#include "stdio.h"  // for printf() -- I know, I know, this is C++...


using namespace NOX;
using namespace NOX::Direction;

Tensor::Tensor(Parameter::List& params) 
{
  predf = NULL;
  vecw = NULL;
  dTensor = NULL;
  dTLambda = NULL;
  dNewton = NULL;
  hess = NULL;
  givensC = NULL;
  givensS = NULL;
  newHessCol = NULL;
  terrvec = NULL; 
  vecg = NULL;
  vecq = NULL;
  vecz = NULL;
  
  reset(params);
}


Tensor::~Tensor()
{
  delete predf;
  
  for (int i=0; i<maxDim; i++)
    delete basisVptr[i];
  if (maxDim > 0)
    delete [] basisVptr;
  delete vecw;
  delete dTensor;
  delete dTLambda;
  delete dNewton;

  delete_matrix(hess);
  delete_matrix(givensC);
  delete_matrix(givensS);
  delete [] newHessCol;
  delete [] terrvec;
  delete [] vecg;
  delete [] vecq;
  delete [] vecz;
}


bool Tensor::reset(Parameter::List& params)
{
  /* 
   *  This code may not be entirely correct for repeatedly resetting
   *  the method.  Memory leaks might be an issue here.
   */

  /*
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
    if (vecq) delete vecq; 
  */
  
  paramsPtr = &params;
    
  Parameter::List& localParams = paramsPtr->sublist("Linear Solver");
  tol = localParams.getParameter("Tolerance", 1.0e-4);
  kmax = localParams.getParameter("Max Iterations", 100);
  outputFreq = localParams.getParameter("Output Frequency", 20);

  string choice = localParams.getParameter("Reorthogonalize", "As Needed");
  if (choice == "Never") {
    reorth = Never;
  }
  else if (choice == "As Needed") {
    reorth = AsNeeded;
  }
  else if (choice == "Always") {
    reorth = Always;
  }
  else {
    cout << "ERROR: NOX::Direction::Tensor::reset() - the choice of "
	 << "\"Reorthogonalization\" parameter is invalid." << endl;
    throw "NOX Error";
  }


  choice = localParams.getParameter("Preconditioning", "None");
  if (choice == "None") {
    precondition = None;
    localParams.setParameter("PreconditioningSide", "None");
  }
  else {
    choice = localParams.getParameter("PreconditioningSide", "Left");
    if (choice == "Left") {
      precondition = Left;
    }
    else if (choice == "Right") {
      precondition = Right;
    }
    else if (choice == "None") {
      precondition = None;
    }
    else {
      cout << "ERROR: NOX::Direction::Tensor::reset() - the choice of "
	   << "\"PreconditioningSide\" parameter is invalid." << endl;
      throw "NOX Error";
    }
  }

  maxp = 3;
  allocate_vector(maxp, terrvec);

  maxDim = 0;   //  updated in compute() once n is known

  arnoldiIters = 0;

  isFreshlyReset = true;
  isTensorCalculated = false;
  isCLParamsCalculated = false;

  return true;
}


bool Tensor::compute(Abstract::Vector& dir, 
		     Abstract::Group& soln, 
		     const Solver::Generic& solver)
{
  if (isFreshlyReset) {

    Parameter::List& localParams = paramsPtr->sublist("Linear Solver");

    // Set more parameters that need info from the Group....
    probSize = soln.getX().length(); 

    // Manipulate some parameters further....
    if (kmax > probSize) {
      kmax = probSize;

      // Update parameter list with the actual value used
      localParams.setParameter("Max Iterations", kmax);
    }
    maxDim = kmax + maxp;

    /* may need to delete these pointers if reset more than once - bwb */
    basisVptr = new Abstract::Vector* [maxDim];
    for (int i=0; i<maxDim; i++)
      basisVptr[i] = soln.getX().clone(ShapeCopy);

    allocate_matrix(maxDim, kmax, hess);
    for (int i=0; i<maxDim; i++)
      for (int j=0; j<kmax; j++)
	hess[i][j] = 0;
    allocate_matrix(maxp+1, maxDim, givensC);
    allocate_matrix(maxp+1, maxDim, givensS);

    allocate_vector(maxDim, vecg);
    allocate_vector(maxDim, vecq);
    allocate_vector(maxDim, vecz);
    allocate_vector(maxDim, newHessCol);

    vecw = soln.getX().clone(ShapeCopy);
    dTensor = soln.getX().clone(ShapeCopy);
    dTLambda = soln.getX().clone(ShapeCopy);
    dNewton = soln.getX().clone(ShapeCopy);
    
    isFreshlyReset = false;
  }

  // Set iteration-specific local parameters....
  bool breakdown = false;
  bool ok = true;
  double w1 = 0;           // temporary variable for Givens rotations
  double w2 = 0;           // temporary variable for Givens rotations
  double qa = 0;           // quadratic equation coefficient
  double qb = 0;           // quadratic equation coefficient
  double qc = 0;           // quadratic equation coefficient
  double qval = 0;         // value of quadratic equation at minimizer
  double errTol = 0;       // stopping tolerance for error 

  y1 = 0;           // smallest magnitude root/minimizer of quadratic
  lambdaBar = 1.0;  // quadratic equation threshhold for a real root
  isTensorCalculated = false;
  isCLParamsCalculated = false;

  // Initialize storage (this may not be needed, but it is good practice)...
  for (int i=0; i<maxDim; i++) {
    vecg[i] = 0;
    vecq[i] = 0;
    newHessCol[i] = 0;
    for (int j=0; j<maxp+1; j++) {
      givensC[j][i] = givensS[j][i] = 0;
    }
  }

  // Compute F at current solution....
  ok = soln.computeF();
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
  normS = basisVptr[0]->norm();

  // Calculate the Newton step if this is the first iteration...
  if (normS == 0) {
    ok = soln.computeNewton(paramsPtr->sublist("Linear Solver"));
    if (!ok) {
      if (Utils::doPrint(Utils::Warning))
	cout << "NOX::Direction::Tensor::compute - "
	     << "Unable to compute Newton direction." << endl;
      return false;
    }
    *dNewton = soln.getNewton();
    dTensor->init(0);
    dTLambda->init(0);
    dir = *dNewton;
    return true;
  }


#ifdef CHECK_RESIDUALS
  // Compute the tensor term ac....
  Abstract::Vector* acPtr = soln.getF().clone(ShapeCopy);
  Abstract::Vector* jcxscPtr = soln.getF().clone(ShapeCopy);
  soln.applyJacobian(*basisVptr[0], *jcxscPtr);
  *acPtr = *jcxscPtr;
  acPtr->update(1.0, solver.getPreviousSolutionGroup().getF(), -1.0);
  acPtr->update(-1.0, soln.getF(), 1.0);    
  acPtr->scale(1/(normS*normS*normS*normS));
#endif


  // Compute the tensor term ac, if needed....
  if (precondition == Right) {
    ok = soln.applyJacobian(*basisVptr[0], *basisVptr[1]);
    if (!ok) {
      if (Utils::doPrint(Utils::Warning))
	cout << "NOX::Direction::Tensor::compute - "
	     << "Unable to apply Jacobian." << endl;
      return false;
    }
    basisVptr[1]->update(1.0, solver.getPreviousSolutionGroup().getF(),
		  -1.0, soln.getF(), -1.0);    
    basisVptr[1]->scale(1/(normS*normS*normS*normS));
  }

  // Continue processing the vector sc for right preoconditioning....
  if (precondition == Right) {
    *vecw = *basisVptr[0];
    ok = soln.applyRightPreconditioning(paramsPtr->sublist("Linear Solver"), 
					*vecw, *basisVptr[0]);
    if (!ok) {
      if (Utils::doPrint(Utils::Warning))
	cout << "NOX::Direction::Tensor::compute - "
	     << "Unable to apply Preconditioning." << endl;
      return false;
    }
    normS = basisVptr[0]->norm();
  }



#if DEBUG_LEVEL > 1
  cout << "Current Point: " << soln.getX().length() 
       << "  " << soln.getX().norm() << "  "; 
  soln.getX().print();
  cout << "Previous direction: " << basisVptr[0]->length() 
       << "  " << normS << "  "; 
  basisVptr[0]->print();
  cout << "Current F: ";
  soln.getF().print();
  cout << "Previous F: "; 
  solver.getPreviousSolutionGroup().getF().print();
#ifdef CHECK_RESIDUALS
  cout << "Jcxsc: " << jcxscPtr->length() << "  " 
       << jcxscPtr->norm() << "  "; 
  jcxscPtr->print();
  cout << "Tensor term ac: " << acPtr->length() 
       << "  " << acPtr->norm() << "  ";
  acPtr->print();
#endif
#endif

  // Construct initial n x 3 matrix...
  if (precondition == Left) {
    ok = soln.applyRightPreconditioning(paramsPtr->sublist("Linear Solver"), 
				solver.getPreviousSolutionGroup().getF(), 
				*basisVptr[1]);
    if (!ok) {
      if (Utils::doPrint(Utils::Warning))
	cout << "NOX::Direction::Tensor::compute - "
	     << "Unable to apply Preconditioning." << endl;
      return false;
    }

    ok = soln.applyRightPreconditioning(paramsPtr->sublist("Linear Solver"), 
					soln.getF(), *basisVptr[2]);
    if (!ok) {
      if (Utils::doPrint(Utils::Warning))
	cout << "NOX::Direction::Tensor::compute - "
	     << "Unable to apply Preconditioning." << endl;
      return false;
    }
  }
  else if (precondition == Right) {
    // *basisVptr[1] is calculated above
    *basisVptr[2] = soln.getF();
  }
  else {
    *basisVptr[1] = solver.getPreviousSolutionGroup().getF();
    *basisVptr[2] = soln.getF();
  }


  // Compute the error tolerance, tol*||Fc||  or  tol*||Minv*Fc||
  errTol = tol * basisVptr[2]->norm();
  
#ifdef BUILD_UP_BASIS
  // Build up an orthonormal basis with up to 3 vectors
  double** factorR = NULL;
  allocate_matrix(maxp, maxp, factorR);
  for (int i=0; i<maxp; i++)
    for (int j=0; j<maxp; j++)
	factorR[i][j] = 0;

  int p = 0;
  for (int j=0; j<3; j++) {
    if (j != p)
      *basisVptr[p] = *basisVptr[j];
    for (int i=0; i<p; i++) {
      factorR[i][j] = basisVptr[j]->dot(*basisVptr[i]);
      basisVptr[j]->update(-factorR[i][j], *basisVptr[i], 1.0);
    }
    factorR[p][j] = basisVptr[p]->norm();
    if (factorR[p][j] == 0.0) {
      cout << p+1 << " vectors are linearly dependent\n";
    }
    else {
      basisVptr[p]->scale(1.0/factorR[p][j]);
      p++;
    }
  }

#if DEBUG_LEVEL > 1
  cout << "p = " << p << endl;
  cout << "factorR matrix:\n";
  print_matrix(maxp,maxp,factorR);
#endif

  normS = factorR[0][0];
  for (int i=0; i<maxp; i++) {
    vecq[i] = factorR[i][1];
    vecg[i] = factorR[i][2];
  }

  delete_matrix(factorR);

#else

  // Perform partial QR factorization of initial n x 3 matrix...
  basisVptr[0]->scale(1.0/normS);

  vecq[0] = basisVptr[1]->dot(*basisVptr[0]);
  basisVptr[1]->update(-vecq[0], *basisVptr[0], 1.0);
  vecq[1] = basisVptr[1]->norm();
  if (vecq[1] == 0.0) {
    cout << "sc and Fp are linearly dependent\n";
  }
  else
    basisVptr[1]->scale(1.0/vecq[1]);


  vecg[0] = basisVptr[2]->dot(*basisVptr[0]);
  basisVptr[2]->update(-vecg[0], *basisVptr[0], 1.0);
  vecg[1] = basisVptr[2]->dot(*basisVptr[1]);
  basisVptr[2]->update(-vecg[1], *basisVptr[1], 1.0);
  vecg[2] = basisVptr[2]->norm();
  if (vecg[2] == 0.0) {
    cout << "sc Fp and Fc are linearly dependent\n";
  }
  else
    basisVptr[2]->scale(1.0/vecg[2]);

  int p = maxp;
#endif  // BUILD_UP_BASIS


#ifdef CHECK_ORTHOGONALITY
  //basisVptr[0]->print();
  for (int i=0; i<3; i++) {
    //printf("norm(%d) = %e\n", i, basisVptr[i]->norm(NOX::Abstract::Vector::TwoNorm));
    for (int j=i; j<3; j++) {
      double temp = basisVptr[i]->dot(*basisVptr[j]);
      printf("<v%d,v%d> = %e\n", i, j, temp);
    }
  }
#endif

#if DEBUG_LEVEL > 1
  cout << "normS = " << Utils::sci(normS) << endl;
  cout << "vecg ";
  print_vector(maxDim, vecg);
  cout << "vecq ";
  print_vector(maxDim, vecq);
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
#endif
  
#ifdef STORE_HESSENBERG
  double** hess2 = NULL;
  allocate_matrix(maxDim, kmax, hess2);
  for (int i=0; i<maxDim; i++)
    for (int j=0; j<kmax; j++)
      hess2[i][j] = 0;
#endif

  
  // Calculate and print initial residual errors...
  double nerr = norm(p,vecg);
  double terr = nerr;         // not quite right if x0 != 0  - bwb
  double error = terr;

  // Begin iterating...
  int k = -1;
  while((error > errTol)  &&  (k+1 < kmax)) {
    k++;
    arnoldiIters++;

    // Print residual errors...
    if (k % outputFreq == 0) 
      printf("Iteration: %2d  Tensor residual: %8e  Newton residual: %8e\n", 
	     k, terr, nerr);

    // Begin Arnoldi process...
    if (precondition == Left) {
      ok = soln.applyJacobian(*basisVptr[k], *vecw);
      if (!ok) {
	if (Utils::doPrint(Utils::Warning))
	  cout << "NOX::Direction::Tensor::compute - "
	       << "Unable to apply Jacobian." << endl;
	return false;
      }

      ok = soln.applyRightPreconditioning(paramsPtr->sublist("Linear Solver"),
					  *vecw, *basisVptr[k+p]);
      if (!ok) {
	if (Utils::doPrint(Utils::Warning))
	  cout << "NOX::Direction::Tensor::compute - "
	       << "Unable to apply Preconditioning." << endl;
	return false;
      }
    }
    else if (precondition == Right) {
      ok = soln.applyRightPreconditioning(paramsPtr->sublist("Linear Solver"),
					  *basisVptr[k], *vecw);
      if (!ok) {
	if (Utils::doPrint(Utils::Warning))
	  cout << "NOX::Direction::Tensor::compute - "
	       << "Unable to apply Preconditioning." << endl;
	return false;
      }

      ok = soln.applyJacobian(*vecw, *basisVptr[k+p]);
      if (!ok) {
	if (Utils::doPrint(Utils::Warning))
	  cout << "NOX::Direction::Tensor::compute - "
	       << "Unable to apply Jacobian." << endl;
	return false;
      }
    }
    else {
      ok = soln.applyJacobian(*basisVptr[k], *basisVptr[k+p]);
      if (!ok) {
	if (Utils::doPrint(Utils::Warning))
	  cout << "NOX::Direction::Tensor::compute - "
	       << "Unable to apply Jacobian." << endl;
	return false;
      }
    }
    double normav = basisVptr[k+p]->norm();
  
    // Perform modified Gram-Schmidt....
    for (int j=0; j<k+p; j++) {
      newHessCol[j] = basisVptr[k+p]->dot(*basisVptr[j]);
      basisVptr[k+p]->update(-newHessCol[j], *basisVptr[j], 1.0);
    }
    newHessCol[k+p] = basisVptr[k+p]->norm();
    double normav2 = newHessCol[k+p];

    // Reorthogonalize against first column only....
    for (int j=0; j<1; j++) {
      double hr = basisVptr[k+p]->dot(*basisVptr[j]);
      basisVptr[k+p]->update(-hr, *basisVptr[j], 1.0);
    }
    newHessCol[k+p] = basisVptr[k+p]->norm();

    // Then reorthogonalize against all other columns, if necessary....
    if ((reorth == Always) ||
	((reorth ==AsNeeded) &&
	 (normav2 == 0 || log10(normav/normav2) > 10))) {
      cout << "Reorthogonalize...\n";
      for (int j=1; j<k+p; j++) {
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
    cout << "newHessCol ";
    print_vector(k+p+1, newHessCol);
#endif 

#ifdef STORE_HESSENBERG
    for (int i=0; i<maxDim; i++)
      hess2[i][k] = newHessCol[i];
#endif


    // Calculate projected tensor term using a shortcut method...
    if (k==0  &&  precondition != Right) {
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
    for (int i=0; i<maxDim; i++)
      hess[i][k] = newHessCol[i];
    
    // Find smallest magnitude root of the quadratic equation (first row)...
    qa = vecq[0] * normS * normS;
    qb = hess[0][0];
    qc = vecg[0];
    double discriminant = qb*qb - 4*qa*qc;
    if (discriminant < 0) {
      y1 = -qb/qa/2;
      qval = (qa*y1*y1 + qb*y1 + qc) * hess[0][0];
      lambdaBar = qb*qb / (4*qa*qc);
#if DEBUG_LEVEL > 0
      cout << "  Discriminant is negative!   (LambdaBar = " 
	   << lambdaBar << ")\n";
#endif
    }
    else {
      qval = 0;
      lambdaBar = 1.0;
      if (fabs(qa/qb) < 1e-8) {
#if DEBUG_LEVEL > 0
	cout << "qa is relatively small\n";
#endif 
	y1 = -qc/qb;
      }
      else {
	double tmp1 = (-qb + sqrt(discriminant)) / (2*qa);
	double tmp2 = (-qb - sqrt(discriminant)) / (2*qa);
	y1 = (fabs(tmp1) < fabs(tmp2)) ? tmp1 : tmp2;
      }
    }	

    // Update the tensor residual norm...
    double beta2 = normS * normS * y1 * y1;
    for (int i=0; i<p; i++) {
      terrvec[i] = -vecg[k+i+1] - vecq[k+i+1] * beta2;
    }
    terr = norm(p,terrvec);

    // Update the Newton residual norm...
    for (int i=0; i<p; i++)
      terrvec[i] = -vecg[k+i+1];
    nerr = norm(p,terrvec);

#ifdef USE_NEWTON_DIR
      error = nerr;
#else
      error = terr;
#endif

  }  // end while loop

  // Print final residual errors...
  printf("Iteration: %2d  Tensor residual: %8e  Newton residual: %8e\n", 
	 k+1, terr, nerr);

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
  
  iterations = k+1;
  
  // Solve linear system for Newton direction...
  /* ...but we can't update Newton vector in group... */
  int* pindex = NULL;
  allocate_vector(iterations, pindex);
  pindex[iterations-1] = 0;
  for (int i=0; i<iterations-1; i++)
    pindex[i] = i+1;
  double* yn = backsolve(hess, vecg, pindex, iterations);
#if DEBUG_LEVEL > 1
  cout << "yn ";
  print_vector(iterations, yn);
#endif
  dNewton->init(0);
  for (int i=0; i<iterations; i++)
    dNewton->update(-yn[i], *basisVptr[i], 1);
  if (precondition == Right) {
    *vecw = *dNewton;
    ok = soln.applyRightPreconditioning(paramsPtr->sublist("Linear Solver"),
					*vecw, *dNewton);
    if (!ok) {
      if (Utils::doPrint(Utils::Warning))
	cout << "NOX::Direction::Tensor::compute - "
	     << "Unable to apply Preconditioning." << endl;
      return false;
    }
  }

#if DEBUG_LEVEL > 1
  cout << "Newton Direction 1: ";
  dNewton->print();
#endif


#if DEBUG_LEVEL > 1
  if (predf == NULL) {
    predf = soln.getF().clone(ShapeCopy);
  }
  soln.applyJacobian(*dNewton, *predf);    
  predf->update(1.0, soln.getF(), 1.0);
  double residual = predf->norm();
  printf("Actual Newton1 residual = %8e\n", residual);
#endif

  
  // Compute the Newton direction directly (temporary code)
  ok = soln.computeNewton(paramsPtr->sublist("Linear Solver"));

  // Computing the Newton step failed, but maybe the step is okay to use...
  if (!ok) {
    double accuracy = soln.getNormNewtonSolveResidual();

    if (accuracy < 0) {
      cerr << "NOX::Direction::Newton::compute " 
           << "- getNormNewtonSolveResidual returned a negative value" 
	   << endl;
    }
 
    if (accuracy < soln.getNormF() ) {
      ok = true;
      double tolerance = paramsPtr->
	sublist("Linear Solver").getParameter("Tolerance", 1.0e-10);
      if (Utils::doPrint(Utils::Warning)) 
	cout << "WARNING: NOX::Direction::Tensor::compute - "
	     << "Newton solve failure.\n" 
	     << "Desired accuracy is " 
	     << Utils::sci(tolerance) << ".\n"
	     << "Using solution with accuracy of " << Utils::sci(accuracy) 
	     << "." << endl;
    }
  }
  if (!ok) {
    if (Utils::doPrint(Utils::Warning))
      cout << "NOX::Direction::Tensor::compute - "
	   << "Unable to compute backup Newton direction." << endl;
    return false;
  }

#if DEBUG_LEVEL > 1
  cout << "Newton Direction 2: ";
  soln.getNewton().print();
  double residual2 = soln.getNormNewtonSolveResidual();
  printf("Actual Newton2 residual = %8e\n", residual2);
  double stdn = basisVptr[0]->dot(soln.getNewton()) * normS;
  printf("sTdn =  %e\n", stdn);
#endif


  // Solve system for tensor direction...
#ifndef USE_NEWTON_DIR
  double beta2 = normS * normS * y1 * y1;
  for (int i=0; i<maxDim; i++) {
    vecz[i] = -vecg[i] - vecq[i] * beta2 - hess[i][0] * y1;
  }
  double* yt = backsolve(hess, vecz, pindex, iterations-1, iterations);
  yt[0] = y1;
  dTensor->init(0);
  for (int i=0; i<iterations; i++)
    dTensor->update(yt[i], *basisVptr[i], 1);
  if (precondition == Right) {
    *vecw = *dTensor;
    ok = soln.applyRightPreconditioning(paramsPtr->sublist("Linear Solver"),
					*vecw, *dTensor);
    if (!ok) {
      if (Utils::doPrint(Utils::Warning))
	cout << "NOX::Direction::Tensor::compute - "
	     << "Unable to apply Preconditioning." << endl;
      return false;
    }
  }
  isTensorCalculated = true;
  isCLParamsCalculated = false;


#if DEBUG_LEVEL > 1
  cout << "Tensor Direction: "; 
  dTensor->print();
#endif
#if DEBUG_LEVEL > 1
  double beta = basisVptr[0]->dot(*dTensor) * normS;
  printf("Beta %e, ||sc||*y1 %e\n", beta, normS*y1);
#ifdef CHECK_RESIDUALS
  if (predf == NULL) {
    predf = soln.getF().clone(ShapeCopy);
  }
  soln.applyJacobian(*dTensor, *predf);    
  predf->update(1.0, soln.getF(), beta*beta, *acPtr, 1.0);
  double residual3 = predf->norm();
  printf("Actual tensor residual = %8e (using beta = %e)\n", residual3, beta);
  if ( lambdaBar == 1.0  &&  fabs((residual3 - terr)/terr) > 1e-3) 
    printf("  *** Warning - check residuals ***\n");
#endif  // CHECK_RESIDUALS
#endif  // DEBUG_LEVEL


  // Set output parameters to pass back to calling function
  NOX::Parameter::List& outputList = paramsPtr->sublist("Output");
  outputList.setParameter("Arnoldi iterations", arnoldiIters);

  // Compute parameters needed for curvilinear linesearch
  stJinvF = normS*qc/qb;
  stJinvA = 2*qa/qb/normS;
#if DEBUG_LEVEL > 1
  printf("stJinvF = %e  ", stJinvF);
  printf("stJinvA = %e\n", stJinvA);
#endif

#endif // USE_NEWTON_DIR

#if DEBUG_LEVEL > 1
  cout << "Hessenberg ";
  print_matrix(maxDim, kmax, hess);
  cout << "vecg ";
  print_vector(maxDim, vecg);
  cout << "vecq ";
  print_vector(maxDim, vecq);
  printf("normS = %e, y1 = %e\n", normS, y1);
  cout << "vecz ";
  print_vector(maxDim, vecz);
  cout << "pindex ";
  print_vector(iterations, pindex);
  cout << "yt ";
  print_vector(iterations, yt);
#endif
  
  // Set search direction...
#ifdef USE_NEWTON_DIR
  dir = *dNewton;
#else
  dir = *dTensor;
#endif

  // For testing linesearch...
  //dir = soln.getNewton();
  //dir.scale(-1.0);

  // Cleanup memory
  delete [] pindex;
  delete [] yn;
#ifndef USE_NEWTON_DIR
  delete [] yt;
#endif
#ifdef CHECK_RESIDUALS
  delete acPtr;
  delete jcxscPtr;
#endif
#ifdef STORE_HESSENBERG
  delete_matrix(hess2);
#endif

  return ok;
}


bool Tensor::computeCurvilinearStep(Abstract::Vector& dir, 
				    Abstract::Group& soln, 
				    const Solver::Generic& solver, 
				    double lambda)
{
  bool ok = true;

  // If the tensor step has not been computed, then return the Newton
  // step because other necessary vectors have not been computed.
  // This is most certainly the case for the first nonlinear
  // iteration, where the Newton step is used.
  if (!isTensorCalculated) {
    dir = *dNewton;
    dir.scale(lambda);
    return ok;
  }

  // If the curvilinear step dt(lambdaBar) and 2 parameters stJinvF
  // and stJinvA have not been computed, then calculate this step and
  // parameters for future use...
  if (!isCLParamsCalculated) {

    if (lambdaBar == 1.0) {
      *dTLambda = *dTensor;
    }
    else {
      printf("lambdaBar = %e\n", lambdaBar);

      // Setup index array
      int* pindex = NULL;
      allocate_vector(iterations, pindex);
      pindex[iterations-1] = 0;
      for (int i=0; i<iterations-1; i++)
	pindex[i] = i+1;

      // Setup new right hand side and solve
      double beta2 = normS * normS * y1 * y1;
      for (int i=0; i<maxDim; i++) {
	vecz[i] = -lambdaBar * vecg[i] - vecq[i] * beta2 - hess[i][0] * y1;
      }
      double* yt = backsolve(hess, vecz, pindex, iterations-1, iterations);
      yt[0] = y1;

      // Compute the step: basis vectors * yt
      dTLambda->init(0);
      for (int i=0; i<iterations; i++)
	dTLambda->update(yt[i], *basisVptr[i], 1);

      // Precondition the answer, if preconditioning from the right...
      if (precondition == Right) {
	*vecw = *dTLambda;
	ok = soln.applyRightPreconditioning(paramsPtr->
					    sublist("Linear Solver"),
					    *vecw, *dTLambda);
	if (!ok) {
	  if (Utils::doPrint(Utils::Warning))
	    cout << "NOX::Direction::Tensor::compute - "
		 << "Unable to apply Preconditioning." << endl;
	  return false;
	}
      }
      delete yt;
      delete pindex;
    }

    // Compute the 2 curvilinear parameters that are needed later.
    // The quantities s'*inv(J)*Fc and s'*inv(J)*a are used in the
    // quadratic formula for calculating beta = s'*dTensor.  Here we
    // back-calculate these quantities.
    vecw->update(1.0, *dTLambda, -lambdaBar, *dNewton, 0);
    double c = basisVptr[0]->dot(*vecw) * normS;
    stJinvF = -basisVptr[0]->dot(*dNewton) * normS;
    beta = -lambdaBar * stJinvF + c;
    stJinvA = -2*c / (beta*beta);

#if DEBUG_LEVEL > 1
    printf("stJinvF = %e  ", stJinvF);
    printf("stJinvA = %e\n", stJinvA);
    printf("s'*dt ~=  %e  %e\n", normS*y1, beta);
    printf("c = %e\n", c);
#endif

    isCLParamsCalculated = true;
  }

  // When computing the curvilinear step, the step may be broken up
  // into the following vectors with associated coefficients
  //   dT(lambda) = lambda*dNewton + betaRatio^2*dQuadratic + alpha*dNonroot 
  // where betaRatio and alpha are both functions of lambda.  The
  // computation of this step may be recast in terms of the 3 vectors
  // dNewton, dTensor, and dTLambda, which have been precomputed.
  // Computing the coefficients of these three vectors will save
  // on the computation and storage of dQuadratic and dNonroot.

  // Compute the alpha coefficient, which multiplies the "non-root"
  // vector.  If lambdaBar < 1, then the tensor step is not a root of
  // the model and alpha is nonzero if lambda in [lambdaBar, 1].
  double alpha = 0;
  if (lambdaBar < 1.0) {
    double temp = (lambda - lambdaBar) / (1 - lambdaBar);
    alpha = (temp > 0) ? temp : 0;
  }

  // Compute betaRatio coefficient, which multiplies the "quadratic"
  // correction vector.
  double betaRatio = 0;
  if (lambda > lambdaBar) 
    betaRatio = 1;
  else {
    double discriminant = 0;
    discriminant = 1 - 2*stJinvF*stJinvA*lambda;
    if (discriminant < 0) {
      cout << "Warning: discriminant is negative ("<< discriminant << ").\n";
      discriminant = 0;
    }
    double tmp1 = (-1 + sqrt(discriminant));
    discriminant = 1 - 2*stJinvF*stJinvA*lambdaBar;
    if (discriminant < 0) {
      discriminant = 0;
    }
    double tmp2 = (-1 + sqrt(discriminant));
    betaRatio = tmp1/tmp2;
  }

  // Now recast the coefficients to work with the vectors dN, dT, dTlambda
  double coefNewt = lambda - betaRatio*betaRatio*lambdaBar - 
    alpha*(1-lambdaBar);
  double coefTlam = betaRatio*betaRatio - alpha;
  double coefTens = alpha;

#if DEBUG_LEVEL > 1
  printf("betaRatio = %e\n", betaRatio);
  printf("lambdaBar = %e\n", lambdaBar);
  printf("lambda = %e\n", lambda);
  printf("coef: %e  %e  %e\n", coefNewt, coefTlam, coefTens);
#endif

  // Compute the curvilinear step
  dir = *dNewton;
  if (coefTens != 0.0) 
    dir.update(coefTlam, *dTLambda, coefTens, *dTensor, coefNewt);
  else
    dir.update(coefTlam, *dTLambda, coefNewt);

  return ok;
}



// private

void** Tensor::allocate_matrix(int rows, int cols, double**& a)
{
  if (a) {
    // delete_matrix(a);
    cout << "Warning: Possibly a previously allocated matrix\n";
  }

  // allocate memory for storing a rectangular cols x rows matrix
  a = new double* [rows];
  if (!a) printf("Allocation error in allocate_matrix()\n");
  a[0] = new double [rows*cols];
  if (!a) printf("Allocation error in allocate_matrix()\n");

  for (int i=1; i<rows; i++) 
    a[i] = a[i-1] + cols;

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
  delete [] A[0];
  delete [] A;
}


void Tensor::print_matrix(int rows, int cols, double** A)
   /*  Print a matrix in conventional form.  */
{
  int i,j;

  for (i=0; i<rows; i++) 
    {
      for (j=0; j<cols;j++)
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
  int dimension;
  
  dimension = (dim>n) ? dim : n;

  //  Allocate memory for solution vector
  allocate_vector(dimension, x);

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
