#ifdef WITH_PRERELEASE
//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

/*  Notes:
**
**  (.) Need to get right preconditioning working with the ifpack
**  preconditioners.
**
**  (.) Need to look into breakdown issues with block method.  A good
**  test for that is Broyden(100,0.99) with x0 = -100, either
**  with/without dir.init(0) line, and forcing 100 inner iterations of
**  the tensor method and limiting the outer iterations to 2.
**
**  (.) Need to look at Feng-Pulliam method with right preconditioning
**  and nonzero initial guess.  Right now it appears to work, but it
**  seems incorrect because of the following 2 lines:
**     //vecz[iterations] = sPtr->dot(dir)/normU;
**     vecz[iterations] = dir0xsc/normS/normU;  // bwb - works with right prec
**  It seems the first should be used based on Dan's code and his defn of ss.
**  Also, if d0_in_subspace_Vm then the line
**     double cc = inner_product(iterations, vecz, yn) - dir0xsc/normS; 
**  cannot use sPtr->dot(dir), presumably because the dot product must
**  use the original sc and not Minv'*sc, hence dir0xsc/normS.  I bet if I
**  work this out on paper, then I would find that d0'*sc should use
**  original sc and not Minv'*sc because d0 added in after right
**  preconditioning dTensor at end.
**
**  (.) Perhaps move vector yt to same level as yn and rewrite
**  backsolve so that allocated vector is passed into subroutine.
**  Then I can move computation of dTensor in F-P step to be outside
**  yt calculation section.
**
**  (.) Reorthogonalization looks like it doesn't include update to
**  Hessenberg.  I added the correct code but left it out so that it
**  won't affect base cases.
**
**  (.) Simplify this whole code by separating out some subroutines
**  for better readability.
**
**  (.) For TensorStep2 and 2+ the value of sc'*dt does not equal
**  normS*y1+dir0xsc.  I don't think this is an error in the code but
**  algorithmically is not equivalent due to the orthogonal matrix Q1.
**  Maybe look into this issue at some point.
*/


/*  Debug Levels:  (not correct yet)
**    0:  No debug information printed.
**    1:  Print residuals
**    2:  Print vectors and matrices
*/

#define DEBUG_LEVEL 0
#undef CALCULATE_GRP_NEWTON
#undef STORE_HESSENBERG
#undef CHECK_ORTHOGONALITY
#undef CHECK_RESIDUALS

#include "NOX_Direction_Tensor.H" // class definition
#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"
#include "stdio.h"  // for printf()


NOX::Direction::Tensor::Tensor(const NOX::Utils& u,
			       NOX::Parameter::List& params) :
  utils(u)

{
  scPtr = NULL;
  acPtr = NULL;
  jcxshatPtr = NULL;
  mtinvscPtr = NULL;
  vecw = NULL;
  dTensor = NULL;
  dTLambda = NULL;
  dNewton = NULL;
  hess = NULL;
  givensC = NULL;
  givensS = NULL;
  hhZ = NULL;
  hq1save = NULL;
  gsave = NULL;
  qsave = NULL;
  qk = NULL;
  q1subdiag = NULL;
  qNorm2 = NULL;
  hk = NULL;
  newHessCol = NULL;
  terrvec = NULL; 
  vecg = NULL;
  vecq = NULL;
  vecz = NULL;
  
  reset(params);
}


NOX::Direction::Tensor::~Tensor()
{
  delete scPtr;
  delete acPtr;
  
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
  if (requestedBaseStep == TensorStep2) {
    delete_matrix(hhZ);
    delete [] hq1save;
    delete [] gsave;
    delete [] qsave;
    delete [] qk;
    delete [] q1subdiag;
    delete [] qNorm2;
    if (precondition == Right && isMinvTransAvailable)
      delete mtinvscPtr;
  }

  if (isSubspaceAugmented) {
    delete [] hk;
    //delete jcxshatPtr;
  }

  delete [] newHessCol;
  delete [] terrvec;
  delete [] vecg;
  delete [] vecq;
  delete [] vecz;
}


bool NOX::Direction::Tensor::reset(NOX::Parameter::List& params)
{
  /*  bwb -
   *  May need code if repeatedly calling reset().
   *  Memory leaks might be an issue here.
   */

  paramsPtr = &params;

  NOX::Parameter::List& p = paramsPtr->sublist("Tensor");
  doRescue = p.getParameter("Rescue Bad Newton Solve", true);
  
  NOX::Parameter::List& localParams = p.sublist("Linear Solver");
  localParamsPtr = &(p.sublist("Linear Solver"));
  tol = localParams.getParameter("Tolerance", 1.0e-4);
  kmax = localParams.getParameter("Size of Krylov Subspace", 30);
  maxRestarts = localParams.getParameter("Max Restarts", 0);
  outputFreq = localParams.getParameter("Output Frequency", 20);
  isSubspaceAugmented = false;

  string choice = localParams.getParameter("Compute Step", "Tensor");
  if (choice == "Tensor") {
    requestedBaseStep = TensorStep3;
  }
  else if (choice == "Tensor2") {
    cout << "Tensor step 2 is requested\n";
    requestedBaseStep = TensorStep2;
  }
  else if (choice == "Tensor2+") {
    cout << "Tensor step 2+ is requested\n";
    requestedBaseStep = TensorStep2;
    isSubspaceAugmented = true;
  }
  else if (choice == "Feng Pulliam") {
    cout << "Feng-Pulliam tensor step is requested\n";
    requestedBaseStep = TensorStepFP;
  }
  else if (choice == "Newton") {
    cout << "Newton step is requested\n";
    requestedBaseStep = NewtonStep;
  }
  else {
    cout << "Warning: NOX::Direction::Tensor::reset() - the choice of "
	 << "\"Compute Step\" parameter is invalid." << endl;
    requestedBaseStep = TensorStep3;
  }


  choice = localParams.getParameter("Reorthogonalize", "As Needed");
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
    NOX::Direction::Tensor::throwError("reset",
	    "the choice of \"Reorthogonalization\" parameter is invalid.");
  }


  choice = localParams.getParameter("Preconditioning", "None");
  if (choice == "None") {
    precondition = None;
    localParams.setParameter("Preconditioning Side", "None");
  }
  else {
    choice = localParams.getParameter("Preconditioning Side", "Left");
    if (choice == "Left") {
      precondition = Left;
    }
    else if (choice == "Right") {
      precondition = Right;   // bwb - does not work yet if M is not symmetric
      cout << "Warning: Right preconditioning does not work yet if "
	   << "the preconditioner is asymmetric." << endl;
    }
    else if (choice == "None") {
      precondition = None;
    }
    else {
      NOX::Direction::Tensor::throwError("reset",
	    "the choice of \"Preconditioning Side\" parameter is invalid.");
    }
  }

  // Determine whether to use a shortcut method for computing rhs V'*ac
  useShortcutMethod = false;
  if (requestedBaseStep == TensorStep3  &&  precondition != Right)
    useShortcutMethod = localParams.getParameter("Use Shortcut Method", true);
  
  isMinvTransAvailable = false;
  
  if (requestedBaseStep == TensorStep3)
    pmax = 3;   // modified block GMRES 
  else if (requestedBaseStep == TensorStep2)
    pmax = 2;
  else
    pmax = 1;   // standard GMRES
  
  allocate_vector(pmax, terrvec);

  maxDim = 0;   //  updated in compute() once n is known

  arnoldiIters = 0;

  isFreshlyReset = true;
  isTensorCalculated = false;
  isCLParamsCalculated = false;

  return true;
}


bool NOX::Direction::Tensor::compute(NOX::Abstract::Vector& dir, 
				     NOX::Abstract::Group& soln, 
				     const NOX::Solver::Generic& solver)
{
  // Continue with the reset() function here now that more information
  // is known via "soln".
  if (isFreshlyReset) {

    // Set more parameters that need info from the Group....
    probSize = soln.getX().length(); 

    // Manipulate some parameters further....
    if (kmax > probSize) {
      kmax = probSize;

      // Update parameter list with the actual value used
      localParamsPtr->setParameter("Max Iterations", kmax);
    }

    maxDim = kmax + pmax;
    if (isSubspaceAugmented) maxDim++;     // bwb - needed for V?
    
    /* may need to delete these pointers if reset more than once - bwb */
    basisVptr = new NOX::Abstract::Vector* [maxDim];
    for (int i=0; i<maxDim; i++)
      basisVptr[i] = soln.getX().clone(ShapeCopy);
    
    scPtr = soln.getX().clone(ShapeCopy);
    acPtr = soln.getF().clone(ShapeCopy);
    if (requestedBaseStep == TensorStep2 &&
	precondition == Right  &&  isMinvTransAvailable)
      mtinvscPtr = soln.getX().clone(ShapeCopy);
    if (isSubspaceAugmented)
      jcxshatPtr = basisVptr[maxDim-1];    // alias
    //jcxshatPtr = soln.getF().clone(ShapeCopy);
    
    int kmaxh = kmax;
    if (requestedBaseStep == TensorStep2)
      kmaxh++;
    
    allocate_matrix(maxDim, kmaxh, hess);
    for (int i=0; i<maxDim; i++)
      for (int j=0; j<kmaxh; j++)
	hess[i][j] = 0;

    allocate_matrix(pmax+1, maxDim, givensC);
    allocate_matrix(pmax+1, maxDim, givensS);

    if (requestedBaseStep == TensorStep2) {
      allocate_matrix(maxDim, pmax+2, hhZ);
      allocate_vector(maxDim, hq1save);
      allocate_vector(maxDim, gsave);
      allocate_vector(maxDim, qsave);
      allocate_vector(maxDim, qk);
      allocate_vector(maxDim, q1subdiag);
      allocate_vector(maxDim, qNorm2);
    }

    if (isSubspaceAugmented) allocate_vector(maxDim, hk);

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

  // Local variables
  NOX::Abstract::Group::ReturnType status;
  double errTol = 0;       // stopping tolerance for error 
  double error = 0;        // actual error of approximate solution
  int restarts = 0;        // counter for the number of restarts
  
  // Initialize class private variables for this step computation...
  y1 = 0;           // smallest magnitude root/minimizer of quadratic
  lambdaBar = 1.0;  // threshold value for a real root in quadratic equation
  isTensorCalculated = false;
  isCLParamsCalculated = false;

  dir.init(0);    // bwb - remove this line if testing restarts
  
  // Compute F at current solution...
  status = soln.computeF();
  if (status != NOX::Abstract::Group::Ok)
    NOX::Direction::Tensor::throwError("compute", "Unable to compute F");

  // Compute Jacobian at current solution....
  status = soln.computeJacobian();
  if (status != NOX::Abstract::Group::Ok)
    NOX::Direction::Tensor::throwError("compute",
				       "Unable to compute Jacobian");

  // Compute the tensor terms sc and ac....
  *scPtr = soln.getX();
  scPtr->update(1.0, solver.getPreviousSolutionGroup().getX(), -1.0);
  normS = scPtr->norm();
  soln.applyJacobian(*scPtr, *acPtr);
  if (isSubspaceAugmented)  *jcxshatPtr = *acPtr;
  acPtr->update(1.0, solver.getPreviousSolutionGroup().getF(), -1.0);
  acPtr->update(-1.0, soln.getF(), 1.0);    
  acPtr->scale(1/(normS*normS*normS*normS));

#if DEBUG_LEVEL > 1
  cout << "Current Point: " << soln.getX().length() 
       << "  " << soln.getX().norm() << "  "; 
  soln.getX().print();
  cout << "Current F: ";
  soln.getF().print();
  cout << "Previous F: "; 
  solver.getPreviousSolutionGroup().getF().print();
  cout << "Previous direction: " << scPtr->length() 
       << "  " << normS << "  "; 
  scPtr->print();
  cout << "Tensor term ac: " << acPtr->length() 
       << "  " << acPtr->norm() << "  ";
  acPtr->print();
#endif // DEBUG_LEVEL

  // Compute the error tolerance, tol*||Fc||  or  tol*||Minv*Fc||
  if (precondition == Left) {
      status = soln.applyRightPreconditioning(false, *localParamsPtr, 
					      soln.getF(), *vecw);
      if (status != NOX::Abstract::Group::Ok)
	NOX::Direction::Tensor::throwError("compute",
					   "Unable to apply preconditioner.");
    errTol = tol * vecw->norm();
  }
  else
    errTol = tol * soln.getF().norm();
  error = errTol + 1;   // to force entry in while loop

  
  // Iterate until solved or max restarts
  while (error > errTol  &&  restarts <= maxRestarts) {
    compute_step(dir, soln, solver, errTol, error);
#if DEBUG_LEVEL > 0
    printf(" Restart: %d --- error = %e\n\n", restarts, error);
#endif
    restarts++;
  }

  // Set output parameters to pass back to calling function
  NOX::Parameter::List& outputList = paramsPtr->sublist("Output");
  outputList.setParameter("Arnoldi iterations", arnoldiIters);
    
  return true;
}


bool NOX::Direction::Tensor::compute_step(NOX::Abstract::Vector& dir, 
					  NOX::Abstract::Group& soln, 
					  const NOX::Solver::Generic& solver,
					  const double errTol,
					  double& error)
{
  // Set iteration-specific local parameters....
  NOX::Abstract::Group::ReturnType status;
  bool breakdown = false;
  bool isDir0 = true;      // flag for noting whether initial guess is nonzero
  double dir0xsc = 0;      // value of dir0'*sc
  double c1 = 0;           // temporary variable for Givens rotation cos value
  double s1 = 0;           // temporary variable for Givens rotation sin value
  double qa = 0;           // quadratic equation coefficient
  double qb = 0;           // quadratic equation coefficient
  double qc = 0;           // quadratic equation coefficient
  double qval = 0;         // value of quadratic equation at minimizer
  double *yn = NULL;       // vector for linear combination of basis vectors
  int* pindex = NULL;      // permutation vector
  int leftPreListSize = 0; // Size of the index list for left preconditioning
  int* leftPreList = NULL; // List of indices for left preconditioning
  int k1 = 0;              // index holding either 0, k, or k+1
  
  // Local variables for TensorStep2 method.
  double jcxshatNorm2 = 0; // Norm-squared of Jc*shat
  double gamma2 = 0;       // Squared value of last entry in Hessenberg 
  double ztz = 0;          // z'*z where z = Jc*shat
  double omega = 0;
  double qkk = 0;
  double qkkm1 = 0;
  double qkm1Norm2 = 0;
  double qkNorm2 = 0;
  
#ifdef CHECK_RESIDUALS
  double residual = 0;
  double fprime = 0;
#endif
  
  // Initialize storage...
  // (This may not be needed, but just to be safe)
  for (int i=0; i<maxDim; i++) {
    vecg[i] = 0;
    vecq[i] = 0;
    newHessCol[i] = 0;
    for (int j=0; j<pmax+1; j++) {
      givensC[j][i] = givensS[j][i] = 0;
    }
  }

  // Use these default values...
  stepType requestedStep = requestedBaseStep;
  bool isAugmentSubspace = isSubspaceAugmented;
  int maxp = pmax;

  // ...unless this is the first iteration, then calculate Newton step...
  if (normS == 0) {
    maxp = 1;
    requestedStep = NewtonStep;
    isAugmentSubspace = false;
  }

  // If initial guess nonzero, then compute d0'*sc and Jc*d0...
  isDir0 = (dir.norm() == 0);
  if (!isDir0) {
    soln.applyJacobian(dir, *vecw);
    dir0xsc = scPtr->dot(dir);
#if DEBUG_LEVEL > 0
    cout << " Initial guess \"dir\" is NOT zero\n";
    cout << " dir0xsc = " << dir0xsc << "\n";
#endif
  }


  // Build up the initial basis matrix...
  if (requestedStep == TensorStep3) {

    // Use sc as first basis vector...
    *basisVptr[0] = *scPtr;

    if (useShortcutMethod) {
      *basisVptr[1] = solver.getPreviousSolutionGroup().getF();
      if (!isDir0)  basisVptr[1]->update(1.0, *vecw, 1.0);
    }
    else {
      *basisVptr[1] = *acPtr;
    }
    *basisVptr[2] = soln.getF();
    if (!isDir0)  basisVptr[2]->update(1.0, *vecw, 1.0);

    // Set up a list for left-preconditioned vectors...
    leftPreListSize = 2;
    allocate_vector(leftPreListSize, leftPreList);
    leftPreList[0] = 1;
    leftPreList[1] = 2;
  }
  else if (requestedStep == TensorStep2) {
    *basisVptr[0] = soln.getF();
    if (!isDir0)  basisVptr[0]->update(1.0, *vecw, 1.0);
    *basisVptr[1] = *acPtr;

    // Set up a list for left-preconditioned vectors...
    leftPreListSize = 3;
    allocate_vector(leftPreListSize, leftPreList);
    leftPreList[0] = 0;
    leftPreList[1] = 1;
    leftPreList[2] = maxDim-1;  // jcxshatPtr
  }
  else {   // requestedStep == NewtonStep or TensorStepFP
    *basisVptr[0] = soln.getF();
    if (!isDir0)  basisVptr[0]->update(1.0, *vecw, 1.0);

    // Set up a list for left-preconditioned vectors...
    leftPreListSize = 1;
    allocate_vector(leftPreListSize, leftPreList);
    leftPreList[0] = 0;
  }  // end if (requestedStep)


  // Precondition the initial basis vectors...
  if (precondition == Left) {
    for (int i=0; i<leftPreListSize;i++) {
      int indx = leftPreList[i];
      *vecw = *basisVptr[indx];
      status = soln.applyRightPreconditioning(false, *localParamsPtr, 
					      *vecw, *basisVptr[indx]);
      if (status != NOX::Abstract::Group::Ok)
	NOX::Direction::Tensor::throwError("compute",
					   "Unable to apply preconditioner.");
    }
  }
  else if (requestedStep == TensorStep3  &&  precondition == Right) {

    // Continue processing the vector sc for right preconditioning....
    *vecw = *basisVptr[0];
    status = soln.applyRightPreconditioning(true, *localParamsPtr,
					    *vecw, *basisVptr[0]);
    if (status != NOX::Abstract::Group::Ok)
      NOX::Direction::Tensor::throwError("compute",
					 "Unable to apply preconditioner.");
  }
  else if (requestedStep == TensorStep2  &&
	   precondition == Right  && isMinvTransAvailable) {

    // Continue processing the vector sc for right preconditioning....
    *vecw = *scPtr;
    status = soln.applyRightPreconditioning(true, *localParamsPtr,
					    *vecw, *mtinvscPtr);
    if (status != NOX::Abstract::Group::Ok)
      NOX::Direction::Tensor::throwError("compute",
					 "Unable to apply preconditioner.");
  }
  delete [] leftPreList;


  // Work with the vector Jc*sc    
  if (isAugmentSubspace) {
    jcxshatPtr->scale(1/normS);
    jcxshatNorm2 = jcxshatPtr->dot(*jcxshatPtr);
    gamma2 = jcxshatNorm2;
  }


  // Build up an orthonormal basis with up to maxp vectors
  double** factorR = NULL;
  allocate_matrix(maxp, maxp, factorR);
  for (int i=0; i<maxp; i++)
    for (int j=0; j<maxp; j++)
	factorR[i][j] = 0;

  int p = 0;
  for (int j=0; j<maxp; j++) {
    if (j != p)
      *basisVptr[p] = *basisVptr[j];
    for (int i=0; i<p; i++) {
      factorR[i][j] = basisVptr[j]->dot(*basisVptr[i]);
      basisVptr[j]->update(-factorR[i][j], *basisVptr[i], 1.0);
    }
    factorR[p][j] = basisVptr[p]->norm();
    if (factorR[p][j] == 0.0) {
      cout << "Note: vector " << p+1 << " is linearly dependent on others "
	   << "and will not be inlcuded\n";
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

  if (requestedStep == TensorStep3) {
    normS = factorR[0][0];       // Modify ||s|| since s was preconditioned. 
    for (int i=0; i<maxp; i++) {
      vecq[i] = factorR[i][1];
      vecg[i] = factorR[i][2];
    }
  }
  else if (requestedStep == TensorStep2) {
    vecg[0] = factorR[0][0];
    vecq[0] = factorR[0][1];
    vecq[1] = factorR[1][1];
    for (int i=0; i<maxDim; i++) {
      gsave[i] = vecg[i];
      qsave[i] = vecq[i];
    }
  }
  else {
    vecg[0] = factorR[0][0];
    vecq[0] = factorR[0][0];    // bwb - line not needed here.
  }

  delete_matrix(factorR);

#ifdef CHECK_ORTHOGONALITY
  //basisVptr[0]->print();
  for (int i=0; i<maxp; i++) {
    printf("norm(v%d) = %e\n", i,
	   basisVptr[i]->norm(NOX::Abstract::Vector::TwoNorm));
    for (int j=i; j<maxp; j++) {
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

  // Form ghost column of H...
  if (isAugmentSubspace) {
    hk[0] = basisVptr[0]->dot(*jcxshatPtr);
    hk[1] = basisVptr[1]->dot(*jcxshatPtr);
    ztz = hk[0]*hk[0] + hk[1]*hk[1];
  }

  // Initialize residual errors...
  double beta2 = dir0xsc;
  beta2 = beta2 * beta2;
  for (int i=0; i<p; i++) {
    terrvec[i] = -vecg[i] - vecq[i] * beta2;
  }
  double terr = norm(p,terrvec);   // bwb - right?  is vecq correct here?
  double nerr = norm(p,vecg);
  error = terr;
  int k = -1;

  // Test for termination on entry (i.e., dir is satisfactory)
  // bwb - probably still needs work (i.e., dTensor/dNewton not set)
  if (error < errTol) {
    goto Cleanup_before_exit;
  }
  

  // Begin iterating...
  //while(k+1 < kmax  &&  !breakdown) {
  while((error > errTol)  &&  (k+1 < kmax)  &&  !breakdown) {
    k++;
    k1 = k;
    arnoldiIters++;

    // Print residual errors...
    if (k % outputFreq == 0) 
      printf("Iteration: %2d  Tensor residual: %8e  Newton residual: %8e\n", 
	     k, terr, nerr);

    // Compute Jacobian-vector products...
    if (precondition == Left) {
      status = soln.applyJacobian(*basisVptr[k], *vecw);
      if (status != NOX::Abstract::Group::Ok)
	NOX::Direction::Tensor::throwError("compute",
					   "Unable to apply Jacobian.");
      
      status = soln.applyRightPreconditioning(false, *localParamsPtr,
					      *vecw, *basisVptr[k+p]);
      if (status != NOX::Abstract::Group::Ok)
	NOX::Direction::Tensor::throwError("compute",
					   "Unable to apply preconditioner.");
    }
    else if (precondition == Right) {
      status = soln.applyRightPreconditioning(false, *localParamsPtr,
					      *basisVptr[k], *vecw);
      if (status != NOX::Abstract::Group::Ok)
	NOX::Direction::Tensor::throwError("compute",
					   "Unable to apply preconditioner.");

      status = soln.applyJacobian(*vecw, *basisVptr[k+p]);
      if (status != NOX::Abstract::Group::Ok)
	NOX::Direction::Tensor::throwError("compute",
					   "Unable to apply Jacobian.");
    }
    else {
      status = soln.applyJacobian(*basisVptr[k], *basisVptr[k+p]);
      if (status != NOX::Abstract::Group::Ok)
	NOX::Direction::Tensor::throwError("compute",
					   "Unable to apply Jacobian.");
    }
    double normav = basisVptr[k+p]->norm();
  
    // Perform modified Gram-Schmidt....
    for (int j=0; j<k+p; j++) {
      newHessCol[j] = basisVptr[k+p]->dot(*basisVptr[j]);
      basisVptr[k+p]->update(-newHessCol[j], *basisVptr[j], 1.0);
    }
    newHessCol[k+p] = basisVptr[k+p]->norm();

    // Reorthogonalize against first column only....
    for (int j=0; j<1; j++) {
      double hr = basisVptr[k+p]->dot(*basisVptr[j]);
      newHessCol[j] += hr;
      basisVptr[k+p]->update(-hr, *basisVptr[j], 1.0);
    }
    newHessCol[k+p] = basisVptr[k+p]->norm();
    double normav2 = newHessCol[k+p];

    // Then reorthogonalize against all other columns, if necessary....
    if ((reorth == Always) ||
	((reorth ==AsNeeded) &&
	 (normav2 == 0 || log10(normav/normav2) > 10))) {
#if DEBUG_LEVEL > 0
      cout << "Reorthogonalize...\n";
#endif
      for (int j=1; j<k+p; j++) {
	double hr = basisVptr[k+p]->dot(*basisVptr[j]);
        newHessCol[j] += hr;
	basisVptr[k+p]->update(-hr, *basisVptr[j], 1.0);
      }
      newHessCol[k+p] = basisVptr[k+p]->norm();
    }

    // Watch out for breakdown of Arnoldi process...
    if (newHessCol[k+p] != 0  && log10(normav/newHessCol[k+p]) < 16)
      basisVptr[k+p]->scale(1.0/newHessCol[k+p]);
    else {
      // bwb - need to look into breakdown issues with block method
      cout << "Breakdown\n";
      breakdown = true;
    }

#ifdef CHECK_ORTHOGONALITY
    printf("%d -- newHessCol[k+p] = %e   log10(normav/normav2) = %e\n",
	   k, newHessCol[k+p], log10(normav/newHessCol[k+p]));
#endif
#if DEBUG_LEVEL > 2
    cout << "Iteration " << k << ":\n";
    cout << "newHessCol ";
    print_vector(k+p+1, newHessCol);
#endif 
#ifdef STORE_HESSENBERG
    for (int i=0; i<maxDim; i++)
      hess2[i][k] = newHessCol[i];
#endif

    // Transform Hessenberg matrix into upper triangular and solve:

    if (requestedStep == TensorStep2) {
      // Special implementation using Householder reflections and Q2,Q1

      double vkts = 0;
      
      // Construct orthogonal Hessenberg Q1 using 3 vectors...
      if (precondition == Right  &&  isMinvTransAvailable) 
	vkts = basisVptr[k]->dot(*mtinvscPtr);
      else if (precondition == Right) {
	status = soln.applyRightPreconditioning(false, *localParamsPtr, 
						*basisVptr[k], *vecw);
	if (status != NOX::Abstract::Group::Ok)
	  NOX::Direction::Tensor::throwError("compute",
			     "Unable to apply preconditioner.");
	vkts = vecw->dot(*scPtr);
      }
      else
	vkts = basisVptr[k]->dot(*scPtr);
      
      qkk = vkts;
      qkkm1 = -qkNorm2 / qkk;
      qkm1Norm2 = qkNorm2 + qkkm1*qkkm1;
      qkNorm2 = qkNorm2 + qkk*qkk;
      qk[k] = qkk;                 // last column of Q1
      qNorm2[k] = qkNorm2;         // Norm-squared of each column of Q1
      if (k>0) {
	q1subdiag[k-1] = qkkm1;    // subdiagonal of Q1
	qNorm2[k-1] = qkm1Norm2;
      }

      // Form ghost column of H and Q1...
      if (isAugmentSubspace) {
	if (gamma2 > 0)
	  hk[k+p] = basisVptr[k+p]->dot(*jcxshatPtr);
	else
	  hk[k+p] = 0;
	ztz += hk[k+p]*hk[k+p];
	gamma2 = jcxshatNorm2 - ztz;
	if (gamma2 <= 0  ||  k+p+1 >= probSize) {
	  printf("sc is now in span{Vm}    (k = %d)\n", k);
	  isAugmentSubspace = 0;       // bwb - temporary solution
	  gamma2 = 0;
	}
	else {
	  double gamma = sqrt(gamma2);
	  hk[k+p+1] = gamma;
	  omega = -qkNorm2 / normS;
	  q1subdiag[k] = omega;
	  qNorm2[k] = qkNorm2 + omega*omega;
	  qk[k+1] = normS;
	  qNorm2[k+1] = qkNorm2 + normS*normS;
	}	  
      }

      // Form H*Q1 efficiently...
      if (k>0) {
	double tmp = sqrt(qNorm2[k-1]);
	for (int i=0; i<=k+p; i++) 
	  hess[i][k-1] = (hq1save[i] + newHessCol[i]*qkkm1) / tmp;
      }
      double tmp = sqrt(qNorm2[k]);
      for (int i=0; i<=k+p; i++) {
	hess[i][k] = hq1save[i] + newHessCol[i]*qkk;
	hq1save[i] = hess[i][k];
	hess[i][k] = hess[i][k] / tmp;
      }
      if (isAugmentSubspace) {
	double tmp1 = sqrt(qNorm2[k]);
	double tmp2 = sqrt(qNorm2[k+1]);
	for (int i=0; i<=k+p+1; i++) {
	  hess[i][k]   = (hq1save[i] + hk[i]*omega) / tmp1;
	  hess[i][k+1] = (hq1save[i] + hk[i]*normS) / tmp2;
	}
      }

      // Transform to triangular system: Q2*H*Q1
      k1 = k;
      if (k>0) {
	for (int j=0; j<=k-2; j++) {
	  applyHouseholder(hhZ[j], hess, k-1, j, j+p+1);
	  applyHouseholder(hhZ[j], hess, k, j, j+p+1);
	  if (isAugmentSubspace)
	    applyHouseholder(hhZ[j], hess, k+1, j, j+p+1);
	}
	k1 = k-1;
	computeHouseholder(hess, k1, k1, k1+p+1, hhZ[k1]);
	applyHouseholder(hhZ[k1], hess, k1+1, k1, k1+p+1);
	if (isAugmentSubspace)
	  applyHouseholder(hhZ[k1], hess, k1+2, k1, k1+p+1);
	applyHouseholder(hhZ[k1], gsave, k1, k1+p+1);
	applyHouseholder(hhZ[k1], qsave, k1, k1+p+1);
	for (int j=0; j<maxDim; j++) {
	  vecg[j] = gsave[j];
	  vecq[j] = qsave[j];
	}
	k1++;
      }

      // Compute and apply Householder on next column of H*Q1 (temporary)
      if (isAugmentSubspace) {
	computeHouseholder(hess, k1, k1, k1+p+1, hhZ[k1+1]);
	applyHouseholder(hhZ[k1+1], hess, k1+1, k1, k1+p+1);
	applyHouseholder(hhZ[k1+1], vecg, k1, k1+p+1);
	applyHouseholder(hhZ[k1+1], vecq, k1, k1+p+1);
	k1++;
      }

      // Compute and apply Householder on last column of H*Q1 (temporary)
      computeHouseholder(hess, k1, k1, k1+p, hhZ[k1+1]);
      applyHouseholder(hhZ[k1+1], vecg, k1, k1+p);
      applyHouseholder(hhZ[k1+1], vecq, k1, k1+p);

      // Find root of q(beta) equation
      qa = vecq[k1] * qNorm2[k1];
      qb = hess[k1][k1] + 2 * vecq[k1] * dir0xsc * sqrt(qNorm2[k1]);
      qc = vecg[k1] + vecq[k1] * dir0xsc * dir0xsc;
      y1 = calculateBeta(qa, qb, qc, qval, lambdaBar, dir0xsc,
			 sqrt(qNorm2[k1]));
      qval = qval * hess[k1][k1];    // bwb - don't think this is right
      
      // Update the tensor residual norm...
      double beta2 = dir0xsc + sqrt(qNorm2[k1]) * y1;
      beta2 = beta2 * beta2;
      for (int i=0; i<p; i++) {
	terrvec[i] = -vecg[k1+i+1] - vecq[k1+i+1] * beta2;
      }
      terr = norm(p,terrvec);
    }
    else {  // Standard implementation using Givens rotations
      
      // Calculate projected tensor term using a shortcut method...
      if (useShortcutMethod  &&  k==0  &&  requestedStep == TensorStep3) {
	double tmp = normS*normS*normS*normS;
	for (int i=0; i<p+1; i++) {
	  vecq[i] = (vecq[i] - vecg[i] - newHessCol[i]*normS) / tmp;
	}
      }

      // Perform Givens Rotations to put Hessenberg in upper triangular form...
      // ... by doing previous rotations on new Hessenberg column
      if (k>0) {
	for (int i=p-1; i>=0; i--) {
	  for (int j=0; j<k; j++)
	    applyGivens(givensC[i][j], givensS[i][j],
			newHessCol[i+j], newHessCol[i+j+1]);
	}
      }
      // ... and then performing current rotations 
      for (int j=k+p-1; j>=k; j--) {
	computeGivens(newHessCol[j], newHessCol[j+1], c1, s1);
	newHessCol[j] = c1 * newHessCol[j] - s1 * newHessCol[j+1];
	newHessCol[j+1] = 0;
	applyGivens(c1, s1, vecg[j], vecg[j+1]);
	if (requestedStep == TensorStep3)
	  applyGivens(c1, s1, vecq[j], vecq[j+1]);
	givensC[j-k][k] = c1;
	givensS[j-k][k] = s1;
      }
  

      if (requestedStep == TensorStep3) {
	// Eliminate elements in first row of Hessenberg...
	// ... by first performing rotations from previous iterations
	if (k>1) {
	  for (int i=1; i<k; i++) 
	    applyGivens(givensC[p][i], givensS[p][i],
			newHessCol[i], newHessCol[0]);
	}
	// ... and then by performing rotations from the current iteration
	if (k>0) {
	  computeGivens(newHessCol[k], newHessCol[0], c1, s1);
	  newHessCol[k] = c1 * newHessCol[k] - s1 * newHessCol[0];
	  newHessCol[0] = 0;

	  applyGivens(c1, s1, hess[k][0], hess[0][0]);
	  applyGivens(c1, s1, vecg[k], vecg[0]);
	  applyGivens(c1, s1, vecq[k], vecq[0]);
	
	  givensC[p][k] = c1;
	  givensS[p][k] = s1;
	}
      }
	
      // Put the working vector into the Hessenberg matrix...
      for (int i=0; i<maxDim; i++)
	hess[i][k] = newHessCol[i];

      if (requestedStep == TensorStep3) {
	// Find smallest magnitude root of quadratic equation (first row)...
	qa = vecq[0] * normS * normS;
	qb = hess[0][0] + 2 * vecq[0] * dir0xsc * normS;
	qc = vecg[0] + vecq[0] * dir0xsc * dir0xsc;
	y1 = calculateBeta(qa, qb, qc, qval, lambdaBar, dir0xsc, normS);
	qval = qval * hess[0][0];
      
	// Update the tensor residual norm...
	double beta2 = dir0xsc + normS * y1;
	beta2 = beta2 * beta2;
	for (int i=0; i<p; i++) {
	  terrvec[i] = -vecg[k+i+1] - vecq[k+i+1] * beta2;
	}
	terr = norm(p,terrvec);
      } // endif (requestedStep)
    } // endif (requestedStep == TensorStep2)

    
    // Update the Newton residual norm...
    for (int i=0; i<p; i++)
      terrvec[i] = -vecg[k1+i+1];
    nerr = norm(p,terrvec);          // bwb - simplify to one line

    if (requestedStep == NewtonStep || requestedStep == TensorStepFP) 
      error = nerr;
    else
      error = terr;

  }  // end while loop

  //iterations = k + 1;
  iterations = k1 + 1;

  // Print final residual errors...
  printf("Iteration: %2d  Tensor residual: %8e  Newton residual: %8e\n", 
	 k+1, terr, nerr);

#ifdef CHECK_ORTHOGONALITY
  for (int i=0; i<iterations; i++) {
    double temp = basisVptr[i]->norm(NOX::Abstract::Vector::TwoNorm);
    if (fabs(temp-1) > 1e-14)  printf("norm(v%d) = %e\n", i, temp);
    for (int j=i; j<iterations; j++) {
      double temp = basisVptr[i]->dot(*basisVptr[j]);
      if (fabs(temp - (i==j?1:0))>1e-8)
	printf("<v%d,v%d> = %e\n", i, j, temp);
    }
  }
#endif

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
  
  
  // Set up permutation array...
  allocate_vector(iterations, pindex);
  if (requestedStep == TensorStep3) {
    pindex[iterations-1] = 0;
    for (int i=0; i<iterations-1; i++)
      pindex[i] = i+1;
  }
  else {
    // No permutations needed
    for (int i=0; i<iterations; i++)
      pindex[i] = i;
  }

  
  // Solve linear system for Newton direction...
  /* ...but we can't update Newton vector in group yet... */
  yn = backsolve(hess, vecg, pindex, iterations);
  if (requestedStep == TensorStep2) {
    // Calculate yn = Q1*yn1 (without explicit Q1)
    double temp = 0;
    for (int i=k1; i>=0; i--) {
      temp += yn[i]/sqrt(qNorm2[i]);
      yn[i] = temp*qk[i];
      if (i>0)
	yn[i] += yn[i-1]*q1subdiag[i-1]/sqrt(qNorm2[i-1]);
    }
  }

  // Compute Newton direction...
  dNewton->init(0);
  for (int i=0; i<=k; i++)
    dNewton->update(-yn[i], *basisVptr[i], 1);
  if (precondition == Right) {
    *vecw = *dNewton;
    status = soln.applyRightPreconditioning(false, *localParamsPtr,
					    *vecw, *dNewton);
    if (status != NOX::Abstract::Group::Ok)
      NOX::Direction::Tensor::throwError("compute",
					 "Unable to apply preconditioner.");
  }
  if (!isDir0) dNewton->update(1, dir, 1);  // add in initial guess
  if (isAugmentSubspace) dNewton->update(-yn[k1]/normS, *scPtr, 1);

#if DEBUG_LEVEL > 1
  cout << "yn ";
  print_vector(iterations, yn);
  cout << "Newton Direction 1: ";
  dNewton->print();
#endif
#ifdef CHECK_RESIDUALS
  residual = getNormModelResidual(*dNewton, soln, solver, false);
  printf(" dNewton model residual = %8e\n", residual);
  fprime = getDirectionalDerivative(*dNewton, soln);
  printf(" dNewton directional derivative = %e\n", fprime);
#endif // CHECK_RESIDUALS

#ifdef CALCULATE_GRP_NEWTON
  /* The following code will compute the Newton step directly and
     store the result in the solution group.  Since at this time we
     can't update the Newton vector in the group with the result
     calculated above, this code is here to do just that (if ever
     necessary).  */
  
  // Compute the Newton direction directly
  status = soln.computeNewton(*localParamsPtr);
  
  // Computing the backup Newton step failed, but maybe it is okay to use...
  if ((status != NOX::Abstract::Group::Ok) && (doRescue == false)) {
    NOX::Direction::Tensor::throwError("compute",
				       "Unable to solve Newton system");
  }
  else if ((status != NOX::Abstract::Group::Ok) &&
	   (doRescue == true)) {
    if (utils.isPrintProcessAndType(NOX::Utils::Warning)) {
      double accuracy = soln.getNormNewtonSolveResidual();
      cout << "WARNING: NOX::Direction::Tensor::compute() - Linear solve "
	   << "failed to achieve convergence - using the step anyway "
	   << "since \"Rescue Bad Newton Solve\" is true. "
	   << "Using solution with accuracy of "
	   << utils.sciformat(accuracy) << "." << endl;
    }
  }
  
#if DEBUG_LEVEL > 1
  cout << "Newton Direction 2: ";
  soln.getNewton().print();
#endif
#ifdef CHECK_RESIDUALS
  residual = soln.getNormNewtonSolveResidual();
  printf(" dNewton* model residual (getNormNewton...) = %8e\n", residual);

  residual = getNormModelResidual(*dNewton, soln, solver, false);
  printf(" dNewton* model residual (manual) = %8e\n", residual);

  fprime = getDirectionalDerivative(soln.getNewton(), soln);
  printf(" dNewton* directional derivative = %e\n", fprime);

  //double stdn = basisVptr[0]->dot(soln.getNewton()) * normS;
  double stdn = scPtr->dot(soln.getNewton());
  printf(" sTdn (dNewton*) =  %e\n", stdn);
#endif
#endif  // CALCULATE_GRP_NEWTON


  if (requestedStep == TensorStep3  ||
      requestedStep == TensorStep2) {

    // Set the row number of quadratic equation
    int indxq = 0;
    if (requestedStep == TensorStep3)
      indxq = 0;
    else if (requestedStep == TensorStep2)
      indxq = k1;
  
    // Solve system for tensor direction...
    double beta2 = dir0xsc + normS * y1;
    if (requestedStep == TensorStep2)
      beta2 = dir0xsc + sqrt(qNorm2[k1]) * y1;
    beta2 = beta2 * beta2;
    for (int i=0; i<maxDim; i++) {
      vecz[i] = -vecg[i] - vecq[i] * beta2 - hess[i][indxq] * y1;
    }
    double* yt = backsolve(hess, vecz, pindex, iterations-1, iterations);
    yt[indxq] = y1;

    if (requestedStep == TensorStep2) {
      // Calculate yt = Q1*yt1 (without explicit Q1)
      double temp = 0;
      for (int i=k1; i>=0; i--) {
      temp += yt[i]/sqrt(qNorm2[i]);
      yt[i] = temp*qk[i];
      if (i>0)
	yt[i] += yt[i-1]*q1subdiag[i-1]/sqrt(qNorm2[i-1]);
      }
    }

    // Compute the tensor direction...
    dTensor->init(0);
    for (int i=0; i<=k; i++)
      dTensor->update(yt[i], *basisVptr[i], 1);
    if (precondition == Right) {
      *vecw = *dTensor;
      status = soln.applyRightPreconditioning(false, *localParamsPtr,
					      *vecw, *dTensor);
      if (status != NOX::Abstract::Group::Ok)
	NOX::Direction::Tensor::throwError("compute",
					   "Unable to apply preconditioner.");
    }
    if (!isDir0) dTensor->update(1, dir, 1);  // add in initial guess
    if (isAugmentSubspace) dTensor->update(yt[k1]/normS, *scPtr, 1);

    isTensorCalculated = true;
    isCLParamsCalculated = false;

#if DEBUG_LEVEL > 1
    cout << "yt ";
    print_vector(iterations, yt);
    cout << "Tensor Direction: "; 
    dTensor->print();
#endif
#ifdef CHECK_RESIDUALS
    printf(" ||sc||*y1+d0'*sc = %e\n", normS*y1+dir0xsc);

    double residual = getNormModelResidual(*dTensor, soln, solver, true);
    printf(" dTensor model residual = %8e\n", residual);
    fprime = getDirectionalDerivative(*dTensor, soln);
    printf(" dTensor directional derivative = %e\n", fprime);

    if ( lambdaBar == 1.0  &&  fabs((residual - terr)/terr) > 1e-3) 
      printf("  *** Warning - check residuals ***\n");
#endif  // CHECK_RESIDUALS

    // Compute parameters needed for curvilinear linesearch
    stJinvF = normS*qc/qb;
    stJinvA = 2*qa/qb/normS;
#if DEBUG_LEVEL > 1
    printf("stJinvF = %e  ", stJinvF);
    printf("stJinvA = %e\n", stJinvA);
#endif

    delete [] yt;
    
  } // endif (requestedStep == TensorStep3 || Step2)


#if DEBUG_LEVEL > 1
  cout << "Hessenberg: \n";
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
#endif
  

  if (requestedStep == TensorStepFP) {

    // Compute the tensor terms....
    // (note that ac is not multiplied by 2)
    NOX::Abstract::Vector* sPtr = soln.getX().clone(ShapeCopy);
    *sPtr = *scPtr;
    sPtr->scale(1/normS);
    NOX::Abstract::Vector* aPtr = soln.getF().clone(ShapeCopy);
    *aPtr = *acPtr;
    aPtr->scale(normS*normS);
    
    if (precondition == Left) {
      *vecw = *aPtr;
      printf("  left preconditioning FP\n");
      status = soln.applyRightPreconditioning(false, *localParamsPtr, *vecw, *aPtr);
      if (status != NOX::Abstract::Group::Ok)
	NOX::Direction::Tensor::throwError("compute",
					   "Unable to apply preconditioner.");
    }
    if (precondition == Right) {
      *vecw = *sPtr;
      printf("  right preconditioning FP\n");
      status = soln.applyRightPreconditioning(true, *localParamsPtr, *vecw, *sPtr);
      if (status != NOX::Abstract::Group::Ok)
	NOX::Direction::Tensor::throwError("compute",
					   "Unable to apply preconditioner.");
    }


    // Solve system for Feng-Pulliam tensor direction...
    if (isDir0) {

      printf(" dir is zero\n");

    dir0_in_subspace_Vm:  // Label in case initial guess d0 is in span{Vm}
      
      // Compute some vectors...
      for (int i=0; i<iterations+1; i++) 
	vecz[i] = aPtr->dot(*basisVptr[i]);
      for (int i=0; i<iterations; i++)
	applyGivens(givensC[0][i], givensS[0][i], vecz[i], vecz[i+1]);
      double* yt = backsolve(hess, vecz, pindex, iterations);
      for (int i=0; i<iterations; i++) 
	vecz[i] = sPtr->dot(*basisVptr[i]);
    
      // Form beta equation...
      double aa = inner_product(iterations, vecz, yt);
      double bb = 1.0;
      double cc = inner_product(iterations, vecz, yn) - dir0xsc/normS;
                                                      //sPtr->dot(dir);
      double qval = 0;
      double beta = calculateBeta(aa, bb, cc, qval, lambdaBar);

      // Calculate basis coefficients...
      for (int i=0; i<iterations; i++) 
	yt[i] = -yn[i] - yt[i]*beta*beta;

      if (qval != 0.0) {
	double* vect1 = backsolve(hess, vecz, pindex, iterations, 0, true);
	double* vect2 = backsolve(hess, vect1, pindex, iterations);
	double tau = qval / inner_product(iterations, vect2, vecz);
	// bwb - question on sign of tau (possible error in F-P paper?)
	
	for (int i=0; i<iterations; i++) 
	  yt[i] += tau*vect2[i];

	delete [] vect1;
	delete [] vect2;
      }

      // Compute the tensor step by summing the basis vectors...
      dTensor->init(0);
      for (int i=0; i<iterations; i++)
	dTensor->update(yt[i], *basisVptr[i], 1);
      if (precondition == Right) {
	*vecw = *dTensor;
	status = soln.applyRightPreconditioning(false, *localParamsPtr,
						*vecw, *dTensor);
	if (status != NOX::Abstract::Group::Ok)
	  NOX::Direction::Tensor::throwError("compute",
					   "Unable to apply preconditioner.");
      }
      if (!isDir0) dTensor->update(1, dir, 1);  // add in initial guess

      delete [] yt;
    }
    else {   // nonzero inital guess

      // Local variables
      double* vectmp = NULL;
      double* vecf = NULL;
      double* yt = NULL;
      
      // Compute unit vectors w = Fc-r0 = J*d0 
      soln.applyJacobian(dir, *vecw);
      if (precondition == Left) {
	*dTensor = *vecw;  // temporary vector
	status = soln.applyRightPreconditioning(false, *localParamsPtr, 
						*dTensor,
						*vecw);
	if (status != NOX::Abstract::Group::Ok)
	  NOX::Direction::Tensor::throwError("compute",
			     "Unable to apply preconditioner.");
      }
      double normU = vecw->norm();
      vecw->scale(1/normU);
      
      // Construct right-most column...
      for (int i=0; i<iterations+1; i++) 
	vecz[i] = vecw->dot(*basisVptr[i]);
      for (int i=0; i<iterations; i++)
	applyGivens(givensC[0][i], givensS[0][i], vecz[i], vecz[i+1]);
      double eta = 1.0 - inner_product(iterations, vecz, vecz);
      printf(" eta = %e\n", eta);
      //eta = 1e-16;   // for testing purposes
      if (eta<=0) goto dir0_in_subspace_Vm;
      double gamma = sqrt(eta);

      // Construct factored form of Jhat'*Jhat...
      double** matC = NULL;
      allocate_matrix(iterations+1, iterations+1, matC);
      for (int i=0; i<iterations; i++) {
	for (int j=0; j<iterations; j++) 
	  matC[i][j] = hess[i][j];
	matC[i][iterations] = vecz[i];
      }
      for (int j=0; j<iterations; j++) 
	matC[iterations][j] = 0;
      matC[iterations][iterations] = gamma;
            
      NOX::Abstract::Vector* bPtr = soln.getF().clone(ShapeCopy);
      *bPtr = soln.getF();
      if (precondition == Left) {
	*dTensor = *bPtr;  // temporary vector
	status = soln.applyRightPreconditioning(false, *localParamsPtr, 
						*dTensor, *bPtr);
	if (status != NOX::Abstract::Group::Ok)
	  NOX::Direction::Tensor::throwError("compute",
			     "Unable to apply preconditioner.");
      }

      // Reset permutation array to be one larger
      delete pindex;  pindex = NULL;
      allocate_vector(iterations+1, pindex);
      for (int i=0; i<iterations+1; i++)
	pindex[i] = i;

      // Calculate tp = C\C'\(inv(Q)*(V'*b))
      for (int i=0; i<iterations+1; i++)          // tp = V(:,1:i+1)'*(-b)
	vecz[i] = basisVptr[i]->dot(*bPtr);
      for (int i=0; i<iterations; i++)            // tp = inv(Q)*tp
	applyGivens(givensC[0][i], givensS[0][i], vecz[i], vecz[i+1]);
      for (int i=iterations-1; i>=0; i--) {       // tp = [R'*tp(1:i);u'*(-b)
	vecz[i] = vecz[i]*hess[i][i];
	for (int j=i-1; j>=0; j--)
	  vecz[i] += vecz[j]*hess[j][i];
      }
      vecz[iterations] = vecw->dot(*bPtr);
      vectmp = backsolve(matC, vecz, pindex, iterations+1, 0, true);
      vecf   = backsolve(matC, vectmp, pindex, iterations+1);
      delete vectmp;  vectmp = NULL;
      
      // Calculate tp = C\C'\(inv(Q)*(V'*ac))
      for (int i=0; i<iterations+1; i++)          // tp = V(:,1:i+1)'*(-b)
	vecz[i] = basisVptr[i]->dot(*aPtr);
      for (int i=0; i<iterations; i++)            // tp = inv(Q)*tp
	applyGivens(givensC[0][i], givensS[0][i], vecz[i], vecz[i+1]);
      for (int i=iterations-1; i>=0; i--) {       // tp = [R'*tp(1:i);u'*(-b)
	vecz[i] = vecz[i]*hess[i][i];
	for (int j=i-1; j>=0; j--)
	  vecz[i] += vecz[j]*hess[j][i];
      }
      vecz[iterations] = vecw->dot(*aPtr);
      vectmp = backsolve(matC, vecz, pindex, iterations+1, 0, true);
      yt     = backsolve(matC, vectmp, pindex, iterations+1);
      delete vectmp;  vectmp = NULL;

      // Calculate sb
      for (int i=0; i<iterations; i++) 
	vecz[i] = sPtr->dot(*basisVptr[i]);
      //vecz[iterations] = sPtr->dot(dir)/normU; // bwb - doesn't work w/ right
      vecz[iterations] = dir0xsc/normS/normU;    // bwb - works with right prec

      // Form beta equation...
      double aa = inner_product(iterations+1, vecz, yt);
      double bb = 1.0;
      double cc = inner_product(iterations+1, vecz, vecf);
      double qval = 0;
      double beta = calculateBeta(aa, bb, cc, qval, lambdaBar);

      // Calculate basis coefficients...
      for (int i=0; i<iterations+1; i++) 
	yt[i] = -vecf[i] - yt[i]*beta*beta;
      if (qval != 0.0) {
	double* vect1 = backsolve(matC, vecz, pindex, iterations+1, 0, true);
	double* vect2 = backsolve(matC, vect1, pindex, iterations+1);
	double tau = qval / inner_product(iterations+1, vect2, vecz);

	printf(" compare values2:  %e  %e\n",
	       inner_product(iterations+1, vect1, vect1),
	       inner_product(iterations+1, vect2, vecz));
	
	for (int i=0; i<iterations+1; i++) 
	  yt[i] += tau*vect2[i];

	delete_matrix(matC);
	delete [] vect1;
	delete [] vect2;
      }

      // Compute the tensor step by summing the basis vectors...
      dTensor->init(0);
      for (int i=0; i<iterations; i++)
	dTensor->update(yt[i], *basisVptr[i], 1);
      if (precondition == Right) {
	*vecw = *dTensor;
	status = soln.applyRightPreconditioning(false, *localParamsPtr,
						*vecw, *dTensor);
	if (status != NOX::Abstract::Group::Ok)
	  NOX::Direction::Tensor::throwError("compute",
					   "Unable to apply preconditioner.");
      }
      dTensor->update(yt[iterations]/normU, dir, 1);  // add in initial guess

      delete [] yt;
      delete [] vecf;
      delete bPtr;
    }  // endif (isDir0)
    
    
#if DEBUG_LEVEL > 1
    cout << "yt ";
    print_vector(iterations, yt);
    cout << "yn ";
    print_vector(iterations, yn);
    printf("qbeta equation: %e  %e  %e\n", aa, bb, cc);
    printf("qbeta results : %e  %e  %e\n", beta, qval, lambdaBar);
#endif
#ifdef CHECK_RESIDUALS
    residual = getNormModelResidual(*dTensor, soln, solver, true);
    printf(" dTensor (F-P) model residual = %8e\n", residual);
    fprime = getDirectionalDerivative(*dTensor, soln);
    printf(" dTensor directional derivative = %e\n", fprime);
#endif  // CHECK_RESIDUALS

    delete sPtr;
    delete aPtr;
    
  } // endif (requestedStep == TensorStepFP)


  // Set search direction...
  if (requestedStep == NewtonStep)
    dir = *dNewton;
  else
    dir = *dTensor;


  // Cleanup memory
  delete [] pindex;
  delete [] yn;

 Cleanup_before_exit:

#ifdef STORE_HESSENBERG
  delete_matrix(hess2);
#endif

  return true;
}


bool NOX::Direction::Tensor::computeCurvilinearStep(NOX::Abstract::Vector& dir,
				    NOX::Abstract::Group& soln,
				    const NOX::Solver::Generic& solver, 
				    double lambda)
{
  bool ok = true;
  NOX::Abstract::Group::ReturnType status;

  // If the tensor step has not been computed, then return the Newton
  // step because other necessary vectors have not been computed.
  // This is most certainly the case for the first nonlinear
  // iteration, where the Newton step is used.
  if (!isTensorCalculated) {
    dir = *dNewton;
    dir.scale(lambda);
    return ok;
  }

  // If the curvilinear step dt(lambdaBar) and the 2 parameters stJinvF
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
	status = soln.applyRightPreconditioning(false, *localParamsPtr,
						*vecw, *dTLambda);
	if (status != NOX::Abstract::Group::Ok)
	  NOX::Direction::Tensor::throwError("compute",
					     "Unable to apply preconditioner.");
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

#if DEBUG_LEVEL > 0
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
    double largeRoot = 0;
    
    // Calculate beta(lambda)
    discriminant = 1 - 2*stJinvF*stJinvA*lambda;
    if (discriminant < 0) {
      cout << "Warning: discriminant is negative ("<< discriminant << ").\n";
      discriminant = 0;
    }
    largeRoot = (-1 - sqrt(discriminant)) / (2*stJinvA);
    double beta1 = (lambda*stJinvF) / (largeRoot*stJinvA);

    // Calculate beta(lambdaBar)
    discriminant = 1 - 2*stJinvF*stJinvA*lambdaBar;
    if (discriminant < 0) {
      // Probably just slightly off due to round-off error, so set to zero
      discriminant = 0;
    }
    largeRoot = (-1 - sqrt(discriminant)) / (2*stJinvA);
    double beta2 = (lambdaBar*stJinvF) / (largeRoot*stJinvA);

    // Set ratio = beta(lambda)/beta(lambdaBar)
    betaRatio = beta1/beta2;
  }

  // Now recast the coefficients to work with the vectors dN, dT, dTlambda
  double coefNewt = lambda - betaRatio*betaRatio*lambdaBar - 
    alpha*(1-lambdaBar);
  double coefTlam = betaRatio*betaRatio - alpha;
  double coefTens = alpha;

#if DEBUG_LEVEL > 0
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


const NOX::Abstract::Vector& NOX::Direction::Tensor::getNewton() const 
{
  return *dNewton;
}


// private

void NOX::Direction::Tensor::throwError(const string& functionName,
					const string& errorMsg)
{
  if (utils.isPrintProcessAndType(NOX::Utils::Error))
    cerr << "NOX::Direction::Tensor::" << functionName << " - " << errorMsg <<
      endl;
  throw "NOX Error";
}


void** NOX::Direction::Tensor::allocate_matrix(int rows, int cols, double**& a)
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


void* NOX::Direction::Tensor::allocate_vector(int n, int*& x)
{
  if (x) {
    // delete x;
    cout << "Warning: Possibly a previously allocated vector\n";
  }

  x = new int [n];
  if (!x) printf("Allocation error with malloc()\n");
  return (void*) x;
}


void* NOX::Direction::Tensor::allocate_vector(int n, double*& x)
{
  if (x) {
    // delete x;
    cout << "Warning: Possibly a previously allocated vector\n";
  }

  x = new double [n];
  if (!x) printf("Memory allocation error in allocate_vector().\n");
  return (void*) x;
}


void NOX::Direction::Tensor::delete_matrix(double** A)
   /* Free the memory previously allocated in allocate_matrix()  */
{
  if (A != NULL) {
    delete [] A[0];
  }
  delete [] A;
}


void NOX::Direction::Tensor::print_matrix(int rows, int cols, double** A)
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


void NOX::Direction::Tensor::print_vector(int n, double* x)
{
  int i;

  printf("Vector = \n");
  for (i=0; i<n; i++)
    printf("   %16.8e\n", x[i]);
  printf("\n");
}


void NOX::Direction::Tensor::print_vector(int n, int* x)
{
  int i;

  printf("Vector = \n");
  for (i=0; i<n; i++)
    printf("   %d\n", x[i]);
  printf("\n");
}


double NOX::Direction::Tensor::inner_product(int n, double* x, double* y)
{
  int i;
  double sum = 0.0;

  for (i=0; i<n; i++)
    {
      sum += x[i]*y[i];
    }

  return sum;
}


double NOX::Direction::Tensor::norm(int n, double* x)
{
  double temp = inner_product(n,x,x);
  return sqrt(temp);
}


void NOX::Direction::Tensor::computeGivens(double a, double b,
					   double& c, double& s)
  /* This function takes the scalars a,b and computes c = cos(theta)
  ** and s = sin(theta) for a Givens rotation matrix G = [c s; -s c].
  ** Premultiplication by G' amounts to a counterclockwise rotation of
  ** theta radians such that [c s; -s c]^T*[a; b] = [r; 0].
  */
{
  double tmp;

  /*    
  tmp = sqrt(a*a + b*b);
  c =  a/tmp;
  s = -b/tmp;
  return;
  */

  if (b == 0) {
    c = 1;
    s = 0;
  }
  else {
    if (fabs(b) > fabs(a)) {
      tmp = -a/b;
      s = 1 / sqrt(1 + tmp*tmp);
      c = s*tmp;
    }
    else {
      tmp = -b/a;
      c = 1 / sqrt(1 + tmp*tmp);
      s = c*tmp;
    }
  }
}
  

void NOX::Direction::Tensor::applyGivens(double c, double s,
					 double& a, double& b)
  /* This function will apply a Givens rotation matrix G = [c s; -s c]
  ** to the vector v = [a; b] and compute the rotation G'*v = inv(G)*v
  */
{
  double tmp;
  
  tmp = c*a - s*b;
  b   = s*a + c*b;
  a   = tmp;
}


void NOX::Direction::Tensor::computeHouseholder(double** a, int j, 
						int i1, int i2, double* z)
  /* This function takes the jth column of matrix a and computes
  ** the householder vector over the indices i1:i2 and returns it in the
  ** vector z.  Also, the jth column of matrix is zeroed (as if the
  ** Householder reflection were applied to this column).
  */
{
  // Compute reflection
  for (int i=0; i<i2-i1+1; i++)
    z[i] = a[i+i1][j];
  double tmp = norm(i2-i1+1, z);
  if (z[0] < 0) tmp = -tmp;
  z[0] += tmp;

  // Apply reflection
  a[i1][j] = -tmp;
  for (int i=i1+1; i<=i2; i++)
    a[i][j] = 0;
}
  

void NOX::Direction::Tensor::applyHouseholder(const double* z, double** a,
					      int j, int i1, int i2)
  /* Apply Householder reflection to the jth column of matrix
  ** a(i1:i2,k).  The rank-one reflector matrix is stored in z,
  ** which is not normalized.  The result is returned in a.  Note that
  ** the indices i1, i2, and k are all 0-based
  */
{
  double ztz = 0;
  for (int i=0; i<i2-i1+1; i++)
    ztz += z[i]*z[i];
  
  double tmp = 0;
  for (int i=i1; i<=i2; i++)
    tmp += z[i-i1]*a[i][j];
  tmp /= ztz;
    
  for (int i=i1; i<=i2; i++)
    a[i][j] -= 2*z[i-i1]*tmp;
}



void NOX::Direction::Tensor::applyHouseholder(const double* z, double* a,
					      int i1, int i2)
  /* Apply Householder reflection to the vector a(i1:i2).
  ** The rank-one reflector matrix is stored in z, which is not
  ** normalized.  The result is returned in a.  Note that the indices i1 and
  ** i2 are 0-based.
  */
{
  double ztz = 0;
  for (int i=0; i<i2-i1+1; i++)
    ztz += z[i]*z[i];
  
  double tmp = 0;
  for (int i=i1; i<=i2; i++)
    tmp += z[i-i1]*a[i];
  tmp /= ztz;

  for (int i=i1; i<=i2; i++)
    a[i] -= 2*z[i-i1]*tmp;
}


double NOX::Direction::Tensor::calculateBeta(double qa, double qb, double qc,
					     double& qval, double& lambdaBar,
					     double dir0xsc,
					     double normS)
{
  double beta = 0;
  double discriminant = qb*qb - 4*qa*qc;

  if (discriminant < 0) {

    // no real root
    beta = -qb/qa/2;
    qval = qa*beta*beta + qb*beta + qc;
    lambdaBar = qb*qb / (4*qa*qc);
#if DEBUG_LEVEL > 0
    cout << "  ####  LambdaBar = " << lambdaBar << "  ####\n";
#endif
  }
  else {
    qval = 0;
    lambdaBar = 1.0;
    if (fabs(qa/qb) < 1e-8) {
#if DEBUG_LEVEL > 0
      cout << "qa is relatively small\n";
#endif 
      beta = -qc/qb;
    }
    else {
      double tmp1 = (-qb + sqrt(discriminant)) / (2*qa);
      double tmp2 = (-qb - sqrt(discriminant)) / (2*qa);
      beta = (fabs(dir0xsc + normS*tmp1) < fabs(dir0xsc + normS*tmp2)) ?
	tmp1 : tmp2;
    }
  }
  return beta;
}


double NOX::Direction::Tensor::getNormModelResidual(
				       const NOX::Abstract::Vector& dir,
				       const NOX::Abstract::Group& soln,
				       const NOX::Solver::Generic& solver,
				       bool isTensorModel)
{
  NOX::Abstract::Vector* residualPtr = NULL;

  // Compute residual of Newton model...
  residualPtr = soln.getF().clone(ShapeCopy);
  soln.applyJacobian(dir, *residualPtr);    
  residualPtr->update(1.0, soln.getF(), 1.0);

  if (isTensorModel) {

    // Compute the tensor term ac....
    Abstract::Vector* sPtr = soln.getX().clone(ShapeCopy);
    Abstract::Vector* aPtr = soln.getF().clone(ShapeCopy);
    *sPtr = soln.getX();
    sPtr->update(1.0, solver.getPreviousSolutionGroup().getX(), -1.0);
    double normS = sPtr->norm();
    soln.applyJacobian(*sPtr, *aPtr);
    aPtr->update(1.0, solver.getPreviousSolutionGroup().getF(), -1.0);
    aPtr->update(-1.0, soln.getF(), 1.0);
    aPtr->scale(1/(normS*normS));
    sPtr->scale(1/normS);

    double beta = sPtr->dot(dir);
    printf(" sc'*dt           = %e\n", beta*normS);
    printf(" norm(dt) = %e\n", dir.norm());
    residualPtr->update(beta*beta, *aPtr, 1.0);

    delete aPtr;
    delete sPtr;
  }
  
  if (precondition == Left) {
    NOX::Abstract::Vector* tmpPtr = soln.getF().clone(ShapeCopy);
    *tmpPtr = *residualPtr;
    NOX::Abstract::Group::ReturnType status =
      soln.applyRightPreconditioning(false, *localParamsPtr,
				     *tmpPtr, *residualPtr);
    if (status != NOX::Abstract::Group::Ok)
      NOX::Direction::Tensor::throwError("compute",
					 "Unable to apply preconditioner.");
    delete tmpPtr;
  }

  double modelNorm = residualPtr->norm();
  delete residualPtr;
  return modelNorm;
}


double NOX::Direction::Tensor::getDirectionalDerivative(
				       const NOX::Abstract::Vector& dir,
				       const NOX::Abstract::Group& soln)
{
  NOX::Abstract::Vector* tmpPtr = soln.getF().clone(ShapeCopy);
  soln.applyJacobian(dir,*tmpPtr);
  double fprime = tmpPtr->dot(soln.getF());
  delete tmpPtr;
  return fprime;
}


double* NOX::Direction::Tensor::backsolve(double** U, double* b, int* perm,
					  int n, int dim,
					  bool isTranspose)
     /* This function solves the triangular system Ux=b when provided
      * an upper triangular matrix U. The array "perm" is a
      * permutation array.  The pointer returned is the newly created
      * solution vector x.  The variable "dim" is the maximum
      * dimension (if nonzero) of the solution vector, which may be
      * different due to the indices in the permutation array.  If the
      * parameter "isTranspose" is true, then the triangular system
      * U'x=b is solved for x.
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

  if (isTranspose) {
    //  Solve U'x=y with forward substitution
    for (i=0; i<n; i++) {
      pi = perm[i];
      temp = b[pi];
      for (j=0; j<i; j++) {
	pj = perm[j];
	temp -= U[pj][pi] * x[pj];
      }
      if (U[pi][pi] == 0) printf("U is singular\n");
      x[pi] = temp / U[pi][pi];
    }
  }
  else {
    //  Solve Ux=y with backward substitution
    for (i=n-1; i>=0; i--) {
      pi = perm[i];
      temp = b[pi];
      for (j=i+1; j<n; j++) {
	pj = perm[j];
	temp -= U[pi][pj] * x[pj];
      }
      if (U[pi][pi] == 0) printf("U is singular\n");
      x[pi] = temp / U[pi][pi];
    }  
  }
  
  return x;
}

#endif  // WITH_PRELEASE
