// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
//
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Andy Salinger (agsalin@sandia.gov) or Eric Phipps
// (etphipp@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#include "LOCA.H"
#include "LOCA_Epetra.H"

#include "LOCA_MultiContinuation_CompositeConstraint.H"
#include "NOX_TestCompare.H"

// Trilinos Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"

// User's application specific files 
#include "LinearConstraint.H"

int main(int argc, char *argv[])
{
  int n = 10;
  int ierr = 0;
  double reltol = 1.0e-14;
  double abstol = 1.0e-14;
  int MyPID;

  try {

    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init(&argc,&argv);
#endif

    // Create a communicator for Epetra objects
#ifdef HAVE_MPI
    Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
    Epetra_SerialComm Comm;
#endif

    MyPID = Comm.MyPID();

    // Create the map
    Epetra_Map map(n, 0, Comm);
    
    bool verbose = false;
    // Check for verbose output
    if (argc>1)
      if (argv[1][0]=='-' && argv[1][1]=='v') 
	verbose = true;

    // Create and initialize the parameter vector
    LOCA::ParameterVector pVector;
    pVector.addParameter("Param 1",  1.69);
    pVector.addParameter("Param 2", -9.7);
    pVector.addParameter("Param 3",  0.35);
    pVector.addParameter("Param 4", -0.78);
    pVector.addParameter("Param 5",  2.53);

    // Create parameter list
    Teuchos::RefCountPtr<Teuchos::ParameterList> paramList =
      Teuchos::rcp(new Teuchos::ParameterList);

    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
    Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.set("MyPID", MyPID);
    if (verbose)
       nlPrintParams.set("Output Information", 
				  NOX::Utils::Error +
				  NOX::Utils::Details +
				  NOX::Utils::OuterIteration + 
				  NOX::Utils::InnerIteration + 
				  NOX::Utils::Warning +
				  NOX::Utils::TestDetails + 
				  NOX::Utils::StepperIteration +
				  NOX::Utils::StepperDetails);
     else
       nlPrintParams.set("Output Information", NOX::Utils::Error);

    // Create global data object
    Teuchos::RefCountPtr<LOCA::GlobalData> globalData =
      LOCA::createGlobalData(paramList);

    Epetra_Vector clone_vec(map);
    NOX::Epetra::Vector nox_clone_vec(clone_vec);

    Teuchos::RefCountPtr<NOX::Abstract::Vector> x = 
      nox_clone_vec.clone(NOX::ShapeCopy);
    x->random();

    Teuchos::RefCountPtr<NOX::Abstract::MultiVector> dx1 = 
      nox_clone_vec.createMultiVector(3);
    Teuchos::RefCountPtr<NOX::Abstract::MultiVector> dx2 = 
      nox_clone_vec.createMultiVector(1);
    Teuchos::RefCountPtr<NOX::Abstract::MultiVector> dx3 = 
      nox_clone_vec.createMultiVector(2);
    Teuchos::RefCountPtr<NOX::Abstract::MultiVector> dx4 = 
      nox_clone_vec.createMultiVector(2);
    dx1->random();
    dx2->random();
    dx3->init(0.0);
    dx4->random();

    Teuchos::RefCountPtr<NOX::Abstract::MultiVector> dx_all = 
      dx1->clone(NOX::DeepCopy);
    dx_all->augment(*dx2);
    dx_all->augment(*dx3);
    dx_all->augment(*dx4);

    NOX::Abstract::MultiVector::DenseMatrix dp1(dx1->numVectors(),
						pVector.length());
    NOX::Abstract::MultiVector::DenseMatrix dp2(dx2->numVectors(),
						pVector.length());
    NOX::Abstract::MultiVector::DenseMatrix dp3(dx3->numVectors(),
						pVector.length());
    NOX::Abstract::MultiVector::DenseMatrix dp4(dx4->numVectors(),
						pVector.length());
    dp1.random();
    dp2.random();
    dp3.random();
    dp4.random();

    NOX::Abstract::MultiVector::DenseMatrix dp_all(dx_all->numVectors(),
						   pVector.length());
    for (int j=0; j<dp_all.numCols(); j++) {
      for (int i=0; i<dp1.numRows(); i++)
	dp_all(i,j) = dp1(i,j);
      for (int i=0; i<dp2.numRows(); i++)
	dp_all(dp1.numRows()+i,j) = dp2(i,j);
      for (int i=0; i<dp3.numRows(); i++)
	dp_all(dp1.numRows()+dp2.numRows()+i,j) = dp3(i,j);
      for (int i=0; i<dp4.numRows(); i++)
	dp_all(dp1.numRows()+dp2.numRows()+dp3.numRows()+i,j) = dp4(i,j);
    }


    vector< Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface> > constraintObjs(4);
    Teuchos::RefCountPtr<LinearConstraint> linear_constraint;

    linear_constraint = Teuchos::rcp(new LinearConstraint(dx1->numVectors(),
							  pVector,
							  nox_clone_vec));
    linear_constraint->setDgDx(*dx1);
    linear_constraint->setDgDp(dp1);
    linear_constraint->setIsZeroDX(false);
    constraintObjs[0] = linear_constraint;

    linear_constraint = Teuchos::rcp(new LinearConstraint(dx2->numVectors(),
							  pVector,
							  nox_clone_vec));
    linear_constraint->setDgDx(*dx2);
    linear_constraint->setDgDp(dp2);
    linear_constraint->setIsZeroDX(false);
    constraintObjs[1] = linear_constraint;

    linear_constraint = Teuchos::rcp(new LinearConstraint(dx3->numVectors(),
							  pVector,
							  nox_clone_vec));
    linear_constraint->setDgDx(*dx3);
    linear_constraint->setDgDp(dp3);
    linear_constraint->setIsZeroDX(true);
    constraintObjs[2] = linear_constraint;

    linear_constraint = Teuchos::rcp(new LinearConstraint(dx4->numVectors(),
							  pVector,
							  nox_clone_vec));
    linear_constraint->setDgDx(*dx4);
    linear_constraint->setDgDp(dp4);
    linear_constraint->setIsZeroDX(false);
    constraintObjs[3] = linear_constraint;

    // Check some statistics on the solution
    NOX::TestCompare testCompare(globalData->locaUtils->out(), 
				 *(globalData->locaUtils));

    LOCA::MultiContinuation::CompositeConstraint composite(globalData,
							   constraintObjs);
    composite.setX(*x);

    LinearConstraint combined(dx_all->numVectors(), pVector, nox_clone_vec);
    combined.setDgDx(*dx_all);
    combined.setDgDp(dp_all);
    combined.setX(*x);

    //
    // test computeConstraints()
    //

    composite.computeConstraints();
    combined.computeConstraints();

    int numConstraints = dx_all->numVectors();
    const NOX::Abstract::MultiVector::DenseMatrix& g_composite =
      composite.getConstraints();
    const NOX::Abstract::MultiVector::DenseMatrix& g_combined =
      combined.getConstraints();

    ierr += testCompare.testMatrix(
				 g_composite, g_combined, reltol, abstol,
				 "CompositeConstraint::computeConstraints()");

    //
    // test computeDP()
    //

    vector<int> paramIDs(3);
    paramIDs[0] = 1;
    paramIDs[1] = 2;
    paramIDs[2] = 4;
    NOX::Abstract::MultiVector::DenseMatrix dgdp_composite(
							numConstraints,
							paramIDs.size()+1);
    NOX::Abstract::MultiVector::DenseMatrix dgdp_combined(
							numConstraints,
							paramIDs.size()+1);
    dgdp_composite.putScalar(0.0);
    dgdp_combined.putScalar(0.0);
    composite.computeDP(paramIDs, dgdp_composite, false);
    combined.computeDP(paramIDs, dgdp_combined, false);

    ierr += testCompare.testMatrix(
				 dgdp_composite, dgdp_combined, reltol, abstol,
				 "CompositeConstraint::computeDP()");

    //
    // test multiplyDX()
    //

    composite.computeDX();
    combined.computeDX();

    int numMultiply = 5;
    Teuchos::RefCountPtr<NOX::Abstract::MultiVector> A = 
      nox_clone_vec.createMultiVector(numMultiply);
    A->random();
    NOX::Abstract::MultiVector::DenseMatrix composite_multiply(numConstraints,
							       numMultiply);
    NOX::Abstract::MultiVector::DenseMatrix combined_multiply(numConstraints,
							      numMultiply);
    composite.multiplyDX(2.65, *A, composite_multiply);
    combined.multiplyDX(2.65, *A, combined_multiply);

    ierr += testCompare.testMatrix(composite_multiply, combined_multiply, 
				    reltol, abstol,
				    "CompositeConstraint::multiplyDX()");

    //
    // test addDX() (No Trans)
    //

    int numAdd = 5;
    NOX::Abstract::MultiVector::DenseMatrix B1(numConstraints, numAdd);
    B1.random();
    NOX::Abstract::MultiVector::DenseMatrix B2(numAdd, numConstraints);
    B2.random();

    Teuchos::RefCountPtr<NOX::Abstract::MultiVector> composite_add1 = 
      nox_clone_vec.createMultiVector(numAdd);
    composite_add1->random();
    Teuchos::RefCountPtr<NOX::Abstract::MultiVector> composite_add2 = 
      nox_clone_vec.createMultiVector(numAdd);
    composite_add2->random();

    Teuchos::RefCountPtr<NOX::Abstract::MultiVector> combined_add1 = 
      composite_add1->clone(NOX::DeepCopy);
    Teuchos::RefCountPtr<NOX::Abstract::MultiVector> combined_add2 = 
      composite_add2->clone(NOX::DeepCopy);

    composite.addDX(Teuchos::NO_TRANS, 1.45, B1, 2.78, *composite_add1);
    combined.addDX(Teuchos::NO_TRANS, 1.45, B1, 2.78, *combined_add1);

    ierr += testCompare.testMultiVector(
			       *composite_add1, *combined_add1, 
			       reltol, abstol, 
			       "CompositeConstraint::addDX() (No Trans)"); 

    // 
    // test addDX() (Trans)
    //

    composite.addDX(Teuchos::TRANS, 1.45, B2, 2.78, *composite_add2);
    combined.addDX(Teuchos::TRANS, 1.45, B2, 2.78, *combined_add2);

    ierr += testCompare.testMultiVector(
			       *composite_add2, *combined_add2, 
			       reltol, abstol, 
			       "CompositeConstraint::addDX() (Trans)");

    destroyGlobalData(globalData);
  }

  catch (std::exception& e) {
    cout << e.what() << endl;
    ierr = 1;
  }
  catch (const char *s) {
    cout << s << endl;
    ierr = 1;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
    ierr = 1;
  }

  if (MyPID == 0) {
    if (ierr == 0)
      cout << "All tests passed!" << endl;
    else
      cout << ierr << " test(s) failed!" << endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return ierr;
}
