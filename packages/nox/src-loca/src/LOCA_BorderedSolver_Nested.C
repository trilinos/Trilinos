// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "LOCA_BorderedSolver_Nested.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Factory.H"
#include "LOCA_MultiContinuation_ConstraintInterface.H"
#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H"
#include "LOCA_BorderedSystem_AbstractGroup.H"
#include "Teuchos_ParameterList.hpp"

LOCA::BorderedSolver::Nested::Nested(
	 const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& slvrParams): 
  globalData(global_data),
  solverParams(slvrParams),
  solver(),
  grp(),
  unbordered_grp(),
  myWidth(0),
  underlyingWidth(0),
  numConstraints(0)
{
  // Get "Nested Solver" sublist
  Teuchos::RefCountPtr<Teuchos::ParameterList> nestedSolverList = 
    Teuchos::rcp(&(solverParams->sublist("Nested Bordered Solver")),false);
  
  // Instantiate underlying solver
  solver = 
    globalData->locaFactory->createBorderedSolverStrategy(topParams,
							  nestedSolverList);
}

LOCA::BorderedSolver::Nested::~Nested()
{
}

void
LOCA::BorderedSolver::Nested::setMatrixBlocks(
         const Teuchos::RefCountPtr<const NOX::Abstract::Group>& group,
	 const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector>& blockA,
	 const Teuchos::RefCountPtr<const LOCA::MultiContinuation::ConstraintInterface>& blockB,
	 const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector::DenseMatrix>& blockC)
{
  string callingFunction = 
    "LOCA::BorderedSolver::Nested::setMatrixBlocks()";

  // Cast group to a bordered group
  grp = 
    Teuchos::rcp_dynamic_cast<const LOCA::BorderedSystem::AbstractGroup>(group);
  if (grp == Teuchos::null) 
    globalData->locaErrorCheck->throwError(
      callingFunction,
      string("Group must implement the LOCA::BorderedSystem::AbstractGroup")
      + string(" interface in order to use nested bordered solver strategy."));

  Teuchos::RefCountPtr<const LOCA::MultiContinuation::ConstraintInterfaceMVDX> con_mvdx = Teuchos::rcp_dynamic_cast<const LOCA::MultiContinuation::ConstraintInterfaceMVDX>(blockB);
  if (con_mvdx == Teuchos::null)
    globalData->locaErrorCheck->throwError(
		 callingFunction,
		 "Constraints object must be of type ConstraintInterfaceMVDX");

  bool isZeroA = (blockA.get() == NULL);
  bool isZeroB = con_mvdx->isDXZero();
  bool isZeroC = (blockC.get() == NULL);
  Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> blockB_dx;
  if (!isZeroB)
    blockB_dx = Teuchos::rcp(con_mvdx->getDX(), false);

  // ensure blocks B and C are not both zero
  if (isZeroB && isZeroC) 
    globalData->locaErrorCheck->throwError(
				        callingFunction,
				        "Blocks B and C cannot both be zero");

  // ensure blocks A and C are not both zero
  if (isZeroA && isZeroC) 
    globalData->locaErrorCheck->throwError(
				         callingFunction,
				         "Blocks A and C cannot both be zero");

  unbordered_grp = grp->getUnborderedGroup();

  // get number of outer constraints
  if (isZeroB)
    numConstraints = blockC->numRows();
  else
    numConstraints = blockB_dx->numVectors();

  // Get total bordered width
  underlyingWidth = grp->getBorderedWidth();
  myWidth = underlyingWidth + numConstraints;
  
  // combine blocks
  bool isCombinedAZero = grp->isCombinedAZero();
  bool isCombinedBZero = grp->isCombinedBZero();
  bool isCombinedCZero = grp->isCombinedCZero();
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> A;
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> B;
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> C;
  
  if (!isCombinedAZero || !isZeroA) {
    A = unbordered_grp->getX().createMultiVector(myWidth);
    A->init(0.0);
  }
  if (!isCombinedBZero || !isZeroB) {
    B = unbordered_grp->getX().createMultiVector(myWidth);
    B->init(0.0);
    
  }
  if (!isCombinedCZero || !isZeroC) {
    C = Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(myWidth,
								 myWidth));
    C->putScalar(0.0);
  }

  std::vector<int> idx1(underlyingWidth);
  for (int i=0; i<underlyingWidth; i++)
    idx1[i] = i;
  if (!isCombinedAZero) {
    Teuchos::RefCountPtr<NOX::Abstract::MultiVector> underlyingA = 
      A->subView(idx1);
    grp->fillA(*underlyingA);
  }
  if (!isCombinedBZero) {
    Teuchos::RefCountPtr<NOX::Abstract::MultiVector> underlyingB = 
      B->subView(idx1);
    grp->fillB(*underlyingB);
  }
  if (!isCombinedCZero) {
    NOX::Abstract::MultiVector::DenseMatrix underlyingC(Teuchos::View, 
							*C, 
							underlyingWidth, 
							underlyingWidth, 
							0, 0);
    grp->fillC(underlyingC);
  }

  std::vector<int> idx2(numConstraints);
  for (int i=0; i<numConstraints; i++)
    idx2[i] = underlyingWidth+i;
  if (!isZeroA) {
    Teuchos::RefCountPtr<NOX::Abstract::MultiVector> my_A_x = A->subView(idx2);
    NOX::Abstract::MultiVector::DenseMatrix my_A_p(Teuchos::View, *C,
						   underlyingWidth, 
						   numConstraints, 0, 
						   underlyingWidth);
    grp->extractSolutionComponent(*blockA, *my_A_x);
    grp->extractParameterComponent(false, *blockA, my_A_p);
  }

  if (!isZeroB) {
    Teuchos::RefCountPtr<NOX::Abstract::MultiVector> my_B_x = B->subView(idx2);
    NOX::Abstract::MultiVector::DenseMatrix my_B_p(Teuchos::View, *C,
						   numConstraints, 
						   underlyingWidth, 
						   underlyingWidth, 0);
    grp->extractSolutionComponent(*blockB_dx, *my_B_x);
    grp->extractParameterComponent(true, *blockB_dx, my_B_p);
  }

  if (!isZeroC) {
    NOX::Abstract::MultiVector::DenseMatrix my_CC(Teuchos::View, *C,
						  numConstraints, 
						  numConstraints, 
						  underlyingWidth, 
						  underlyingWidth);
    my_CC.assign(*blockC);
  }
    
  // set blocks in solver
  solver->setMatrixBlocksMultiVecConstraint(unbordered_grp, A, B, C);
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::Nested::initForSolve()
{
  return solver->initForSolve();
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::Nested::initForTransposeSolve()
{
  return solver->initForSolve();
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::Nested::apply(
			  const NOX::Abstract::MultiVector& X,
			  const NOX::Abstract::MultiVector::DenseMatrix& Y,
			  NOX::Abstract::MultiVector& U,
			  NOX::Abstract::MultiVector::DenseMatrix& V) const
{
  int num_cols = X.numVectors();
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> XX = 
    unbordered_grp->getX().createMultiVector(num_cols);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> UU = 
    unbordered_grp->getX().createMultiVector(num_cols);
  NOX::Abstract::MultiVector::DenseMatrix YY(myWidth, num_cols);
  NOX::Abstract::MultiVector::DenseMatrix VV(myWidth, num_cols);
  NOX::Abstract::MultiVector::DenseMatrix YY1(Teuchos::View, YY,
					      underlyingWidth, num_cols, 
					      0, 0);
  NOX::Abstract::MultiVector::DenseMatrix YY2(Teuchos::View, YY,
					      numConstraints, num_cols, 
					      underlyingWidth, 0);
  NOX::Abstract::MultiVector::DenseMatrix VV1(Teuchos::View, VV,
					      underlyingWidth, num_cols, 
					      0, 0);
  NOX::Abstract::MultiVector::DenseMatrix VV2(Teuchos::View, VV,
					      numConstraints, num_cols, 
					      underlyingWidth, 0);
  
  grp->extractSolutionComponent(X, *XX);
  grp->extractParameterComponent(false, X, YY1);
  YY2.assign(Y);

  NOX::Abstract::Group::ReturnType status = 
    solver->apply(*XX, YY, *UU, VV);

  V.assign(VV2);
  grp->loadNestedComponents(*UU, VV1, U);

  return status;
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::Nested::applyTranspose(
			  const NOX::Abstract::MultiVector& X,
			  const NOX::Abstract::MultiVector::DenseMatrix& Y,
			  NOX::Abstract::MultiVector& U,
			  NOX::Abstract::MultiVector::DenseMatrix& V) const
{
  int num_cols = X.numVectors();
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> XX = 
    unbordered_grp->getX().createMultiVector(num_cols);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> UU = 
    unbordered_grp->getX().createMultiVector(num_cols);
  NOX::Abstract::MultiVector::DenseMatrix YY(myWidth, num_cols);
  NOX::Abstract::MultiVector::DenseMatrix VV(myWidth, num_cols);
  NOX::Abstract::MultiVector::DenseMatrix YY1(Teuchos::View, YY,
					      underlyingWidth, num_cols, 
					      0, 0);
  NOX::Abstract::MultiVector::DenseMatrix YY2(Teuchos::View, YY,
					      numConstraints, num_cols, 
					      underlyingWidth, 0);
  NOX::Abstract::MultiVector::DenseMatrix VV1(Teuchos::View, VV,
					      underlyingWidth, num_cols, 
					      0, 0);
  NOX::Abstract::MultiVector::DenseMatrix VV2(Teuchos::View, VV,
					      numConstraints, num_cols, 
					      underlyingWidth, 0);
  
  grp->extractSolutionComponent(X, *XX);
  grp->extractParameterComponent(false, X, YY1);
  YY2.assign(Y);

  NOX::Abstract::Group::ReturnType status = 
    solver->applyTranspose(*XX, YY, *UU, VV);

  V.assign(VV2);
  grp->loadNestedComponents(*UU, VV1, U);

  return status;
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::Nested::applyInverse(
			      Teuchos::ParameterList& params,
			      const NOX::Abstract::MultiVector* F,
			      const NOX::Abstract::MultiVector::DenseMatrix* G,
			      NOX::Abstract::MultiVector& X,
			      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);

  if (isZeroF && isZeroG) {
    X.init(0.0);
    Y.putScalar(0.0);
  }

  int num_cols = X.numVectors();  
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> FF;
  if (!isZeroF) 
    FF = unbordered_grp->getX().createMultiVector(num_cols);
  NOX::Abstract::MultiVector::DenseMatrix GG(myWidth, num_cols);  
  GG.putScalar(0.0);
  
  if (!isZeroF) {
    NOX::Abstract::MultiVector::DenseMatrix GG1(Teuchos::View, GG,
						underlyingWidth, num_cols, 
						0, 0);
    grp->extractSolutionComponent(*F, *FF);
    grp->extractParameterComponent(false, *F, GG1);
  }
  if (!isZeroG) {
    NOX::Abstract::MultiVector::DenseMatrix GG2(Teuchos::View, GG,
						numConstraints, num_cols, 
						underlyingWidth, 0);
    GG2.assign(*G);
  }

  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> XX = 
    unbordered_grp->getX().createMultiVector(num_cols);
  NOX::Abstract::MultiVector::DenseMatrix YY(myWidth, num_cols);
  NOX::Abstract::MultiVector::DenseMatrix YY1(Teuchos::View, YY,
					      underlyingWidth, num_cols, 
					      0, 0);
  NOX::Abstract::MultiVector::DenseMatrix YY2(Teuchos::View, YY,
					      numConstraints, num_cols, 
					      underlyingWidth, 0);

  NOX::Abstract::Group::ReturnType status = 
    solver->applyInverse(params, FF.get(), &GG, *XX, YY);

  Y.assign(YY2);
  grp->loadNestedComponents(*XX, YY1, X);

  return status;
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::Nested::applyInverseTranspose(
			      Teuchos::ParameterList& params,
			      const NOX::Abstract::MultiVector* F,
			      const NOX::Abstract::MultiVector::DenseMatrix* G,
			      NOX::Abstract::MultiVector& X,
			      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);

  if (isZeroF && isZeroG) {
    X.init(0.0);
    Y.putScalar(0.0);
  }

  int num_cols = X.numVectors();  
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> FF;
  if (!isZeroF) 
    FF = unbordered_grp->getX().createMultiVector(num_cols);
  NOX::Abstract::MultiVector::DenseMatrix GG(myWidth, num_cols);  
  GG.putScalar(0.0);
  
  if (!isZeroF) {
    NOX::Abstract::MultiVector::DenseMatrix GG1(Teuchos::View, GG,
						underlyingWidth, num_cols, 
						0, 0);
    grp->extractSolutionComponent(*F, *FF);
    grp->extractParameterComponent(false, *F, GG1);
  }
  if (!isZeroG) {
    NOX::Abstract::MultiVector::DenseMatrix GG2(Teuchos::View, GG,
						numConstraints, num_cols, 
						underlyingWidth, 0);
    GG2.assign(*G);
  }

  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> XX = 
    unbordered_grp->getX().createMultiVector(num_cols);
  NOX::Abstract::MultiVector::DenseMatrix YY(myWidth, num_cols);
  NOX::Abstract::MultiVector::DenseMatrix YY1(Teuchos::View, YY,
					      underlyingWidth, num_cols, 
					      0, 0);
  NOX::Abstract::MultiVector::DenseMatrix YY2(Teuchos::View, YY,
					      numConstraints, num_cols, 
					      underlyingWidth, 0);

  NOX::Abstract::Group::ReturnType status = 
    solver->applyInverseTranspose(params, FF.get(), &GG, *XX, YY);

  Y.assign(YY2);
  grp->loadNestedComponents(*XX, YY1, X);

  return status;
}

