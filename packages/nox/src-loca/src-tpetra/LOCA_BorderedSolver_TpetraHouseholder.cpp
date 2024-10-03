// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_BorderedSolver_TpetraHouseholder.hpp"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H"

#include "LOCA_Thyra_Group.H"
#include "NOX_TpetraTypedefs.hpp"
#include "NOX_Thyra_MultiVector.H"
#include "Thyra_TpetraVectorSpace.hpp"
#include "Thyra_TpetraMultiVector.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "LOCA_Tpetra_LowRankUpdateRowMatrix.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "LOCA_BorderedSolver_LowerTriangularBlockElimination.H"
#include "LOCA_BorderedSolver_UpperTriangularBlockElimination.H"
#include "LOCA_Abstract_TransposeSolveGroup.H"
#include "LOCA_BorderedSolver_JacobianOperator.H"
#include "LOCA_BorderedSolver_ComplexOperator.H"
#include "LOCA_Hopf_ComplexMultiVector.H"

// To suppress unreachable return warnings on cuda
#include "Teuchos_CompilerCodeTweakMacros.hpp"

// For debugging output
#include <fstream>

// Forward declaration needed for ParameterList validation
namespace LOCA {
  namespace MultiContinuation {
    class ConstraintInterface;
  }
}

// Utility for extracting tpetra vector from nox vector
using ST = NOX::Scalar;
using LO = NOX::LocalOrdinal;
using GO = NOX::GlobalOrdinal;
using NT = NOX::NodeType;
namespace NOX {
  namespace Tpetra {

    TMultiVector& getTpetraMultiVector(NOX::Abstract::MultiVector& v) {
      auto& v_thyra = *(dynamic_cast<NOX::Thyra::MultiVector&>(v).getThyraMultiVector());
      auto& v_tpetra = *(dynamic_cast<::Thyra::TpetraMultiVector<ST,LO,GO,NT>&>(v_thyra).getTpetraMultiVector());
      return v_tpetra;
    }

    const TMultiVector& getTpetraMultiVector(const NOX::Abstract::MultiVector& v) {
      const auto& v_thyra = *(dynamic_cast<const NOX::Thyra::MultiVector&>(v).getThyraMultiVector());
      const auto& v_tpetra = *(dynamic_cast<const ::Thyra::TpetraMultiVector<ST,LO,GO,NT>&>(v_thyra).getConstTpetraMultiVector());
      return v_tpetra;
    }

    Teuchos::RCP<TMultiVector> getTpetraMultiVector(const Teuchos::RCP<NOX::Abstract::MultiVector>& v) {
      auto v_thyra = Teuchos::rcp_dynamic_cast<NOX::Thyra::MultiVector>(v)->getThyraMultiVector();
      auto v_tpetra = Teuchos::rcp_dynamic_cast<::Thyra::TpetraMultiVector<ST,LO,GO,NT>>(v_thyra)->getTpetraMultiVector();
      return v_tpetra;
    }

    Teuchos::RCP<const TMultiVector> getTpetraMultiVector(const Teuchos::RCP<const NOX::Abstract::MultiVector>& v) {
      auto v_thyra = Teuchos::rcp_dynamic_cast<const NOX::Thyra::MultiVector>(v)->getThyraMultiVector();
      auto v_tpetra = Teuchos::rcp_dynamic_cast<const ::Thyra::TpetraMultiVector<ST,LO,GO,NT>>(v_thyra)->getConstTpetraMultiVector();
      return v_tpetra;
    }

  }
}

LOCA::BorderedSolver::TpetraHouseholder::
TpetraHouseholder(const Teuchos::RCP<LOCA::GlobalData>& global_data,
                  const Teuchos::RCP<LOCA::Parameter::SublistParser>& /* topParams */,
                  const Teuchos::RCP<Teuchos::ParameterList>& slvrParams):
  globalData(global_data),
  solverParams(slvrParams),
  grp(),
  op(),
  A(),
  B(),
  C(),
  constraints(),
  qrFact(globalData),
  house_x(),
  house_p(),
  T(),
  R(),
  U(),
  V(),
  house_x_trans(),
  house_p_trans(),
  T_trans(),
  R_trans(),
  U_trans(),
  V_trans(),
  Ablock(),
  Bblock(),
  Ascaled(),
  Bscaled(),
  Cscaled(),
  // linSys(),
  tpetraOp(),
  // baseMap(),
  // globalMap(),
  numConstraints(0),
  isZeroA(true),
  isZeroB(true),
  isZeroC(true),
  isValidForSolve(false),
  isValidForTransposeSolve(false),
  dblas(),
  scale_rows(true),
  scale_vals(),
  precMethod(JACOBIAN),
  includeUV(false),
  use_P_For_Prec(false),
  isComplex(false),
  omega(0.0)
{
  Teuchos::ParameterList validParams;
  validParams.set("Bordered Solver Method", "Householder");
  validParams.set("Constraint Object",Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>(Teuchos::null));
  validParams.set("Constraint Parameter Names",Teuchos::RCP<std::vector<std::string>>(Teuchos::null));
  validParams.set("Scale Augmented Rows", true);
  validParams.set("Preconditioner Method", "Jacobian", "Matrix to use for Preconditioning",
    rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("Jacobian","SWM"))));
  validParams.set("Include UV In Preconditioner", false);
  validParams.set("Use P For Preconditioner", false);
  solverParams->validateParametersAndSetDefaults(validParams);

  scale_rows = solverParams->get<bool>("Scale Augmented Rows");
  std::string prec_method =
    solverParams->get<std::string>("Preconditioner Method");
  if (prec_method == "Jacobian")
    precMethod = JACOBIAN;
  else if (prec_method == "SMW")
    precMethod = SMW;
  else
    globalData->locaErrorCheck->throwError(
        "LOCA::BorderedSolver::TpetraHouseholder::TpetraHouseholder()",
        "Unknown preconditioner method!  Choices are Jacobian, SMW");

  includeUV =
    solverParams->get<bool>("Include UV In Preconditioner");
  use_P_For_Prec =
    solverParams->get<bool>("Use P For Preconditioner");
}

LOCA::BorderedSolver::TpetraHouseholder::~TpetraHouseholder()
{
}

void
LOCA::BorderedSolver::TpetraHouseholder::
setMatrixBlocks(const Teuchos::RCP<const LOCA::BorderedSolver::AbstractOperator>& op_,
                const Teuchos::RCP<const NOX::Abstract::MultiVector>& blockA,
                const Teuchos::RCP<const LOCA::MultiContinuation::ConstraintInterface>& blockB,
                const Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>& blockC)
{
  std::string callingFunction =
    "LOCA::BorderedSolver::TpetraHouseholder::setMatrixBlocks";

  op = op_;
  A = blockA;

  // Cast constraints to a ConstraintInterfaceMVDX
  constraints = Teuchos::rcp_dynamic_cast<const LOCA::MultiContinuation::ConstraintInterfaceMVDX>(blockB);
  if (constraints.get() == NULL)
    globalData->locaErrorCheck->throwError(
         callingFunction,
         "Constraints object must be of type ConstraintInterfaceMVDX");
  C = blockC;

  // Determine which blocks are zero
  isZeroA = (A.get() == NULL);
  isZeroB = constraints->isDXZero();
  isZeroC = (C.get() == NULL);

  // Get multivec constraint if it is nonzero
  if (isZeroB)
    B = Teuchos::null;
  else
    B = Teuchos::rcp(constraints->getDX(), false);

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

  if (isZeroB)
    numConstraints = C->numRows();
  else
    numConstraints = B->numVectors();

  // We only use the Householder technique if A and B are nonzero
  // (otherwise we use a block elimination).  If C is zero, we just create
  // a matrix of zeros and apply the standard algorithm
  if (isZeroC && !isZeroA && !isZeroB) {

    Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> tmpC =
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(
                               B->numVectors(),
                               B->numVectors()));
    tmpC->putScalar(0.0);
    C = tmpC;
    isZeroC = false;
  }

  // Set flags indicating we have to compute constraint factorizations
  isValidForSolve = false;
  isValidForTransposeSolve = false;

  Teuchos::RCP<const LOCA::BorderedSolver::JacobianOperator> jacOp =
    Teuchos::rcp_dynamic_cast<const LOCA::BorderedSolver::JacobianOperator>(op);
  Teuchos::RCP<const LOCA::BorderedSolver::ComplexOperator> complexOp =
    Teuchos::rcp_dynamic_cast<const LOCA::BorderedSolver::ComplexOperator>(op);

  if (jacOp != Teuchos::null) {

    isComplex = false;

    // Group must be a Thyra (Tpetra) group
    Teuchos::RCP<const LOCA::Thyra::Group> constGrp =
      Teuchos::rcp_dynamic_cast<const LOCA::Thyra::Group>(jacOp->getGroup());
    if (constGrp.get() == NULL)
      globalData->locaErrorCheck->throwError(
                    callingFunction,
                    "Group object must be an Epetra group");

    // Get A block
    Ablock = A;

    // Get B block
    Bblock = B;

    // Cast to non-const group
    grp = Teuchos::rcp_const_cast<LOCA::Thyra::Group>(constGrp);

    // Get linear system, and jacobian
    //ROGER linSys = grp->getLinearSystem();
    using ttlop = ::Thyra::TpetraLinearOp<NOX::Scalar,NOX::LocalOrdinal,NOX::GlobalOrdinal,NOX::NodeType>;
    auto constTpetraOperator = Teuchos::rcp_dynamic_cast<const ttlop>(grp->getJacobianOperator(),true);
    tpetraOp = Teuchos::rcp_const_cast<ttlop>(constTpetraOperator)->getTpetraOperator();
  }
  else if (complexOp != Teuchos::null) {

    // Complex not yet supported
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                               "ERROR: Thyra/Tpetra Group doesn't support COMPLEX bordered solve!");

    /*
    isComplex = true;

    omega = complexOp->getFrequency();

    // Group must be an Epetra group
    Teuchos::RCP<const LOCA::Epetra::Group> constGrp =
      Teuchos::rcp_dynamic_cast<const LOCA::Epetra::Group>(complexOp->getGroup());
    if (constGrp.get() == NULL)
      globalData->locaErrorCheck->throwError(
                    callingFunction,
                    "Group object must be an Epetra group");

    // Cast to non-const group
    grp = Teuchos::rcp_const_cast<LOCA::Epetra::Group>(constGrp);

    // Get linear system and complex
    linSys = grp->getComplexLinearSystem();
    epetraOp = linSys->getJacobianOperator();

    // Get maps for complex vectors
    grp->getComplexMaps(baseMap, globalMap);

    // Get A block
    if (!isZeroA)
      Ablock = createBlockMV(*A);

    // Get A block
    if (!isZeroB)
      Bblock = createBlockMV(*B);
    */
  }
  else {
    globalData->locaErrorCheck->throwError(
              callingFunction,
              std::string("Op argument must be of type !\n") +
                  std::string("LOCA::BorderedSolver::JacobianOperator or \n") +
              std::string("LOCA::BorderedSolver::ComplexOperator."));
  }

  Ascaled = Teuchos::null;
  Bscaled = Teuchos::null;
  Cscaled = Teuchos::null;

  if (Ablock != Teuchos::null)
    Ascaled = Ablock->clone();
  if (Bblock != Teuchos::null)
    Bscaled = Bblock->clone();
  if (C != Teuchos::null) {
    Cscaled =
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(
                                C->numRows(),
                                C->numCols()));
    Cscaled->assign(*C);
  }

}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::TpetraHouseholder::initForSolve()
{
  NOX::Abstract::Group::ReturnType res = NOX::Abstract::Group::Ok;

  // Check if we need to compute anything
  if (isValidForSolve)
    return res;

  // Only compute QR factorization if A and B are nonzero
  if (!isZeroA && !isZeroB) {

    // Scale rows
    if (scale_rows) {
      scale_vals.resize(numConstraints);
      double sn = std::sqrt( static_cast<double>(Bblock->length()) );
      for (int i=0; i<numConstraints; i++) {
        scale_vals[i] = (*Bblock)[i].norm() / sn;
        double t = 0.0;
        for (int j=0; j<numConstraints; j++)
          t += (*C)(i,j) * (*C)(i,j);
        scale_vals[i] += std::sqrt(t);
        scale_vals[i] = 1.0 / scale_vals[i];
        (*Bscaled)[i].scale(scale_vals[i]);
        for (int j=0; j<numConstraints; j++)
          (*Cscaled)(i,j) *= scale_vals[i];
      }
    }

    // Allocate vectors and matrices for factorization
    if (house_x.get() == NULL || house_x->numVectors() != numConstraints) {
      house_x = Bblock->clone(NOX::ShapeCopy);
      U = house_x->clone(NOX::ShapeCopy);
      V = house_x->clone(NOX::ShapeCopy);
      house_p.reshape(numConstraints, numConstraints);
      T.reshape(numConstraints, numConstraints);
      R.reshape(numConstraints, numConstraints);
    }

    // Factor constraints
    qrFact.computeQR(*Cscaled, *Bscaled, true, house_p, *house_x, T, R);

    // Compute U & V in P operator
    res = computeUV(house_p, *house_x, T, *Ablock, *U, *V, false);
    globalData->locaErrorCheck->checkReturnType(res,
      "LOCA::BorderedSolver::TpetraHouseholder::initForSolve()");
  }

  isValidForSolve = true;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::TpetraHouseholder::initForTransposeSolve()
{
  NOX::Abstract::Group::ReturnType res = NOX::Abstract::Group::Ok;

  // Check if we need to compute anything
  if (isValidForTransposeSolve)
    return res;

  // Only compute QR factorization if A and B are nonzero
  if (!isZeroA && !isZeroB) {

    // Scale rows
    if (scale_rows) {
      scale_vals.resize(numConstraints);
      double sn = std::sqrt( static_cast<double>(Ablock->length()) );
      for (int i=0; i<numConstraints; i++) {
    scale_vals[i] = (*Ablock)[i].norm() / sn;
    double t = 0.0;
    for (int j=0; j<numConstraints; j++)
      t += (*C)(j,i) * (*C)(j,i);
    scale_vals[i] += std::sqrt(t);
    scale_vals[i] = 1.0 / scale_vals[i];
    (*Ascaled)[i].scale(scale_vals[i]);
    for (int j=0; j<numConstraints; j++)
      (*Cscaled)(j,i) *= scale_vals[i];
      }
    }

    // Allocate vectors and matrices for factorization
    if (house_x_trans.get() == NULL ||
    house_x_trans->numVectors() != numConstraints) {
      house_x_trans = Ablock->clone(NOX::ShapeCopy);
      U_trans = house_x_trans->clone(NOX::ShapeCopy);
      V_trans = house_x_trans->clone(NOX::ShapeCopy);
      house_p_trans.reshape(numConstraints, numConstraints);
      T_trans.reshape(numConstraints, numConstraints);
      R_trans.reshape(numConstraints, numConstraints);
    }

    // Factor constraints for transposed system
    qrFact.computeQR(*Cscaled, *Ascaled, false, house_p_trans, *house_x_trans,
             T_trans, R_trans);

    // Compute U & V in transposed P operator
    res = computeUV(house_p_trans, *house_x_trans, T_trans, *Bblock, *U_trans, *V_trans, true);
    globalData->locaErrorCheck->checkReturnType(res,
      "LOCA::BorderedSolver::TpetraHouseholder::initForTransposeSolve()");
  }

  isValidForTransposeSolve = true;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::TpetraHouseholder::apply(
              const NOX::Abstract::MultiVector& X,
              const NOX::Abstract::MultiVector::DenseMatrix& Y,
              NOX::Abstract::MultiVector& inU,
              NOX::Abstract::MultiVector::DenseMatrix& inV) const
{
  // Compute J*X
  NOX::Abstract::Group::ReturnType status = op->apply(X, inU);

  // Compute J*X + A*Y
  if (!isZeroA)
    inU.update(Teuchos::NO_TRANS, 1.0, *A, Y, 1.0);

  // Compute B^T*X
  if (!isZeroB)
    constraints->multiplyDX(1.0, X, inV);

  // Compute B^T*X + C*Y
  if (!isZeroC) {
    int e;
    if (isZeroB)
      e = inV.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, *C, Y, 0.0);
    else
      e = inV.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, *C, Y, 1.0);
    if (e < 0)
      status = NOX::Abstract::Group::Failed;
  }

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::TpetraHouseholder::applyTranspose(
              const NOX::Abstract::MultiVector& X,
              const NOX::Abstract::MultiVector::DenseMatrix& Y,
              NOX::Abstract::MultiVector& inU,
              NOX::Abstract::MultiVector::DenseMatrix& inV) const
{
  // Compute J*X
  NOX::Abstract::Group::ReturnType status = op->applyTranspose(X, inU);

  // Compute J*X + B*Y
  if (!isZeroA)
    constraints->addDX(Teuchos::NO_TRANS, 1.0, Y, 1.0, inU);

  // Compute A^T*X
  if (!isZeroB)
    X.multiply(1.0, *A, inV);

  // Compute A^T*X + C^T*Y
  if (!isZeroC) {
    int e;
    if (isZeroB)
      e = inV.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, *C, Y, 0.0);
    else
      e = inV.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, *C, Y, 1.0);
    if (e < 0)
      status = NOX::Abstract::Group::Failed;
  }

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::TpetraHouseholder::applyInverse(
                  Teuchos::ParameterList& params,
                  const NOX::Abstract::MultiVector* F,
                  const NOX::Abstract::MultiVector::DenseMatrix* G,
                  NOX::Abstract::MultiVector& X,
                  NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);

  // If F and G are zero, return zero
  if (isZeroF && isZeroG) {
    X.init(0.0);
    Y.putScalar(0.0);
    return NOX::Abstract::Group::Ok;
  }

  if (isZeroA) {
    LOCA::BorderedSolver::LowerTriangularBlockElimination ltbe(globalData);
    return ltbe.solve(params, *op, *constraints, *C, F, G, X, Y);
  }

  if (isZeroB) {
    LOCA::BorderedSolver::UpperTriangularBlockElimination utbe(globalData);
    return utbe.solve(params, *op, A.get(), *C, F, G, X, Y);
  }

  // Scale G
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> Gscaled;
  if (G != NULL) {
    Gscaled =
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(G->numRows(),
                                   G->numCols()));
    if (scale_rows)
      for (int i=0; i<G->numRows(); i++)
    for (int j=0; j<G->numCols(); j++)
      (*Gscaled)(i,j) = scale_vals[i]*(*G)(i,j);
    else
      Gscaled->assign(*G);
  }

  NOX::Abstract::Group::ReturnType res;
  if (!isComplex)
    res = solve(params, F, Gscaled.get(), X, Y);
  else {
    Teuchos::RCP<NOX::Abstract::MultiVector> blockF;
    if (!isZeroF)
      blockF = createBlockMV(*F);
    Teuchos::RCP<NOX::Abstract::MultiVector> blockX = createBlockMV(X);
    res = solve(params, blockF.get(), Gscaled.get(), *blockX, Y);
    setBlockMV(*blockX, X);
  }

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::TpetraHouseholder::applyInverseTranspose(
                  Teuchos::ParameterList& params,
                  const NOX::Abstract::MultiVector* F,
                  const NOX::Abstract::MultiVector::DenseMatrix* G,
                  NOX::Abstract::MultiVector& X,
                  NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);

  // If F and G are zero, return zero
  if (isZeroF && isZeroG) {
    X.init(0.0);
    Y.putScalar(0.0);
    return NOX::Abstract::Group::Ok;
  }

  // If A or B is zero, we use bordering, which requires a transpose solve
  Teuchos::RCP<const LOCA::Abstract::TransposeSolveGroup> ts_grp;
  if (isZeroA || isZeroB) {
    ts_grp =
      Teuchos::rcp_dynamic_cast<const LOCA::Abstract::TransposeSolveGroup>(grp);
    if (ts_grp == Teuchos::null)
      globalData->locaErrorCheck->throwError(
    "LOCA::BorderedSolver::TpetraHouseholder::applyInverseTranspose()",
    std::string("Group must implement the LOCA::Abstract::TransposeSolveGroup")
    + std::string(" interface in order to solve the transpose of the bordered")
    + std::string(" system via bordering."));
  }

  if (isZeroA) {
    LOCA::BorderedSolver::UpperTriangularBlockElimination utbe(globalData);
    return utbe.solveTranspose(params, *op, B.get(), *C, F, G, X, Y);
  }

  else if (isZeroB) {
    LOCA::BorderedSolver::LowerTriangularBlockElimination ltbe(globalData);
    return ltbe.solveTranspose(params, *op, *A, *C, F, G, X, Y);
  }

  // Scale G
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> Gscaled;
  if (G != NULL) {
    Gscaled =
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(G->numRows(),
                                   G->numCols()));
    if (scale_rows)
      for (int i=0; i<G->numRows(); i++)
    for (int j=0; j<G->numCols(); j++)
      (*Gscaled)(i,j) = scale_vals[i]*(*G)(i,j);
    else
      Gscaled->assign(*G);
  }

  NOX::Abstract::Group::ReturnType res;
  if (!isComplex)
    res = solveTranspose(params, F, Gscaled.get(), X, Y);
  else {
    Teuchos::RCP<NOX::Abstract::MultiVector> blockF;
    if (!isZeroF)
      blockF = createBlockMV(*F);
    Teuchos::RCP<NOX::Abstract::MultiVector> blockX = createBlockMV(X);
    res = solveTranspose(params, blockF.get(), Gscaled.get(), *blockX, Y);
    setBlockMV(*blockX, X);
  }

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::TpetraHouseholder::solve(
                  Teuchos::ParameterList& params,
                  const NOX::Abstract::MultiVector* F,
                  const NOX::Abstract::MultiVector::DenseMatrix* G,
                  NOX::Abstract::MultiVector& X,
                  NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  std::string callingFunction =
    "LOCA::BorderedSolver::TpetraHouseholder::applyInverse()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);

  // Make sure constraint factorization has been computed
  if (!isValidForSolve)
    globalData->locaErrorCheck->throwError(
           callingFunction,
           std::string("applyInverse() called with invalid constraint") +
           std::string(" factorizations.  Call initForSolve() first."));

  Teuchos::RCP<const NOX::Abstract::MultiVector> cRHS;

  if (!isZeroG) {

    Teuchos::RCP<NOX::Abstract::MultiVector> tmp_x = X.clone(NOX::ShapeCopy);
    Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> tmp_y =
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(G->numRows(),
                                                               G->numCols()));
    Teuchos::RCP<NOX::Abstract::MultiVector> RHS = X.clone(NOX::ShapeCopy);

    // Compute Z_y = R^-T * G
    Y.assign(*G);
    dblas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS,
               Teuchos::NON_UNIT_DIAG,
               G->numRows(), G->numCols(), 1.0, R.values(),
               R.numRows(), Y.values(), Y.numRows());

    // Compute P*[Z_y; 0]
    qrFact.applyCompactWY(house_p, *house_x, T, &Y, NULL, *tmp_y, *tmp_x, false);

    // Compute -[A J]*P*[Z_y 0]
    tpetraOp->apply(NOX::Tpetra::getTpetraMultiVector(*tmp_x),
                    NOX::Tpetra::getTpetraMultiVector(*RHS));

    // Left over from check to make sure above apply was successful. Can probably remove.
    status = NOX::Abstract::Group::Ok;

    RHS->update(Teuchos::NO_TRANS, -1.0, *Ablock, *tmp_y, -1.0);

    // Compute F - [A J]*P*[Z_y 0]
    if (!isZeroF)
      RHS->update(1.0, *F, 1.0);

    cRHS = RHS;
  }
  else
    cRHS = Teuchos::rcp(F, false);

  // Create operator for P = J + U*V^T
  auto tpetra_U = NOX::Tpetra::getTpetraMultiVector(U);
  auto tpetra_V = NOX::Tpetra::getTpetraMultiVector(V);

  // Create a row-matrix version of P if J is a row matrix
  Teuchos::RCP<NOX::TRowMatrix> jac_rowmatrix =
    Teuchos::rcp_dynamic_cast<NOX::TRowMatrix>(tpetraOp);
  Teuchos::RCP<NOX::TOperator> tmpOp;
  if (jac_rowmatrix != Teuchos::null) {
    tmpOp = Teuchos::rcp(new LOCA::Tpetra::LowRankUpdateRowMatrix(globalData,
                                                                  jac_rowmatrix,
                                                                  tpetra_U,
                                                                  tpetra_V,
                                                                  false,
                                                                  includeUV));
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                               "ERROR - LOCA::BorderedSolver::TpetraHouseholder::solve(): LowRankUpdateOp NOT IMPLEMENTED YET!");
    // tmpOp = Teuchos::rcp(new LOCA::Epetra::LowRankUpdateOp(globalData,
    //                                                        epetraOp,
    //                                                        epetra_U,
    //                                                        epetra_V,
    //                                                        false));
  }

  // Allocate a separate matrix for the preconditioner. Don't want to
  // corrupt J with U*V^T if not using P for Prec (we can't for
  // tpetra).  Copy J values and add in U*V^T if it's a CRS matrix
  // and we aren't using P for the preconditioner
  if (includeUV && !use_P_For_Prec) {
    auto jac_crs = Teuchos::rcp_dynamic_cast<NOX::TCrsMatrix>(tpetraOp,true);
    if (tpetraPrecMatrix.is_null()) {
      tpetraPrecMatrix = Teuchos::rcp(new NOX::TCrsMatrix(*jac_crs,Teuchos::Copy));
      Teuchos::RCP<::Thyra::VectorSpaceBase<double>> domain = ::Thyra::tpetraVectorSpace<double>(tpetraPrecMatrix->getDomainMap());
      Teuchos::RCP<::Thyra::VectorSpaceBase<double>> range = ::Thyra::tpetraVectorSpace<double>(tpetraPrecMatrix->getRangeMap());
      auto prec_thyra_op = Teuchos::rcp(new ::Thyra::TpetraLinearOp<ST,LO,GO,NT>);
      prec_thyra_op->initialize(range,domain,tpetraPrecMatrix);
      Teuchos::RCP<::Thyra::LinearOpBase<double>> tmp_for_thyra_ambiguity = prec_thyra_op;
      prec_losb = Teuchos::rcp(new ::Thyra::DefaultLinearOpSource<double>(tmp_for_thyra_ambiguity));
    }

    // Copy J values into preconditioner matrix
    //*tpetraPrecMatrix = *jac_crs;
    {
      tpetraPrecMatrix->resumeFill();
      auto jac_view = jac_crs->getLocalMatrixDevice().values;
      auto prec_view = tpetraPrecMatrix->getLocalMatrixDevice().values;
      Kokkos::deep_copy(prec_view,jac_view);
      tpetraPrecMatrix->fillComplete();
    }

    bool print_debug = false;
    if (print_debug) {
      std::fstream fsj("jac_matrix_before.out",std::fstream::out|std::fstream::trunc);
      Teuchos::FancyOStream tfsj(Teuchos::rcpFromRef(fsj));
      jac_crs->describe(tfsj,Teuchos::VERB_EXTREME);
      std::fstream fsp("prec_matrix_before.out",std::fstream::out|std::fstream::trunc);
      Teuchos::FancyOStream tfsp(Teuchos::rcpFromRef(fsp));
      tpetraPrecMatrix->describe(tfsp,Teuchos::VERB_EXTREME);
      std::fstream fsu("u_vector.out",std::fstream::out|std::fstream::trunc);
      Teuchos::FancyOStream tfsu(Teuchos::rcpFromRef(fsu));
      tpetra_U->describe(tfsu,Teuchos::VERB_EXTREME);
      std::fstream fsv("v_vector.out",std::fstream::out|std::fstream::trunc);
      Teuchos::FancyOStream tfsv(Teuchos::rcpFromRef(fsv));
      tpetra_V->describe(tfsv,Teuchos::VERB_EXTREME);
    }

    // Update locally owned non-zero values for U*V^T
    updateCrsMatrixForPreconditioner(*U, *V, *tpetraPrecMatrix);

    if (print_debug) {
      std::fstream fsj("jac_matrix_after.out",std::fstream::out|std::fstream::trunc);
      Teuchos::FancyOStream tfsj(Teuchos::rcpFromRef(fsj));
      jac_crs->describe(tfsj,Teuchos::VERB_EXTREME);
      std::fstream fsp("prec_matrix_after.out",std::fstream::out|std::fstream::trunc);
      Teuchos::FancyOStream tfsp(Teuchos::rcpFromRef(fsp));
      tpetraPrecMatrix->describe(tfsp,Teuchos::VERB_EXTREME);
    }

    grp->setPreconditionerMatrix(prec_losb);
  }

  // Set operator in solver to compute preconditioner
  // if (use_P_For_Prec)
  //   linSys->setJacobianOperatorForSolve(tmpOp);
  // else
  //   linSys->setJacobianOperatorForSolve(tpetraOp);

  // Now compute preconditioner
  // linSys->destroyPreconditioner();
  // linSys->createPreconditioner(solution_vec, params, false);

  // Now recompute J if we modified it
  // if (includeUV && !use_P_For_Prec && jac_crs != Teuchos::null) {
  //   grp->setX(solution_vec);
  //   if (isComplex)
  //     grp->computeComplex(omega);
  //   else
  //     grp->computeJacobian();
  // }

  // Set operator for P as Jacobian in solver
  auto save_original_jacobian = grp->getNonconstJacobianOperator();
  auto P_thyra = Teuchos::rcp(new ::Thyra::TpetraLinearOp<ST,LO,GO,NT>);
  auto range_space = ::Thyra::tpetraVectorSpace<ST,LO,GO,NT>(tmpOp->getRangeMap());
  auto domain_space = ::Thyra::tpetraVectorSpace<ST,LO,GO,NT>(tmpOp->getDomainMap());
  P_thyra->initialize(range_space,domain_space,tmpOp);
  grp->setJacobianOperator(P_thyra);

  /*
  // Set preconditioner
  Teuchos::RCP<Epetra_Operator> prec_op;
  Teuchos::RCP<Epetra_Operator> epetraPrecOp;
  if (precMethod == SMW) {
    epetraPrecOp = linSys->getGeneratedPrecOperator();
    prec_op = Teuchos::rcp(new LOCA::Epetra::LowRankUpdateOp(globalData,
                                                             epetraPrecOp,
                                                             epetra_U,
                                                             epetra_V,
                                                             true));
    linSys->setPrecOperatorForSolve(prec_op);
  }
  */

  // Solve for each RHS
  int m = X.numVectors();
  X.init(0.0);
  for (int i=0; i<m; i++) {
    status =
      grp->applyJacobianInverse(
             params,
             dynamic_cast<const NOX::Thyra::Vector&>((*cRHS)[i]),
             dynamic_cast<NOX::Thyra::Vector&>(X[i]));

    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                                             finalStatus,
                                                             callingFunction);
  }

  qrFact.applyCompactWY(house_p, *house_x, T, Y, X, isZeroG, false, false);

  // Set original operators in linear system
  // if (precMethod == SMW)
  //   linSys->setPrecOperatorForSolve(epetraPrecOp);
  // linSys->destroyPreconditioner();

  grp->setJacobianOperator(save_original_jacobian);

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::TpetraHouseholder::solveTranspose(
                  Teuchos::ParameterList& params,
                  const NOX::Abstract::MultiVector* F,
                  const NOX::Abstract::MultiVector::DenseMatrix* G,
                  NOX::Abstract::MultiVector& X,
                  NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                             "ERROR - LOCA::BorderedSolver::TpetraHouseholder::solveTranspose(): TRANSPOSE SOLVE NOT SUPPORTED YET!");

  /*
  std::string callingFunction =
    "LOCA::BorderedSolver::TpetraHouseholder::applyInverseTranspose()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);

  // Make sure constraint factorization has been computed
  if (!isValidForTransposeSolve)
    globalData->locaErrorCheck->throwError(
        callingFunction,
        std::string("applyInverseTranspose() called with invalid constraint") +
        std::string(" factorizations.  Call initForSolve() first."));

  Teuchos::RCP<const NOX::Abstract::MultiVector> cRHS;

  if (!isZeroG) {

    Teuchos::RCP<NOX::Epetra::MultiVector> tmp_x =
      Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(X.clone(NOX::ShapeCopy));
    Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> tmp_y =
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(G->numRows(),
                                   G->numCols()));
    Teuchos::RCP<NOX::Epetra::MultiVector> RHS =
      Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(X.clone(NOX::ShapeCopy));

    // Compute Z_y = R^-T * G
    Y.assign(*G);
    dblas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS,
           Teuchos::NON_UNIT_DIAG,
           G->numRows(), G->numCols(), 1.0, R_trans.values(),
           R_trans.numRows(), Y.values(), Y.numRows());

    // Compute P*[Z_y; 0]
    qrFact.applyCompactWY(house_p_trans, *house_x_trans, T_trans, &Y, NULL,
              *tmp_y, *tmp_x, false);

    // Compute -[B J^T]*P*[Z_y 0]
    epetraOp->SetUseTranspose(true);
    int res = epetraOp->Apply(tmp_x->getEpetraMultiVector(),
                  RHS->getEpetraMultiVector());
    epetraOp->SetUseTranspose(false);
    if (res == 0)
      status = NOX::Abstract::Group::Ok;
    else
      status = NOX::Abstract::Group::Failed;
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
    RHS->update(Teuchos::NO_TRANS, -1.0, *Bblock, *tmp_y, -1.0);

    // Compute F - [B J^T]*P*[Z_y 0]
    if (!isZeroF)
      RHS->update(1.0, *F, 1.0);

    cRHS = RHS;
  }
  else
    cRHS = Teuchos::rcp(F, false);

  // Get solution vector
  const NOX::Epetra::Vector& solution_vec =
    dynamic_cast<const NOX::Epetra::Vector&>(grp->getX());

  // Instantiate transpose solver
  LOCA::Epetra::TransposeLinearSystem::Factory tls_factory(globalData);
//   Teuchos::RCP<LOCA::Epetra::TransposeLinearSystem::AbstractStrategy> tls_strategy = tls_factory.create(Teuchos::rcp(&params, false), linSys);
  Teuchos::RCP<LOCA::Epetra::TransposeLinearSystem::AbstractStrategy> tls_strategy = tls_factory.create(solverParams, linSys);

  // Compute Jacobian transpose (J^T)
  tls_strategy->createJacobianTranspose();
  Teuchos::RCP<Epetra_Operator> jac_trans =
    tls_strategy->getJacobianTransposeOperator();

  // Create operator for P = J^T + U*V^T
  Teuchos::RCP<NOX::Epetra::MultiVector> nox_epetra_U =
    Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(U_trans);
  Teuchos::RCP<Epetra_MultiVector> epetra_U =
    Teuchos::rcp(&(nox_epetra_U->getEpetraMultiVector()), false);
  Teuchos::RCP<NOX::Epetra::MultiVector> nox_epetra_V =
    Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(V_trans);
  Teuchos::RCP<Epetra_MultiVector> epetra_V =
    Teuchos::rcp(&(nox_epetra_V->getEpetraMultiVector()), false);

  // Create a row-matrix version of P if J^T is a row matrix
  Teuchos::RCP<Epetra_RowMatrix> jac_trans_rowmatrix =
    Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(jac_trans);
  Teuchos::RCP<Epetra_Operator> tmpOp;
  if (jac_trans_rowmatrix != Teuchos::null)
    tmpOp = Teuchos::rcp(new LOCA::Epetra::LowRankUpdateRowMatrix(
                              globalData,
                              jac_trans_rowmatrix,
                              epetra_U,
                              epetra_V,
                              false,
                              includeUV));
  else
    tmpOp = Teuchos::rcp(new LOCA::Epetra::LowRankUpdateOp(globalData,
                            jac_trans,
                            epetra_U,
                            epetra_V,
                            false));

  // Overwrite J^T with J^T + U*V^T if it's a CRS matrix and we aren't
  // using P for the preconditioner
  Teuchos::RCP<Epetra_CrsMatrix> jac_trans_crs;
  if (includeUV && !use_P_For_Prec) {
    jac_trans_crs = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(jac_trans);
    if (jac_trans_crs != Teuchos::null) {
      updateJacobianForPreconditioner(*U_trans, *V_trans, *jac_trans_crs);
    }
  }

   // Set operator in solver to compute preconditioner
  if (use_P_For_Prec)
    tls_strategy->setJacobianTransposeOperator(tmpOp);
  else
    tls_strategy->setJacobianTransposeOperator(jac_trans);

  // Now compute preconditioner
  tls_strategy->createTransposePreconditioner(solution_vec, params);

  // Now recompute J^T if we modified it
  if (includeUV && !use_P_For_Prec && jac_trans_crs != Teuchos::null) {
    grp->setX(solution_vec);
    if (isComplex)
      grp->computeComplex(omega);
    else
      grp->computeJacobian();
    tls_strategy->createJacobianTranspose();
  }

  // Set operator for P in transpose solver
  tls_strategy->setJacobianTransposeOperator(tmpOp);

  // Set preconditioner
  Teuchos::RCP<Epetra_Operator> prec_op;
  Teuchos::RCP<Epetra_Operator> epetraPrecOp;
  if (precMethod == SMW) {
    epetraPrecOp = tls_strategy->getTransposePreconditioner();
    prec_op = Teuchos::rcp(new LOCA::Epetra::LowRankUpdateOp(globalData,
                                 epetraPrecOp,
                                 epetra_U,
                                 epetra_V,
                                 true));
    tls_strategy->setTransposePreconditioner(prec_op);
  }

  // Solve for each RHS
  int m = X.numVectors();
  X.init(0.0);
  for (int i=0; i<m; i++) {
    bool stat =
      tls_strategy->applyJacobianTransposeInverse(
              params,
              dynamic_cast<const NOX::Epetra::Vector&>((*cRHS)[i]),
              dynamic_cast<NOX::Epetra::Vector&>(X[i]));
    if (stat == true)
      status = NOX::Abstract::Group::Ok;
    else
      status = NOX::Abstract::Group::NotConverged;
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);


  }

  qrFact.applyCompactWY(house_p_trans, *house_x_trans, T_trans, Y, X, isZeroG,
            false, false);

  epetraOp->SetUseTranspose(false);

  // Set original operators in linear system
  linSys->setJacobianOperatorForSolve(epetraOp);
  linSys->destroyPreconditioner();

  return finalStatus;
  */

  // Temporary dummy return to quiet warnings. Remove when the code
  // above gets fixed.
  TEUCHOS_UNREACHABLE_RETURN(NOX::Abstract::Group::Failed);
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::TpetraHouseholder::computeUV(
                const NOX::Abstract::MultiVector::DenseMatrix& Y1,
                const NOX::Abstract::MultiVector& Y2,
                const NOX::Abstract::MultiVector::DenseMatrix& inT,
                const NOX::Abstract::MultiVector& inA,
                NOX::Abstract::MultiVector& UU,
                NOX::Abstract::MultiVector& VV,
                bool use_jac_transpose)
{
  // Now compute U and V in P = J + U*V^T, U = A*Y_1 + J*Y_2, V = Y_2*T^T
  {
    const auto& Y2_tpetra = NOX::Tpetra::getTpetraMultiVector(Y2);
    auto& UU_tpetra = NOX::Tpetra::getTpetraMultiVector(UU);
    if (use_jac_transpose)
      tpetraOp->apply(Y2_tpetra,UU_tpetra,Teuchos::TRANS);
    else
      tpetraOp->apply(Y2_tpetra,UU_tpetra,Teuchos::NO_TRANS);
  }

  UU.update(Teuchos::NO_TRANS, 1.0, inA, Y1, 1.0);
  VV.update(Teuchos::TRANS, 1.0, Y2, inT, 0.0);

  return NOX::Abstract::Group::Ok;
}

void
LOCA::BorderedSolver::TpetraHouseholder::
updateCrsMatrixForPreconditioner(const NOX::Abstract::MultiVector& UU,
                                 const NOX::Abstract::MultiVector& VV,
                                 NOX::TCrsMatrix& matrix) const
{
  matrix.resumeFill();

  auto& UU_tpetra = NOX::Tpetra::getTpetraMultiVector(UU);
  auto& VV_tpetra = NOX::Tpetra::getTpetraMultiVector(VV);
  const auto uu = UU_tpetra.getLocalViewDevice(::Tpetra::Access::ReadOnly);
  const auto vv = VV_tpetra.getLocalViewDevice(::Tpetra::Access::ReadOnly);

  const auto numRows = matrix.getLocalNumRows();
  const auto rowMap = matrix.getRowMap()->getLocalMap();
  const auto colMap = matrix.getColMap()->getLocalMap();
  const auto uMap = UU_tpetra.getMap()->getLocalMap();
  const auto vMap = VV_tpetra.getMap()->getLocalMap();
  auto J_view = matrix.getLocalMatrixDevice();
  auto numConstraintsLocal = numConstraints; // for cuda lambda capture

  TEUCHOS_ASSERT(static_cast<size_t>(matrix.getRowMap()->getLocalNumElements()) == uu.extent(0));
  TEUCHOS_ASSERT(static_cast<size_t>(matrix.getRowMap()->getLocalNumElements()) == vv.extent(0));
  TEUCHOS_ASSERT(numConstraintsLocal == static_cast<int>(uu.extent(1)));
  TEUCHOS_ASSERT(numConstraintsLocal == static_cast<int>(vv.extent(1)));

  Kokkos::parallel_for("Add UV^T to M",Kokkos::RangePolicy<NOX::DeviceSpace>(0,numRows),KOKKOS_LAMBDA (const int row) {
    const GO row_gid = rowMap.getGlobalElement(row);
    const LO u_row_lid = uMap.getLocalElement(row_gid);
    auto rowView = J_view.row(row);

    const auto numEntries = rowView.length;
    for (int col=0; col<numEntries; ++col) {

      // Only included contributions from U*V^T on this proc
      const GO col_gid = colMap.getGlobalElement(rowView.colidx(col));
      int v_row_lid = vMap.getLocalElement(col_gid);
      if (v_row_lid != ::Tpetra::Details::OrdinalTraits<LO>::invalid()) {

        // val = sum_{k=0}^m U(i,k)*V(j,k)
        ST val = 0.0;
        for (int k=0; k<numConstraintsLocal; ++k)
          val += uu(u_row_lid,k) * vv(v_row_lid,k);

        // replace J(row,col) with J(row,col) + U*V^T
        rowView.value(col) += val;
      }
    }
  });
  Kokkos::fence();
  matrix.fillComplete();

  /*
  // Get number of rows on this processor
  int numMyRows = jac.NumMyRows();

  int numEntries;
  int U_LDA;
  int V_LDA;
  int *row_indices;
  double *row_values;
  double *U_values;
  double *V_values;
  double val;

  // Get pointers to U & V values
  const NOX::Epetra::MultiVector nox_epetra_U =
    dynamic_cast<const NOX::Epetra::MultiVector&>(UU);
  const Epetra_MultiVector& epetra_U = nox_epetra_U.getEpetraMultiVector();
  const NOX::Epetra::MultiVector& nox_epetra_V =
    dynamic_cast<const NOX::Epetra::MultiVector&>(VV);
  const Epetra_MultiVector& epetra_V = nox_epetra_V.getEpetraMultiVector();
  epetra_U.ExtractView(&U_values, &U_LDA);
  epetra_V.ExtractView(&V_values, &V_LDA);

  const Epetra_BlockMap& v_map = epetra_V.Map();
  const Epetra_BlockMap& u_map = epetra_V.Map();
  const Epetra_BlockMap& row_map = jac.RowMap();
  const Epetra_BlockMap& col_map = jac.ColMap();

  for (int row=0; row<numMyRows; row++) {

    // extract view of row
    jac.ExtractMyRowView(row, numEntries, row_values, row_indices);

    int row_gid = row_map.GID(row);
    int u_row_lid = u_map.LID(row_gid);

    for (int col=0; col<numEntries; col++) {

      // Only included contributions from U*V^T on this proc
      int col_gid = col_map.GID(row_indices[col]);
      if (v_map.MyGID(col_gid)) {
    int v_row_lid = v_map.LID(col_gid);

    // val = sum_{k=0}^m U(i,k)*V(j,k)
    val = 0;
    for (int k=0; k<numConstraints; k++)
      val += U_values[k*U_LDA+u_row_lid]*V_values[k*V_LDA+v_row_lid];

    // replace J(row,col) with J(row,col) + val
    row_values[col] += val;

      }

    }
  }
  */
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::BorderedSolver::TpetraHouseholder::createBlockMV(
                    const NOX::Abstract::MultiVector& v) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                             "ERROR: LOCA::BorderedSolver::TpetraHouseholder::createBlockMV - Support for COMPLEX systems not supported!");

  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);

  /*
#ifdef HAVE_NOX_EPETRAEXT
  const LOCA::Hopf::ComplexMultiVector& cv =
    dynamic_cast<const LOCA::Hopf::ComplexMultiVector&>(v);
  Teuchos::RCP<const NOX::Abstract::MultiVector> v_real =
    cv.getRealMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector> v_imag =
    cv.getImagMultiVec();
  Teuchos::RCP<const NOX::Epetra::MultiVector> nox_epetra_v_real =
    Teuchos::rcp_dynamic_cast<const NOX::Epetra::MultiVector>(v_real);
  Teuchos::RCP<const NOX::Epetra::MultiVector> nox_epetra_v_imag =
    Teuchos::rcp_dynamic_cast<const NOX::Epetra::MultiVector>(v_imag);

  Teuchos::RCP<EpetraExt::BlockMultiVector> epetra_v =
    Teuchos::rcp(new EpetraExt::BlockMultiVector(*baseMap, *globalMap,
                         v.numVectors()));
  epetra_v->LoadBlockValues(nox_epetra_v_real->getEpetraMultiVector(), 0);
  epetra_v->LoadBlockValues(nox_epetra_v_imag->getEpetraMultiVector(), 1);

  return Teuchos::rcp(new NOX::Epetra::MultiVector(
                      epetra_v,
                      NOX::DeepCopy,
                      NOX::Epetra::MultiVector::CreateView));
#else
  globalData->locaErrorCheck->throwError("LOCA::BorderedSolver::TpetraHouseholder::createBlockMV()",
                     "Must have EpetraExt support for Hopf tracking.  Configure trilinos with --enable-epetraext");
  return Teuchos::null;
#endif
  */
}

void
LOCA::BorderedSolver::TpetraHouseholder::setBlockMV(
                       const NOX::Abstract::MultiVector& bv,
                       NOX::Abstract::MultiVector& v) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                             "ERROR: LOCA::BorderedSolver::TpetraHouseholder::setBlockMV - Support for COMPLEX systems not supported!");

  /*
#ifdef HAVE_NOX_EPETRAEXT
  LOCA::Hopf::ComplexMultiVector& cv =
    dynamic_cast<LOCA::Hopf::ComplexMultiVector&>(v);
  Teuchos::RCP<NOX::Abstract::MultiVector> v_real =
    cv.getRealMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector> v_imag =
    cv.getImagMultiVec();
  Teuchos::RCP<NOX::Epetra::MultiVector> nox_epetra_v_real =
    Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(v_real);
  Teuchos::RCP<NOX::Epetra::MultiVector> nox_epetra_v_imag =
    Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(v_imag);

  const NOX::Epetra::MultiVector& nox_epetra_bv =
    dynamic_cast<const NOX::Epetra::MultiVector&>(bv);
  const Epetra_MultiVector& epetra_bv = nox_epetra_bv.getEpetraMultiVector();
  const EpetraExt::BlockMultiVector& block_bv =
    dynamic_cast<const EpetraExt::BlockMultiVector&>(epetra_bv);

  block_bv.ExtractBlockValues(nox_epetra_v_real->getEpetraMultiVector(), 0);
  block_bv.ExtractBlockValues(nox_epetra_v_imag->getEpetraMultiVector(), 1);
#else
  globalData->locaErrorCheck->throwError("LOCA::BorderedSolver::TpetraHouseholder::setBlockMV()",
                     "Must have EpetraExt support for Hopf tracking.  Configure trilinos with --enable-epetraext");
#endif
  */
}
