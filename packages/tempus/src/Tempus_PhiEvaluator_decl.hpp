//@HEADER
// *****************************************************************************
// TODO
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluator_decl_hpp
#define Tempus_PhiEvaluator_decl_hpp

#include "Tempus_config.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"

#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_ModelEvaluator.hpp"

#include "Teuchos_DefaultComm.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Teuchos_CommHelpers.hpp"

namespace Tempus {

enum class PhiInitialization {
    ONLY_MASS,
    JACOBIAN_AND_MASS,
};

/** \brief PhiLinearSolver
 * Uses the ModelEvaluator to compute and hold the Mass matrix and Jacobian
 */
template <class Scalar>
class PhiLinearSolver {
 public:
  PhiLinearSolver(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>> appModel, bool lumpMass=false)
    : appModel_(appModel), lumpMass_(lumpMass) {
  }

  ~PhiLinearSolver() {}

  void setLumpMassMatrix(const bool lump);
  void clearMemory();

  /** \brief Set the eigensolver parameter list used by computeJacobianSpectrumBounds.
   *
   *  The ParameterList should contain the entries defined in
   *  PhiEvaluator::getValidParametersBasic() under the "Eigensolver" sublist.
   *  When null, all eigensolver parameters fall back to compiled-in defaults.
   */
  void setEigensolverParams(Teuchos::RCP<Teuchos::ParameterList> pl);

  /** \brief Initialize PhiSolver
   *
   *  This function will check if mass matrix and Jacobian have been computed.
   *  This function does not make member data consistent, but just checks it.
   *  This ensures it is inexpensive.
   */
  /// Return if PhiSolver is initialized.
  void checkInitialized(const PhiInitialization& mode) const;

  void computeMassMatrix(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs);
  void applyMass(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> Mf,
                 const Teuchos::RCP<const Thyra::VectorBase<Scalar>> f) const;
  void solveMass(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> f,
                 const Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mf) const;

  void computeJacobian(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs);
  void applyJacobian(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> Jf,
                     const Teuchos::RCP<const Thyra::VectorBase<Scalar>> f) const;

  // build and return a dt-scaled LinOp
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> buildL(const Scalar dt) const;

  // TODO: make that one public function
  // build and return extended LinOp
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> buildATilde(
      const Scalar dt,
      const Teuchos::ArrayView<const Teuchos::RCP<const Thyra::VectorBase<Scalar>>> &rhs_B);

  // build and return the right hand side for the extended LinOp
  Teuchos::RCP<Thyra::ProductVectorBase<Scalar>> buildv(
      const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> space,
      const Teuchos::RCP<const Thyra::VectorBase<Scalar>> x0) const;

  /** \brief Compute the extrema of the the scaled system Jacobian.
  *
  *   Computes the minimum real, maximum real and maximum imaginary components
  *   of the Jacobian spectrum using Anasazi block Krylov-Schur.
  *   These values maybe used to set hyperparameters of the PhiEvaluators.
  *
   @param a    The minimum real spectrum bound
   @param b    The maximum real spectrum bound
   @param c    The maximum imaginary spectrum bound
  */
  void computeJacobianSpectrumBounds(double& a, double& b, double& c);

  // Solve Mass plus Jacobian, for given inArgs (recompute matrices from from ModelEvaluator)
  Thyra::SolveStatus<Scalar> assembleAndsolveMpJ(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
                                                 const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
                                                 const Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mf,
                                                 Scalar alpha = 1., Scalar beta = 0.) const;

  // Solve Mass plus Jacobian, for precomputed Mass and Jacobian
  Thyra::SolveStatus<Scalar> solveMpJ(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
                                      const Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mf,
                                      Scalar alpha=1., Scalar beta=0.) const;

 private:
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > appModel_;
  bool lumpMass_;
  bool isInitialized_;

  // massMatrix_ and inverseMassMatrix_ are either lumped, or not, depending on bool lumpMass_
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> massMatrix_;
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> inverseMassMatrix_;

  Teuchos::RCP<Thyra::LinearOpBase<Scalar>> jacobianMatrix_;

  // internal variables for extended matrix strategy
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> Atilde_;
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> KMatrix_;
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> bMatrix_;

  // internal storage for eigenvectors of the M^{-1}*J LinOp
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar>> evecsMinvJ_;

  // eigensolver parameters forwarded from the owning PhiEvaluator
  Teuchos::RCP<Teuchos::ParameterList> eigensolverPL_;
  // internal methods for the ATilde LinOp
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> buildK(const Thyra::Ordinal max_phi_order);
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> buildB(
      const Teuchos::ArrayView<const Teuchos::RCP<const Thyra::VectorBase<Scalar>>> &rhs_B);
};


/** \brief PhiEvaluator evaluates
 *
 *  \f$[x = \varphi_k(J) b]\f$, where
 *
 *   - b is a right hand side vector
 *   - J is a linear operator
 */
template <class Scalar>
class PhiEvaluator
  : virtual public Teuchos::Describable,
    virtual public Teuchos::VerboseObject<PhiEvaluator<Scalar> > {
 public:
  PhiEvaluator();

  PhiEvaluator(
      std::string name);

  ~PhiEvaluator() {}

  /// \name Basic PhiEvaluator Methods
  //@{

  /// Make a shallow copy of PhiEvaluator (i.e., only RCPs).
  void copy(Teuchos::RCP<const PhiEvaluator<Scalar> > sh);
  //@}

  /// \name Accessor methods
  //@{
  /// Get this PhiEvaluator's name
  std::string getName() const { return name_; }

  /// Set this PhiEvaluator's name
  void setName(std::string name) { name_ = name; }

  /// Return a valid ParameterList with current settings.
  virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  /// Return a valid ParameterList with basic settings.
  Teuchos::RCP<Teuchos::ParameterList> getValidParametersBasic() const;

  /// Return a valid non-const ParameterList with current settings.
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();

  /// Set the parameters from a ParameterList
  virtual void setPhiEvaluatorValues(Teuchos::RCP<Teuchos::ParameterList> pl);

  /// \name Overridden from Teuchos::Describable
  //@{
  virtual std::string description() const;
  virtual void describe(Teuchos::FancyOStream& out,
                        const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  /** \brief Initialize PhiEvaluator
   *
   *  This function will check if all member data is initialized
   *  and is consistent.  This function does not make member data
   *  consistent, but just checks it.  This ensures it is inexpensive.
   */
  void initialize();

  /// Return if PhiEvaluator is initialized.
  bool isInitialized() const { return isInitialized_; }

  void checkInitialized() const;

  void setLumpMassMatrix(bool lumpMassMatrix);
  void setConstantMassMatrix(bool constantMassMatrix);

  /// set the ModelEvaluator
  void setModel(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > appModel);

  /** \brief   Set the linearization point for the Jacobian calculation
   *
   *  The linearization point x and time t are taken from inArgs.
   *  This computes the Mass and Jacobian matrix for future use.
   */
  void setLinearizationPoint(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
                             const PhiInitialization& mode = PhiInitialization::JACOBIAN_AND_MASS);

  /** \brief   Adapt internal PhiEvaluator hyperparameters to the current Jacobian.
   *
   *  Called after setLinearizationPoint.  Default implementation is a null-op.
   *  This method should be overridden by inheriting PhiEvaluator implementations when
   *  there are hyperparameters of the method that require analysis of the Jacobian.
   */
  virtual void adaptEvaluator() { this->isInitialized(); };

  /** \brief  Compute the Phi_k function of cdt*Jacobian for right hand side Mrhs_b
   *
   *  phi_order is the index of the phi-function Phi_k.
   *  For an implicit model, the right hand side contains the mass matrix M,
   *  which is solved as part of this method.
   */
  virtual Thyra::SolveStatus<Scalar> computePhi(
      const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
      const int phi_order, const Scalar cdt,
      const Teuchos::RCP<const Thyra::VectorBase<Scalar>> &Mrhs_b);

  /** \brief  Compute the Phi function of cdt times Jacobian for a linear combination of vectors Mrhs_B
   *
   *  The vectors in Mrhs_B are at the index of the vector corresponding to the phi_order of the
   *  respective Phi function, Mrhs_b[0] is the rhs for the matrix exponential, Mrhs_b[1] is the
   *  right-hand side for the phi_1 function. A Teuchos::null RCP is interpreted as a zero vector.
   *  For an implicit model, the right hand side contains a multiplication with the mass matrix M,
   *  which is solved as part of this method.
   */
  virtual Thyra::SolveStatus<Scalar> computePhis(
      const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
      const Scalar cdt,
      const Teuchos::ArrayView<const Teuchos::RCP<const Thyra::VectorBase<Scalar>>> &Mrhs_B);

  // Multiply the mass matrix (lumped or not) with right hand side f
  void applyMass(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> Mf,
                 const Teuchos::RCP<const Thyra::VectorBase<Scalar>> f) const;

  // Invert the mass matrix (lumped or not) with right hand side Mf
  void solveMass(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> f,
                 const Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mf) const;

  // Multiply the MassJacobian matrix with right hand side Mf
  void applyJacobian(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> MJf,
                     const Teuchos::RCP<const Thyra::VectorBase<Scalar>> f) const;

 protected:
  std::string name_;
  bool lumpMassMatrix_;
  bool constantMassMatrix_;
  bool useAtildeForSingleRHS_;

  bool isInitialized_;  ///< Bool if PhiEvaluator is initialized.

  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>> appModel_;
  Teuchos::RCP<Tempus::PhiLinearSolver<Scalar>> phiLinSolv_;

  /// Eigenvalue solver settings
  Teuchos::RCP<Teuchos::ParameterList> eigensolverPL_;

  //mutable
  Teuchos::RCP<const Thyra::ModelEvaluatorBase::InArgs<Scalar>> inArgs_lin_;

  /** \brief  Internal method for a LinOp, used for default impl. of computePhi/computePhis
   *
   *  Computes v := phi_{phi_order}(L)v in place (overwriting the rhs with the result)
   */
  virtual Thyra::SolveStatus<Scalar> computeLinOpPhi(
      const int phi_order,
      const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> L,
      const Teuchos::Ptr<Thyra::VectorBase<Scalar>> v,
      const Scalar cdt=1.0
    ) = 0;
};

/// Nonmember helper to convert a Thyra LinearOp to a SerialDenseMatrix
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::SerialDenseMatrix<int, Scalar>
operatorToDense(Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> lop)
{
  const int numCols = static_cast<int>(lop->domain()->dim());
  const int numRows = static_cast<int>(lop->range()->dim());

  auto rangeSpmd =
    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<Scalar>>(
        lop->range(), true);

  auto domainSpmd =
    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<Scalar>>(
        lop->domain(), true);

  auto comm = rangeSpmd->getComm();
  const int rank = comm.is_null() ? 0 : comm->getRank();
  const int numProcs = comm.is_null() ? 1 : comm->getSize();

  const int rowOffset   = static_cast<int>(rangeSpmd->localOffset());
  const int rowLocalDim = static_cast<int>(rangeSpmd->localSubDim());

  const int colOffset   = static_cast<int>(domainSpmd->localOffset());
  const int colLocalDim = static_cast<int>(domainSpmd->localSubDim());

  Teuchos::RCP<Thyra::MultiVectorBase<Scalar>> Id =
      Thyra::createMembers(lop->domain(), numCols);

  Thyra::assign(Id.ptr(), Scalar(0));

  // Build distributed identity. Each rank sets only the identity entries
  // whose domain rows it owns.
  for (int jLocal = 0; jLocal < colLocalDim; ++jLocal) {
    const int jGlobal = colOffset + jLocal;

    Teuchos::RCP<Thyra::VectorBase<Scalar>> col_j = Id->col(jGlobal);
    Thyra::set_ele(jGlobal, Scalar(1), col_j.ptr());
  }

  Teuchos::RCP<Thyra::MultiVectorBase<Scalar>> Y =
      Thyra::createMembers(lop->range(), numCols);

  Thyra::apply(*lop, Thyra::NOTRANS, *Id, Y.ptr(), Scalar(1), Scalar(0));

  // Local contribution to the dense matrix.
  // Column-major indexing: i + j*numRows.
  std::vector<Scalar> localDense(numRows * numCols, Scalar(0));
  std::vector<Scalar> globalDense(numRows * numCols, Scalar(0));

  for (int j = 0; j < numCols; ++j) {
    Teuchos::RCP<const Thyra::VectorBase<Scalar>> colY = Y->col(j);

    if (rowLocalDim > 0) {
      Thyra::Range1D localRange(rowOffset, rowOffset + rowLocalDim - 1);

      Thyra::ConstDetachedVectorView<Scalar> localView(*colY, localRange);

      for (int iLocal = 0; iLocal < rowLocalDim; ++iLocal) {
        const int iGlobal = rowOffset + iLocal;
        localDense[iGlobal + j * numRows] = localView[iGlobal];
      }
    }
  }

  // Combine local pieces. Since each row is owned by one rank,
  // sum-reduction reconstructs the full dense matrix on every rank.
  // Ignore the logic if comm is null.
  if (comm.is_null()) {
    globalDense = localDense;
  }
  else {
    Teuchos::reduceAll(
        *comm,
        Teuchos::REDUCE_SUM,
        static_cast<long int>(localDense.size()),
        localDense.data(),
        globalDense.data());
  }

  Teuchos::SerialDenseMatrix<int, Scalar> globalDenseMat(numRows, numCols);

  for (int j = 0; j < numCols; ++j) {
    for (int i = 0; i < numRows; ++i) {
      globalDenseMat(i, j) = globalDense[i + j * numRows];
    }
  }

  return globalDenseMat;
}

/// Nonmember helper to compute eigenvalues of a SerialDenseMatrix
// ------------------------------------------------------------------------
template<typename OrdinalType, typename Scalar>
int denseEigenvalues(const Teuchos::SerialDenseMatrix<OrdinalType, Scalar>& A,
                           Teuchos::Array<Scalar>& eigs_re,
                           Teuchos::Array<Scalar>& eigs_im) {
    OrdinalType n = A.numRows();
    Teuchos::SerialDenseMatrix<OrdinalType, Scalar> A_copy(Teuchos::Copy, A);
    Teuchos::LAPACK<OrdinalType, Scalar> lapack;
    // Quick workspace query
    Scalar lworkQuery = 0.0;
    OrdinalType info = 0;
    lapack.GEEV('N', 'N', n, A_copy.values(), A_copy.stride(), eigs_re.getRawPtr(), eigs_im.getRawPtr(),
                nullptr, 1, nullptr, 1, &lworkQuery, -1, &info);
    // Solve
    OrdinalType lwork = static_cast<OrdinalType>(lworkQuery);
    Teuchos::Array<Scalar> work(lwork);
    lapack.GEEV('N', 'N', n, A_copy.values(), A_copy.stride(), eigs_re.getRawPtr(), eigs_im.getRawPtr(),
                nullptr, 1, nullptr, 1, work.getRawPtr(), lwork, &info);
    return info;
}

}  // namespace Tempus

#endif  // Tempus_PhiEvaluator_decl_hpp
