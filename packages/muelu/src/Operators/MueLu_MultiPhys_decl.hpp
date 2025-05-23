// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_MULTIPHYS_DECL_HPP
#define MUELU_MULTIPHYS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"

#include "MueLu_SaPFactory_fwd.hpp"

#include "MueLu_SmootherFactory_fwd.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Utilities_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_Hierarchy_fwd.hpp"
#include "MueLu_RAPFactory_fwd.hpp"
#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_SmootherBase.hpp"

#include "Xpetra_Map_fwd.hpp"
#include "Xpetra_Matrix_fwd.hpp"
#include "Xpetra_MatrixFactory_fwd.hpp"
#include "Xpetra_MultiVectorFactory_fwd.hpp"
#include "Xpetra_VectorFactory_fwd.hpp"
#include "Xpetra_CrsMatrixWrap_fwd.hpp"

namespace MueLu {

/*!
  @brief Preconditioner (wrapped as a Xpetra::Operator) for solving MultiPhysics PDEs.


  @ingroup MueLuAdapters
*/
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class MultiPhys : public VerboseObject, public Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_MULTIPHYS_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;
  typedef typename Teuchos::ScalarTraits<Scalar>::coordinateType coordinateType;
  typedef typename Xpetra::MultiVector<coordinateType, LO, GO, NO> RealValuedMultiVector;

  //! Constructor
  MultiPhys()
    : AmatMultiphysics_(Teuchos::null)
    , hierarchyMultiphysics_(Teuchos::null)
    , nBlks_(0) {
  }

  /** Constructor
   *
   * \param[in] AmatMultiPhysics      Multiphysics discretization matrix
   * \param[in] arrayOfAuxMatrices    Auxiliary matrices used to generate subblock prolongators for multiphysics system
   * \param[in] arrayOfNullspaces     Nullspace multivectors used to generate subblock prolongators for multiphysics system
   * \param[in] arrayOfCoords         Coordinate multivectors used to generate subblock prolongators for multiphysics system
   * \param[in] nBlks                 nBlks x nBlks gives the block dimensions of the multiphysics operator
   * \param[in] List Parameter list
   * \param[in] ComputePrec If true, compute the preconditioner immediately
   * \param[in] arrayOfMaterials      Material multivectors used to generate subblock prolongators for multiphysics system
   * \param[in] OmitSubblockSmoother If true, omit construction of subblock-level smoothers
   */
  MultiPhys(const Teuchos::RCP<Matrix>& AmatMultiPhysics,
            const Teuchos::ArrayRCP<RCP<Matrix>> arrayOfAuxMatrices,
            const Teuchos::ArrayRCP<Teuchos::RCP<MultiVector>> arrayOfNullspaces,
            const Teuchos::ArrayRCP<Teuchos::RCP<RealValuedMultiVector>> arrayOfCoords,
            const int nBlks,
            Teuchos::ParameterList& List,
            bool ComputePrec                                                    = true,
            const Teuchos::ArrayRCP<Teuchos::RCP<MultiVector>> arrayOfMaterials = Teuchos::null,
            bool OmitSubblockSmoother                                           = true)
    : AmatMultiphysics_(AmatMultiPhysics)
    , arrayOfAuxMatrices_(arrayOfAuxMatrices)
    , arrayOfNullspaces_(arrayOfNullspaces)
    , arrayOfCoords_(arrayOfCoords)
    , arrayOfMaterials_(arrayOfMaterials)
    , OmitSubblockSmoother_(OmitSubblockSmoother)
    , nBlks_(nBlks) {
    initialize(AmatMultiPhysics, arrayOfAuxMatrices, arrayOfNullspaces, arrayOfCoords, nBlks, List, arrayOfMaterials);
    compute(false);
  }

  //! Destructor.
  virtual ~MultiPhys() {}

  //! Returns the Xpetra::Map object associated with the domain of this operator.
  const Teuchos::RCP<const Map> getDomainMap() const;

  //! Returns the Xpetra::Map object associated with the range of this operator.
  const Teuchos::RCP<const Map> getRangeMap() const;

  //! Set parameters
  void setParameters(Teuchos::ParameterList& list);

  //! Setup the preconditioner
  void compute(bool reuse = false);

  //! Reset system matrix
  void resetMatrix(Teuchos::RCP<Matrix> SM_Matrix_new, bool ComputePrec = true);

  //! Returns in Y the result of a Xpetra::Operator applied to a Xpetra::MultiVector X.
  //! \param[in]  X - MultiVector of dimension NumVectors to multiply with matrix.
  //! \param[out] Y - MultiVector of dimension NumVectors containing result.
  void apply(const MultiVector& X, MultiVector& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS,
             Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
             Scalar beta           = Teuchos::ScalarTraits<Scalar>::zero()) const;

  //! Indicates whether this operator supports applying the adjoint operator.
  bool hasTransposeApply() const;

  void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_HIGH) const;

  //! Compute a residual R = B - (*this) * X
  void residual(const MultiVector& X,
                const MultiVector& B,
                MultiVector& R) const {
    using STS = Teuchos::ScalarTraits<Scalar>;
    R.update(STS::one(), B, STS::zero());
    this->apply(X, R, Teuchos::NO_TRANS, -STS::one(), STS::one());
  }

  //! Return the multiphysics hierarchy
  Teuchos::RCP<Hierarchy> multiphysicsHierarchy() {
    return hierarchyMultiphysics_;
  }

  //! Return an array of hierarchies corresponding to each diagonal subblock
  Teuchos::ArrayRCP<Teuchos::RCP<Hierarchy>> subblockHierarchies() {
    return arrayOfHierarchies_;
  }

 private:
  /** Initialize with matrices except the Jacobian (don't compute the preconditioner)
   *
   * \param[in] AmatMultiPhysics      Multiphysics discretization matrix
   * \param[in] arrayOfAuxMatrices    Array of auxiliary matrices used to generate subblock prolongators for multiphysics system
   * \param[in] arrayOfNullspaces     Array of nullspace  multivectors used to generate subblock prolongators for multiphysics system
   * \param[in] arrayOfCoords         Array of coordinate multivectors used to generate subblock prolongators for multiphysics system
   * \param[in] nBlks                 nBlks x nBlks gives the block dimensions of the multiphysics operator
   * \param[in] List Parameter list
   * \param[in] arrayOfMaterials      Array of material multivectors used to generate subblock prolongators for multiphysics system
   */
  void initialize(const Teuchos::RCP<Matrix>& AmatMultiPhysics,
                  const Teuchos::ArrayRCP<RCP<Matrix>> arrayOfAuxMatrices,
                  const Teuchos::ArrayRCP<Teuchos::RCP<MultiVector>> arrayOfNullspaces,
                  const Teuchos::ArrayRCP<Teuchos::RCP<RealValuedMultiVector>> arrayOfCoords,
                  const int nBlks,
                  Teuchos::ParameterList& List,
                  const Teuchos::ArrayRCP<Teuchos::RCP<MultiVector>> arrayOfMaterials);

  //! apply standard MultiPhys cycle
  void applyInverse(const MultiVector& RHS, MultiVector& X) const;

  //! allocate multivectors for solve
  void allocateMemory(int numVectors) const;

  //! get a (synced) timer
  Teuchos::RCP<Teuchos::TimeMonitor> getTimer(std::string name, RCP<const Teuchos::Comm<int>> comm = Teuchos::null) const;

  //! ParameterLists
  mutable Teuchos::ParameterList parameterList_;

  //! Hierarchies: used to define P for (0,0)-block, .... (nBlks_-1,nBlks_-1) block

  Teuchos::RCP<Matrix> AmatMultiphysics_;                       // multiphysics discretization matrix
  Teuchos::RCP<Teuchos::ParameterList> paramListMultiphysics_;  // array of parameter lists directing MueLu's construct of subblock P operators
  Teuchos::RCP<Hierarchy> hierarchyMultiphysics_;               // multiphysics discretization matrix

  Teuchos::ArrayRCP<Teuchos::RCP<Teuchos::ParameterList>> arrayOfParamLists_;  // array of parameter lists directing MueLu's construct of subblock P operators
  Teuchos::ArrayRCP<Teuchos::RCP<Hierarchy>> arrayOfHierarchies_;
  Teuchos::ArrayRCP<Teuchos::RCP<Matrix>> arrayOfAuxMatrices_;            // array of discretization/auxiliary matrices used to generate subblock prolongators
  Teuchos::ArrayRCP<Teuchos::RCP<MultiVector>> arrayOfNullspaces_;        // array of nullspaces for smoothed aggregation.
  Teuchos::ArrayRCP<Teuchos::RCP<RealValuedMultiVector>> arrayOfCoords_;  // array of coordinates for smoothed aggregation/rebalancing.
  Teuchos::ArrayRCP<Teuchos::RCP<MultiVector>> arrayOfMaterials_;         // array of materials for smoothed aggregation.

  bool OmitSubblockSmoother_;

  int nBlks_;  // number of PDE sub-systems within multiphysics system
  bool useKokkos_, enable_reuse_, syncTimers_;
};

}  // namespace MueLu

#define MUELU_MULTIPHYS_SHORT
#endif  // MUELU_MULTIPHYS_DECL_HPP
