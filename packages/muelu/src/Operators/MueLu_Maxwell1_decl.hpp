// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_MAXWELL1_DECL_HPP
#define MUELU_MAXWELL1_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"

#include "MueLu_ReitzingerPFactory_fwd.hpp"

#include "MueLu_Utilities_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_Hierarchy_fwd.hpp"
#include "MueLu_RAPFactory_fwd.hpp"
#include "MueLu_RebalanceAcFactory_fwd.hpp"
#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_SmootherBase_fwd.hpp"

#include "Xpetra_Map_fwd.hpp"
#include "Xpetra_Matrix_fwd.hpp"
#include "Xpetra_MatrixFactory_fwd.hpp"
#include "Xpetra_MultiVectorFactory_fwd.hpp"
#include "Xpetra_VectorFactory_fwd.hpp"
#include "Xpetra_CrsMatrixWrap_fwd.hpp"

namespace MueLu {

/*!
  @brief Preconditioner (wrapped as a Xpetra::Operator) for Maxwell's equations in curl-curl form.


  @ingroup MueLuAdapters
*/
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class Maxwell1 : public VerboseObject, public Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_MAXWELL1_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;
  typedef typename Teuchos::ScalarTraits<Scalar>::coordinateType coordinateType;
  typedef typename Xpetra::MultiVector<coordinateType, LO, GO, NO> RealValuedMultiVector;

  //! Constructor
  Maxwell1()
    : Hierarchy11_(Teuchos::null)
    , Hierarchy22_(Teuchos::null)
    , HierarchyGmhd_(Teuchos::null)
    , mode_(MODE_STANDARD) {
  }

  //! Constructor with Hierarchies
  Maxwell1(Teuchos::RCP<Hierarchy> H11, Teuchos::RCP<Hierarchy> H22)
    : Hierarchy11_(H11)
    , Hierarchy22_(H22)
    , HierarchyGmhd_(Teuchos::null)
    , mode_(MODE_STANDARD) {
  }

  /** Constructor with Jacobian
   *
   * \param[in] SM_Matrix Jacobian
   * \param[in] D0_Matrix Discrete Gradient
   * \param[in] Nullspace Null space (needed for periodic)
   * \param[in] Coords Nodal coordinates
   * \param[in] List Parameter list
   * \param[in] ComputePrec If true, compute the preconditioner immediately
   */
  Maxwell1(const Teuchos::RCP<Matrix>& SM_Matrix,
           const Teuchos::RCP<Matrix>& D0_Matrix,
           const Teuchos::RCP<MultiVector>& Nullspace,
           const Teuchos::RCP<RealValuedMultiVector>& Coords,
           Teuchos::ParameterList& List,
           bool ComputePrec = true)
    : mode_(MODE_STANDARD) {
    RCP<Matrix> Kn_Matrix;
    initialize(D0_Matrix, Kn_Matrix, Nullspace, Coords, List);
    resetMatrix(SM_Matrix, ComputePrec);
  }

  /** Constructor with Jacobian and nodal matrix
   *
   * \param[in] SM_Matrix Jacobian
   * \param[in] D0_Matrix Discrete Gradient
   * \param[in] Kn_Matrix Nodal Laplacian
   * \param[in] Coords Nodal coordinates
   * \param[in] List Parameter list
   * \param[in] ComputePrec If true, compute the preconditioner immediately
   */
  Maxwell1(const Teuchos::RCP<Matrix>& SM_Matrix,
           const Teuchos::RCP<Matrix>& D0_Matrix,
           const Teuchos::RCP<Matrix>& Kn_Matrix,
           const Teuchos::RCP<MultiVector>& Nullspace,
           const Teuchos::RCP<RealValuedMultiVector>& Coords,
           Teuchos::ParameterList& List,
           bool ComputePrec = true)
    : mode_(MODE_STANDARD) {
    initialize(D0_Matrix, Kn_Matrix, Nullspace, Coords, List);
    resetMatrix(SM_Matrix, ComputePrec);
  }

  /** Gmhd GMHD Constructor with Jacobian and nodal matrix AND Gmhd matrix
   *
   * \param[in] SM_Matrix Jacobian
   * \param[in] D0_Matrix Discrete Gradient
   * \param[in] Kn_Matrix Nodal Laplacian
   * \param[in] Coords Nodal coordinates
   * \param[in] List Parameter list
   * \param[in] GmhdA_Matrix Gmhd matrix including generalized Ohms law equations
   * \param[in] ComputePrec If true, compute the preconditioner immediately
   */
  Maxwell1(const Teuchos::RCP<Matrix>& SM_Matrix,
           const Teuchos::RCP<Matrix>& D0_Matrix,
           const Teuchos::RCP<Matrix>& Kn_Matrix,
           const Teuchos::RCP<MultiVector>& Nullspace,
           const Teuchos::RCP<RealValuedMultiVector>& Coords,
           Teuchos::ParameterList& List, const Teuchos::RCP<Matrix>& GmhdA_Matrix,
           bool ComputePrec = true)
    : mode_(MODE_GMHD_STANDARD) {
    initialize(D0_Matrix, Kn_Matrix, Nullspace, Coords, List);
    resetMatrix(SM_Matrix, ComputePrec);
    GmhdA_Matrix_  = GmhdA_Matrix;
    HierarchyGmhd_ = rcp(new Hierarchy("HierarchyGmhd"));
    GMHDSetupHierarchy(List);
  }

  /** Constructor with parameter list
   *
   * \param[in] SM_Matrix Jacobian
   * \param[in] List Parameter list
   * \param[in] ComputePrec If true, compute the preconditioner immediately
   */
  Maxwell1(const Teuchos::RCP<Matrix>& SM_Matrix,
           Teuchos::ParameterList& List,
           bool ComputePrec = true)
    : mode_(MODE_STANDARD) {
    RCP<MultiVector> Nullspace        = List.get<RCP<MultiVector> >("Nullspace", Teuchos::null);
    RCP<RealValuedMultiVector> Coords = List.get<RCP<RealValuedMultiVector> >("Coordinates", Teuchos::null);
    RCP<Matrix> D0_Matrix             = List.get<RCP<Matrix> >("D0");
    RCP<Matrix> Kn_Matrix;
    if (List.isType<RCP<Matrix> >("Kn"))
      Kn_Matrix = List.get<RCP<Matrix> >("Kn");

    initialize(D0_Matrix, Kn_Matrix, Nullspace, Coords, List);

    if (SM_Matrix != Teuchos::null)
      resetMatrix(SM_Matrix, ComputePrec);
  }

  //! Destructor.
  virtual ~Maxwell1() {}

  //! Returns the Xpetra::Map object associated with the domain of this operator.
  const Teuchos::RCP<const Map> getDomainMap() const;

  //! Returns the Xpetra::Map object associated with the range of this operator.
  const Teuchos::RCP<const Map> getRangeMap() const;

  //! Returns Jacobian matrix SM
  const Teuchos::RCP<Matrix>& getJacobian() const {
    return SM_Matrix_;
  }

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

 private:
  //! Generates the Kn matrix
  Teuchos::RCP<Matrix> generate_kn() const;

  //! Sets up hiearchy for GMHD matrices that include generalized Ohms law equations
  void GMHDSetupHierarchy(Teuchos::ParameterList& List) const;

  /** Initialize with matrices except the Jacobian (don't compute the preconditioner)
   *
   * \param[in] D0_Matrix Discrete Gradient
   * \param[in] Kn_Matrix Kn nodal matrix
   * \param[in] Nullspace Null space (needed for periodic)
   * \param[in] Coords Nodal coordinates
   * \param[in] List Parameter list
   */
  void initialize(const Teuchos::RCP<Matrix>& D0_Matrix,
                  const Teuchos::RCP<Matrix>& Kn_Matrix,
                  const Teuchos::RCP<MultiVector>& Nullspace,
                  const Teuchos::RCP<RealValuedMultiVector>& Coords,
                  Teuchos::ParameterList& List);

  //! apply RefMaxwell additive 2x2 style cycle
  void applyInverseRefMaxwellAdditive(const MultiVector& RHS, MultiVector& X) const;

  //! apply standard Maxwell1 cycle
  void applyInverseStandard(const MultiVector& RHS, MultiVector& X) const;

  //! allocate multivectors for solve
  void allocateMemory(int numVectors) const;

  //! dump out matrix
  void dump(const Matrix& A, std::string name) const;

  //! dump out multivector
  void dump(const MultiVector& X, std::string name) const;

  //! dump out real-valued multivector
  void dumpCoords(const RealValuedMultiVector& X, std::string name) const;

  //! dump out boolean ArrayView
  void dump(const Teuchos::ArrayRCP<bool>& v, std::string name) const;

  //! dump out boolean Kokkos::View
  void dump(const Kokkos::View<bool*, typename Node::device_type>& v, std::string name) const;

  //! get a (synced) timer
  Teuchos::RCP<Teuchos::TimeMonitor> getTimer(std::string name, RCP<const Teuchos::Comm<int> > comm = Teuchos::null) const;

  //! ParameterLists
  mutable Teuchos::ParameterList parameterList_, precList11_, precList22_;

  //! Two hierarchies: one for the (1,1)-block, another for the (2,2)-block
  Teuchos::RCP<Hierarchy> Hierarchy11_, Hierarchy22_, HierarchyGmhd_;

  //! Various matrices
  Teuchos::RCP<Matrix> SM_Matrix_, D0_Matrix_, Kn_Matrix_, GmhdA_Matrix_;

  //! Vectors for BCs
  Kokkos::View<bool*, typename Node::device_type::memory_space> BCrowsKokkos_, BCcolsKokkos_, BCdomainKokkos_;
  int BCedges_, BCnodes_;
  Teuchos::ArrayRCP<bool> BCrows_, BCcols_, BCdomain_;
  //! Nullspace
  Teuchos::RCP<MultiVector> Nullspace_;
  //! Coordinates
  Teuchos::RCP<RealValuedMultiVector> Coords_;
  //! Some options
  bool useKokkos_, allEdgesBoundary_, allNodesBoundary_, dump_matrices_, enable_reuse_, syncTimers_;
  bool applyBCsTo22_;

  //! Execution modes
  typedef enum { MODE_STANDARD = 0,
                 MODE_REFMAXWELL,
                 MODE_EDGE_ONLY,
                 MODE_GMHD_STANDARD } mode_type;
  mode_type mode_;

  //! Temporary memory (cached vectors for RefMaxwell-style)
  RCP<Matrix> P11_;
  mutable Teuchos::RCP<MultiVector> residualFine_, residual11c_, residual22_, update11c_, update22_;
};

}  // namespace MueLu

#define MUELU_MAXWELL1_SHORT
#endif  // MUELU_MAXWELL1_DECL_HPP
