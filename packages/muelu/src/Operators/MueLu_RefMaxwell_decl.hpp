// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_REFMAXWELL_DECL_HPP
#define MUELU_REFMAXWELL_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"

#include "MueLu_CoalesceDropFactory_fwd.hpp"
#include "MueLu_CoarseMapFactory_fwd.hpp"
#include "MueLu_CoordinatesTransferFactory_fwd.hpp"
#include "MueLu_TentativePFactory_fwd.hpp"
#include "MueLu_SaPFactory_fwd.hpp"
#include "MueLu_UncoupledAggregationFactory_fwd.hpp"
#include "MueLu_AggregationExportFactory_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"

#include "MueLu_MasterList.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_SmootherBase_fwd.hpp"
#include "MueLu_SmootherPrototype_fwd.hpp"
#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_RAPFactory_fwd.hpp"

#include "MueLu_CoalesceDropFactory_kokkos_fwd.hpp"
#include "MueLu_TentativePFactory_kokkos_fwd.hpp"
#include "MueLu_SaPFactory_kokkos_fwd.hpp"
#include "MueLu_UncoupledAggregationFactory_kokkos_fwd.hpp"

#include "MueLu_ZoltanInterface_fwd.hpp"
#include "MueLu_Zoltan2Interface_fwd.hpp"
#include "MueLu_RepartitionHeuristicFactory_fwd.hpp"
#include "MueLu_RepartitionFactory_fwd.hpp"
#include "MueLu_RebalanceAcFactory_fwd.hpp"
#include "MueLu_RebalanceTransferFactory_fwd.hpp"

#include "MueLu_SmootherFactory_fwd.hpp"
#include "MueLu_TrilinosSmoother_fwd.hpp"
#include "MueLu_Hierarchy_fwd.hpp"

#include "Xpetra_Operator.hpp"
#include "Xpetra_Map_fwd.hpp"
#include "Xpetra_Matrix_fwd.hpp"
#include "Xpetra_MapFactory_fwd.hpp"
#include "Xpetra_MatrixFactory_fwd.hpp"
#include "Xpetra_MultiVectorFactory_fwd.hpp"
#include "Xpetra_VectorFactory_fwd.hpp"
#include "Xpetra_CrsMatrixWrap_fwd.hpp"

namespace MueLu {

/*!
  @brief Preconditioner (wrapped as a Xpetra::Operator) for Maxwell's equations in curl-curl form.

  Let

    M_k(gamma)

  be the mass matrix with a weight coefficient gamma on the k-th
  space in the deRham complex

    k=0     ---> k=1     ----> k=2    ----> k=3

    H(grad) ---> H(curl) ----> H(div) ----> L^2

            D_0          D_1          D_2

            grad         curl         div

  and let

    D_k

  be the discretized derivative from the k-th to the (k+1)-th space.

  For example, in 3D, D_0 = grad, D_1 = curl and D_2 = div.

  We want to solve the system

    S e_k = b

  where

    S = M_k(alpha) + D_k^T * M_{k+1}(beta) * D_k.

  In 3D and for k=1, S corresponds to a discretization of

    alpha * identity + beta * curl curl.

  In 3D and for k=2, S corresponds to

    alpha * identity + beta * grad div.

  We precondition S by constructing a block diagonal preconditioner
  for the equivalent 2x2 block system

    e_k = a_k + D_{k-1}*p_{k-1}

    ( A11                      M_k(alpha) * D_{k-1} ) ( a_k     )   ( b             )
    (                                               } (         ) = (               )
    ( D_{k-1}^T * M_k(alpha)   A22                  ) ( p_{k-1} )   ( D_{k-1}^T * b )

  where

    A11     = S                        + addon11
    A22     = D_{k-1} ^T * S * D_{k-1} + addon22
    addon11 = M_k(1)     * D_{k-1} * M_{k-1}(1/beta)^-1  * D_{k-1}^T * M_k(1)
    addon22 = M_{k-1}(1) * D_{k-2} * M_{k-2}(1/alpha)^-1 * D_{k-2}^T * M_{k-1}(1)

  We note that due to D_k * D_{k-1} = 0 we have that

    A22 = D_{k-1} ^T * M_k(alpha) * D_{k-1} + addon22

  The addon terms mimic the term that would be needed to augment to
  a Hodge Laplacian. In practice it seems that the addon terms are
  rarely strictly required. For k=1 (Maxwell) the addon22 term (in
  H(grad)) is always zero.

  We cannot directly apply AMG to A11 and A22 in case they are not
  expressed in terms of nodal bases. Let Pi_k be the projection from
  the first order vectorial nodal finite element space into the k-th
  space. We will apply AMG to

    Pi_k^T     * A11 * Pi_k
    Pi_{k-1}^T * A22 * Pi_{k-1}

  respectively. For k=1 A22 is already given in a nodal
  discretization and this step is omitted.

  It can be beneficial to directly coarsen the projected diagonal
  blocks and skip any smoothing on them.

  To obtain prolongators for the vectorial nodal problems, we
  construct the auxiliary operators

    A11_nodal = D_{0}^T * M_1(beta)  * D_{0}
    A22_nodal = D_{0}^T * M_1(alpha) * D_{0}

  then constuct typical nodal scalar HGrad prolongators.

  We then replicate them into prolongators for vectorial H(grad)
  problem. Finally, we multiply from the left with the projection
  Pi_k onto the k-th FE space.

  When k=1 this only needs to be done for the A11 block, as the A22
  block is already given wrt a scalar nodal discretization.

  The following input matrices need to be supplied by the user:
  - S
  - D_{k-1}

  If addon11 is used we need:
  - M_{k-1}(1/beta)^-1   - typically, this is the inverse of the lumped M_{k-1}(1/beta), not the actual inverse
  - M_k(1)

  If addon22 is used we need:
  - M_{k-2}(1/alpha)^-1  - typically, this is the inverse of the lumped M_{k-2}(1/alpha), not the actual inverse
  - M_{k-1}(1)

  To construct the special prolongators we need
  - D_{0}
  - M_1(beta)
  - M_1(alpha)
  If these mass matrices are not given, but M_1(1) is, then that matrix can be used instead.
  Alternatively, A11_nodal and A22_nodal can be passed directly.


  We are using the following variable names:

    | variable                         | matrix                | note
    |----------------------------------|-----------------------|----------------------
    | <tt>SM               </tt>       | S                     |
    | <tt>Dk_1             </tt>       | D_{k-1}               | same as D0 for k=1
    | <tt>Dk_2             </tt>       | D_{k-2}               |
    | <tt>D0               </tt>       | D_0                   | same as Dk_1 for k=1
    | <tt>Mk_one           </tt>       | M_k(1)                |
    | <tt>Mk_1_one         </tt>       | M_{k-1}(1)            |
    | <tt>M1_beta          </tt>       | M_1(beta)             |
    | <tt>M1_alpha         </tt>       | M_1(alpha)            |
    | <tt>invMk_1_invBeta  </tt>       | M_{k-1}(1/beta)^{-1}  |
    | <tt>invMk_2_invAlpha </tt>       | M_{k-2}(1/alpha)^{-1} |

  For backwards compatibility the interfaces also allow

    | variable          | matrix           | note
    |-------------------|------------------|---------------------------------
    | <tt>Ms    </tt>   | M_1(beta)        | alias for M1_beta
    | <tt>M1    </tt>   | M_1(1)           | alias for Mk_one when k=1
    | <tt>M0inv </tt>   | M_0(1/beta)      | alias for Mk_1_invBeta when k=1


  Reference:
  P. Bochev, J. Hu, C. Siefert, and R. Tuminaro. "An algebraic multigrid approach based on
  a compatible gauge reformulation of Maxwell's equations." SIAM Journal on Scientific
  Computing, 31(1), 557-583.

  Parameter list options:
  - <tt>refmaxwell: mode</tt> - a <tt>string</tt> specifying the order of solve of the block system.
                                Allowed values are: "additive" (default), "121", "212", "1", "2"
  - <tt>refmaxwell: disable addon</tt> - <tt>bool</tt> specifing whether the addon should be built for stabilization.
                                         Default: "true"
  - <tt>refmaxwell: use as preconditioner</tt> - <tt>bool</tt> specifing whether RefMaxwell is used as a preconditioner or as a solver.
  - <tt>refmaxwell: dump matrices</tt> - <tt>bool</tt> specifing whether the matrices should be dumped.
                                         Default: "false"
  - <tt>refmaxwell: 11list</tt> and <tt>refmaxwell: 22list</tt> - parameter list for the multigrid hierarchies on 11 and 22 blocks
  - <tt>refmaxwell: subsolves on subcommunicators</tt> - <tt>bool</tt> redistribute the two subsolves to disjoint sub-communicators (so that the additive solve can occur in parallel)
                                         Default: "false"

  @ingroup MueLuAdapters
*/
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class RefMaxwell : public VerboseObject, public Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_REFMAXWELL_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;
  typedef typename Teuchos::ScalarTraits<Scalar>::coordinateType coordinateType;
  typedef typename Xpetra::MultiVector<coordinateType, LO, GO, NO> RealValuedMultiVector;

  //! Constructor
  RefMaxwell()
    : HierarchyCoarse11_(Teuchos::null)
    , Hierarchy22_(Teuchos::null)
    , disable_addon_(MasterList::getDefault<bool>("refmaxwell: disable addon"))
    , mode_(MasterList::getDefault<std::string>("refmaxwell: mode")) {
  }

  //! Constructor with Hierarchies
  RefMaxwell(Teuchos::RCP<Hierarchy> HH, Teuchos::RCP<Hierarchy> H22)
    : HierarchyCoarse11_(HH)
    , Hierarchy22_(H22)
    , disable_addon_(MasterList::getDefault<bool>("refmaxwell: disable addon"))
    , mode_(MasterList::getDefault<std::string>("refmaxwell: mode")) {
  }

  RefMaxwell(const Teuchos::RCP<Matrix> &SM_Matrix,
             const Teuchos::RCP<Matrix> &Dk_1,
             const Teuchos::RCP<Matrix> &Dk_2,
             const Teuchos::RCP<Matrix> &D0,
             const Teuchos::RCP<Matrix> &M1_beta,
             const Teuchos::RCP<Matrix> &M1_alpha,
             const Teuchos::RCP<Matrix> &Mk_one,
             const Teuchos::RCP<Matrix> &Mk_1_one,
             const Teuchos::RCP<Matrix> &invMk_1_invBeta,
             const Teuchos::RCP<Matrix> &invMk_2_invAlpha,
             const Teuchos::RCP<MultiVector> &Nullspace11,
             const Teuchos::RCP<MultiVector> &Nullspace22,
             const Teuchos::RCP<RealValuedMultiVector> &NodalCoords,
             Teuchos::ParameterList &List,
             bool ComputePrec = true) {
    int spaceNumber = List.get<int>("refmaxwell: space number", 1);
    initialize(spaceNumber,
               Dk_1, Dk_2, D0,
               M1_beta, M1_alpha,
               Mk_one, Mk_1_one,
               invMk_1_invBeta, invMk_2_invAlpha,
               Nullspace11, Nullspace22, NodalCoords,
               List);
    resetMatrix(SM_Matrix, ComputePrec);
  }

  /** Constructor with Jacobian (with add on)
   *
   * \param[in] SM_Matrix Jacobian
   * \param[in] D0_Matrix Discrete Gradient
   * \param[in] Ms_Matrix Edge mass matrix for the nodal aggregates
   * \param[in] M0inv_Matrix Inverse of lumped nodal mass matrix (add on only)
   * \param[in] M1_Matrix Edge mass matrix for the add on
   * \param[in] Nullspace11 Null space (needed for periodic)
   * \param[in] NodalCoords Nodal coordinates
   * \param[in] List Parameter list
   * \param[in] ComputePrec If true, compute the preconditioner immediately
   */
  RefMaxwell(const Teuchos::RCP<Matrix> &SM_Matrix,
             const Teuchos::RCP<Matrix> &D0_Matrix,
             const Teuchos::RCP<Matrix> &Ms_Matrix,
             const Teuchos::RCP<Matrix> &M0inv_Matrix,
             const Teuchos::RCP<Matrix> &M1_Matrix,
             const Teuchos::RCP<MultiVector> &Nullspace11,
             const Teuchos::RCP<RealValuedMultiVector> &NodalCoords,
             Teuchos::ParameterList &List,
             bool ComputePrec = true) {
    initialize(D0_Matrix, Ms_Matrix, M0inv_Matrix, M1_Matrix, Nullspace11, NodalCoords, List);
    resetMatrix(SM_Matrix, ComputePrec);
  }

  /** Constructor with Jacobian (with add on)
   *
   * \param[in] SM_Matrix Jacobian
   * \param[in] D0_Matrix Discrete Gradient
   * \param[in] M0inv_Matrix Inverse of lumped nodal mass matrix (add on only)
   * \param[in] M1_Matrix Edge mass matrix for the
   * \param[in] Nullspace Null space (needed for periodic)
   * \param[in] Coords Nodal coordinates
   * \param[in] List Parameter list
   * \param[in] ComputePrec If true, compute the preconditioner immediately
   */
  RefMaxwell(const Teuchos::RCP<Matrix> &SM_Matrix,
             const Teuchos::RCP<Matrix> &D0_Matrix,
             const Teuchos::RCP<Matrix> &M0inv_Matrix,
             const Teuchos::RCP<Matrix> &M1_Matrix,
             const Teuchos::RCP<MultiVector> &Nullspace11,
             const Teuchos::RCP<RealValuedMultiVector> &NodalCoords,
             Teuchos::ParameterList &List,
             bool ComputePrec = true) {
    initialize(D0_Matrix, M1_Matrix, M0inv_Matrix, M1_Matrix, Nullspace11, NodalCoords, List);
    resetMatrix(SM_Matrix, ComputePrec);
  }

  /** Constructor without Jacobian (with add on)
   *
   * \param[in] D0_Matrix Discrete Gradient
   * \param[in] M0inv_Matrix Inverse of lumped nodal mass matrix (add on only)
   * \param[in] M1_Matrix Edge mass matrix for the
   * \param[in] Nullspace Null space (needed for periodic)
   * \param[in] Coords Nodal coordinates
   * \param[in] List Parameter list
   */
  RefMaxwell(const Teuchos::RCP<Matrix> &D0_Matrix,
             const Teuchos::RCP<Matrix> &M0inv_Matrix,
             const Teuchos::RCP<Matrix> &M1_Matrix,
             const Teuchos::RCP<MultiVector> &Nullspace11,
             const Teuchos::RCP<RealValuedMultiVector> &NodalCoords,
             Teuchos::ParameterList &List)
    : SM_Matrix_(Teuchos::null) {
    initialize(D0_Matrix, M1_Matrix, M0inv_Matrix, M1_Matrix, Nullspace11, NodalCoords, List);
  }

  /** Constructor with Jacobian (no add on)
   *
   * \param[in] SM_Matrix Jacobian
   * \param[in] D0_Matrix Discrete Gradient
   * \param[in] M1_Matrix Edge mass matrix for the
   * \param[in] Nullspace Null space (needed for periodic)
   * \param[in] Coords Nodal coordinates
   * \param[in] List Parameter list
   * \param[in] ComputePrec If true, compute the preconditioner immediately
   */
  RefMaxwell(const Teuchos::RCP<Matrix> &SM_Matrix,
             const Teuchos::RCP<Matrix> &D0_Matrix,
             const Teuchos::RCP<Matrix> &M1_Matrix,
             const Teuchos::RCP<MultiVector> &Nullspace11,
             const Teuchos::RCP<RealValuedMultiVector> &NodalCoords,
             Teuchos::ParameterList &List,
             bool ComputePrec) {
    initialize(D0_Matrix, M1_Matrix, Teuchos::null, M1_Matrix, Nullspace11, NodalCoords, List);
    resetMatrix(SM_Matrix, ComputePrec);
  }

  /** Constructor without Jacobian (no add on)
   *
   * \param[in] D0_Matrix Discrete Gradient
   * \param[in] M1_Matrix Edge mass matrix for the
   * \param[in] Nullspace Null space (needed for periodic)
   * \param[in] Coords Nodal coordinates
   * \param[in] List Parameter list
   */
  RefMaxwell(const Teuchos::RCP<Matrix> &D0_Matrix,
             const Teuchos::RCP<Matrix> &M1_Matrix,
             const Teuchos::RCP<MultiVector> &Nullspace11,
             const Teuchos::RCP<RealValuedMultiVector> &NodalCoords,
             Teuchos::ParameterList &List)
    : SM_Matrix_(Teuchos::null) {
    initialize(D0_Matrix, M1_Matrix, Teuchos::null, M1_Matrix, Nullspace11, NodalCoords, List);
  }

  /** Constructor with parameter list
   *
   * \param[in] SM_Matrix Jacobian
   * \param[in] List Parameter list
   * \param[in] ComputePrec If true, compute the preconditioner immediately
   */
  RefMaxwell(const Teuchos::RCP<Matrix> &SM_Matrix,
             Teuchos::ParameterList &List,
             bool ComputePrec = true);

  //! Destructor.
  virtual ~RefMaxwell() {}

  //! Returns the Xpetra::Map object associated with the domain of this operator.
  Teuchos::RCP<const Map> getDomainMap() const;

  //! Returns the Xpetra::Map object associated with the range of this operator.
  Teuchos::RCP<const Map> getRangeMap() const;

  //! Returns Jacobian matrix SM
  const Teuchos::RCP<Matrix> &getJacobian() const {
    return SM_Matrix_;
  }

  //! Set parameters
  void setParameters(Teuchos::ParameterList &list);

  //! Setup the preconditioner
  void compute(bool reuse = false);

  //! Reset system matrix
  void resetMatrix(Teuchos::RCP<Matrix> SM_Matrix_new, bool ComputePrec = true);

  //! Returns in Y the result of a Xpetra::Operator applied to a Xpetra::MultiVector X.
  //! \param[in]  X - MultiVector of dimension NumVectors to multiply with matrix.
  //! \param[out] Y - MultiVector of dimension NumVectors containing result.
  void apply(const MultiVector &X, MultiVector &Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS,
             Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
             Scalar beta           = Teuchos::ScalarTraits<Scalar>::zero()) const;

  //! Indicates whether this operator supports applying the adjoint operator.
  bool hasTransposeApply() const;

  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_HIGH) const;

  //! Compute a residual R = B - (*this) * X
  void residual(const MultiVector &X,
                const MultiVector &B,
                MultiVector &R) const {
    using STS = Teuchos::ScalarTraits<Scalar>;
    R.update(STS::one(), B, STS::zero());
    this->apply(X, R, Teuchos::NO_TRANS, -STS::one(), STS::one());
  }

 private:
  Teuchos::RCP<Teuchos::ParameterList> getValidParamterList();

  /** Initialize with matrices except the Jacobian (don't compute the preconditioner)
   *
   * Note: This uses old notation that only makes sense for curl-curl problems.
   *
   * \param[in] D0_Matrix Discrete Gradient
   * \param[in] Ms_Matrix Edge mass matrix for nodal aggregates
   * \param[in] M0inv_Matrix Inverse of lumped nodal mass matrix (add on only)
   * \param[in] M1_Matrix Edge mass matrix for add on
   * \param[in] Nullspace Null space (needed for periodic)
   * \param[in] Coords Nodal coordinates
   * \param[in] List Parameter list
   */
  void initialize(const Teuchos::RCP<Matrix> &D0_Matrix,
                  const Teuchos::RCP<Matrix> &Ms_Matrix,
                  const Teuchos::RCP<Matrix> &M0inv_Matrix,
                  const Teuchos::RCP<Matrix> &M1_Matrix,
                  const Teuchos::RCP<MultiVector> &Nullspace11,
                  const Teuchos::RCP<RealValuedMultiVector> &NodalCoords,
                  Teuchos::ParameterList &List);

  /** Initialize with matrices except the Jacobian (don't compute the preconditioner)
   *
   * \param[in] k number of the space in the deRham sequence of the problem to be solved
   * \param[in] Dk_1 Discrete derivative from (k-1)-th to k-th space
   * \param[in] Dk_2 Discrete derivative from (k-2)-th to (k-1)-th space
   * \param[in] D0 Discrete Gradient
   * \param[in] M1_beta Mass matrix on 1-st space with weight beta for nodal aggregates
   * \param[in] M1_alpha Mass matrix on 1-st space with weight alpha for nodal aggregates
   * \param[in] Mk_one Mass matrix on k-th space with unit weight for addon11
   * \param[in] Mk_1_one Mass matrix on (k-1)-th space with unit weight for addon22
   * \param[in] invMk_1_invBeta Approximate inverse of mass matrix on (k-1)-th space with weight 1/beta (addon11 only)
   * \param[in] invMk_2_invAlpha Approximate inverse of mass matrix on (k-2)-th space with weight 1/alpha (addon22 only)
   * \param[in] Nullspace11 Null space (needed for periodic)
   * \param[in] Nullspace22 Null space (needed for periodic)
   * \param[in] NodalCoords Nodal coordinates
   * \param[in] List Parameter list
   */
  void initialize(const int k,
                  const Teuchos::RCP<Matrix> &Dk_1,
                  const Teuchos::RCP<Matrix> &Dk_2,
                  const Teuchos::RCP<Matrix> &D0,
                  const Teuchos::RCP<Matrix> &M1_beta,
                  const Teuchos::RCP<Matrix> &M1_alpha,
                  const Teuchos::RCP<Matrix> &Mk_one,
                  const Teuchos::RCP<Matrix> &Mk_1_one,
                  const Teuchos::RCP<Matrix> &invMk_1_invBeta,
                  const Teuchos::RCP<Matrix> &invMk_2_invAlpha,
                  const Teuchos::RCP<MultiVector> &Nullspace11,
                  const Teuchos::RCP<MultiVector> &Nullspace22,
                  const Teuchos::RCP<RealValuedMultiVector> &NodalCoords,
                  Teuchos::ParameterList &List);

  //! Determine how large the sub-communicators for the two hierarchies should be
  void determineSubHierarchyCommSizes(bool &doRebalancing, int &rebalanceStriding, int &numProcsCoarseA11, int &numProcsA22);

  //! Compute coarseA11 = P11^{T}*SM*P11 + addon efficiently
  void buildCoarse11Matrix();

  //! rebalance the coarse A11 matrix, as well as P11, CoordsCoarse11 and Addon11
  void rebalanceCoarse11Matrix(const int rebalanceStriding, const int numProcsCoarseA11);

  //! Setup A22 = D0^T SM D0 and rebalance it, as well as D0 and Coords_
  void build22Matrix(const bool reuse, const bool doRebalancing, const int rebalanceStriding, const int numProcsA22);

 public:  // due to Cuda build errors otherwise
  //! Builds a nullspace
  RCP<MultiVector> buildNullspace(const int spaceNumber, const Kokkos::View<bool *, typename Node::device_type> &bcs, const bool applyBCs);

  //! Builds a projection from a vector values space into a vector valued nodal space
  Teuchos::RCP<Matrix> buildProjection(const int spaceNumber, const RCP<MultiVector> &EdgeNullspace) const;

  /** Setup an auxiliary nodal prolongator
   *
   * \param[in]  A_nodal
   * \param[out] P_nodal
   * \param[out] Nullspace_nodal
   */
  void buildNodalProlongator(const Teuchos::RCP<Matrix> &A_nodal,
                             Teuchos::RCP<Matrix> &P_nodal,
                             Teuchos::RCP<MultiVector> &Nullspace_nodal,
                             Teuchos::RCP<RealValuedMultiVector> &Coords_nodal) const;

  /** Setup a vectorial nodal prolongator
   *
   * \param[in]  P_nodal_scalar
   * \param[out] P_nodal_vector
   */
  RCP<Matrix> buildVectorNodalProlongator(const Teuchos::RCP<Matrix> &P_nodal) const;

  /** Setup a special prolongator from vectorial nodal to edge or face discretization
   *
   * \param[in]  spaceNumber the type of target discretization
   * \param[in]  A_nodal_Matrix used for the nodal aggregation
   * \param[in]  EdgeNullspace edge nullspace
   * \param[out] edgeProlongator edge prolongator
   * \param[out] coarseEdgeNullspace coarse edge nullspace
   * \param[out] coarseNodalCoords coarse nodal coordinates
   */
  void buildProlongator(const int spaceNumber,
                        const Teuchos::RCP<Matrix> &A_nodal_Matrix,
                        const RCP<MultiVector> &EdgeNullspace,
                        Teuchos::RCP<Matrix> &edgeProlongator,
                        Teuchos::RCP<MultiVector> &coarseEdgeNullspace,
                        Teuchos::RCP<RealValuedMultiVector> &coarseNodalCoords) const;

  /** Construct an addon matrix
   *
   * \param[in]  spaceNumber the type of target discretization
   */
  RCP<Matrix> buildAddon(const int spaceNumber);

 private:
  //! Setup a subsolve
  void setupSubSolve(Teuchos::RCP<Hierarchy> &hierarchy,
                     Teuchos::RCP<Operator> &thyraPrecOp,
                     const Teuchos::RCP<Matrix> &A,
                     const Teuchos::RCP<MultiVector> &Nullspace,
                     const Teuchos::RCP<RealValuedMultiVector> &Coords,
                     Teuchos::ParameterList &params,
                     std::string &label,
                     const bool reuse,
                     const bool isSingular = false);

  //! Set the fine level smoother
  void setFineLevelSmoother11();

  //! apply additive algorithm for 2x2 solve
  void applyInverseAdditive(const MultiVector &RHS, MultiVector &X) const;

  //! apply solve to 1-1 block only
  void solveH(const MultiVector &RHS, MultiVector &X) const;

  //! apply solve to 2-2 block only
  void solve22(const MultiVector &RHS, MultiVector &X) const;

  //! allocate multivectors for solve
  void allocateMemory(int numVectors) const;

  //! dump out matrix
  void dump(const RCP<Matrix> &A, std::string name) const;

  //! dump out multivector
  void dump(const RCP<MultiVector> &X, std::string name) const;

  //! dump out real-valued multivector
  void dumpCoords(const RCP<RealValuedMultiVector> &X, std::string name) const;

  //! dump out boolean ArrayView
  void dump(const Teuchos::ArrayRCP<bool> &v, std::string name) const;

  //! dump out boolean Kokkos::View
  void dump(const Kokkos::View<bool *, typename Node::device_type> &v, std::string name) const;

  //! get a (synced) timer
  Teuchos::RCP<Teuchos::TimeMonitor> getTimer(std::string name, RCP<const Teuchos::Comm<int> > comm = Teuchos::null) const;

  //! Two hierarchies: one for the coarse (1,1)-block, another for the (2,2)-block
  Teuchos::RCP<Hierarchy> HierarchyCoarse11_, Hierarchy22_;
  Teuchos::RCP<SmootherBase> PreSmoother11_, PostSmoother11_;
  Teuchos::RCP<SmootherPrototype> PreSmootherData11_, PostSmootherData11_;
  RCP<Operator> thyraPrecOpH_, thyraPrecOp22_;
  //! The number of the space in the deRham complex
  int spaceNumber_;
  //! The name of the solver
  std::string solverName_;
  //! The spatial dimension
  size_t dim_;
  //! The system that is getting preconditioned
  Teuchos::RCP<Matrix> SM_Matrix_;
  //! D_{k-1} matrix and its transpose
  Teuchos::RCP<Matrix> Dk_1_, Dk_1_T_;
  //! D_{k-2} matrix
  Teuchos::RCP<Matrix> Dk_2_;
  //! D_0 matrix
  Teuchos::RCP<Matrix> D0_;
  //! inverse of mass matrices on (k-1)-th and (k-2)-th space with weights 1/beta and 1/alpha respectively
  Teuchos::RCP<Matrix> invMk_1_invBeta_, invMk_2_invAlpha_;
  //! mass matrices with unit weight on k-th and (k-1)-th spaces
  Teuchos::RCP<Matrix> Mk_one_, Mk_1_one_;
  //! mass matrices on first space with weights beta and alpha respectively
  Teuchos::RCP<Matrix> M1_beta_, M1_alpha_;
  //! special prolongator for 11 block and its transpose
  Teuchos::RCP<Matrix> P11_, R11_;
  //! special prolongator for 22 block and its transpose
  Teuchos::RCP<Matrix> P22_, R22_;
  //! coarse 11, 22 and coarse 22 blocks
  Teuchos::RCP<Matrix> coarseA11_, A22_, coarseA22_;
  //! the addon for the 11 block
  Teuchos::RCP<Matrix> Addon11_;
  //! the addon for the 22 block
  Teuchos::RCP<Matrix> Addon22_;
  Teuchos::RCP<const Map> DorigDomainMap_;
  Teuchos::RCP<const Import> DorigImporter_;
  //! Vectors for BCs
  Kokkos::View<bool *, typename Node::device_type> BCrows11_, BCcols22_, BCdomain22_;
  int globalNumberBoundaryUnknowns11_, globalNumberBoundaryUnknowns22_;
  //! Nullspace for (1.1) block
  Teuchos::RCP<MultiVector> Nullspace11_, Nullspace22_;
  //! Coordinates
  Teuchos::RCP<RealValuedMultiVector> NodalCoords_, CoordsCoarse11_, Coords22_;
  //! Nullspace for coarse (1,1) problem
  Teuchos::RCP<MultiVector> NullspaceCoarse11_;
  //! Nullspace for coarse (2,2) problem
  Teuchos::RCP<MultiVector> CoarseNullspace22_;
  //! Importer to coarse (1,1) hierarchy
  Teuchos::RCP<const Import> ImporterCoarse11_, Importer22_;
  bool Dk_1_T_R11_colMapsMatch_;
  bool onlyBoundary11_, onlyBoundary22_;
  //! Parameter lists
  Teuchos::ParameterList parameterList_, precList11_, precList22_;
  Teuchos::RCP<Teuchos::ParameterList> coarseA11_AP_reuse_data_, coarseA11_RAP_reuse_data_;
  Teuchos::RCP<Teuchos::ParameterList> A22_AP_reuse_data_, A22_RAP_reuse_data_;
  //! Some options
  bool disable_addon_, disable_addon_22_, dump_matrices_, useKokkos_, use_as_preconditioner_, implicitTranspose_, fuseProlongationAndUpdate_, syncTimers_, enable_reuse_, skipFirst11Level_, skipFirst22Level_;
  bool applyBCsToAnodal_, applyBCsToCoarse11_, applyBCsTo22_;
  int numItersCoarse11_, numIters22_;
  std::string mode_;

  //! Temporary memory
  mutable Teuchos::RCP<MultiVector> P11res_, P11x_, P11resSubComm_, P11xSubComm_, DresIntermediate_, Dres_, DxIntermediate_, Dx_, DresSubComm_, DxSubComm_, residual_, P11resTmp_, DresTmp_, DTR11Tmp_;
};

}  // namespace MueLu

#define MUELU_REFMAXWELL_SHORT
#endif  // MUELU_REFMAXWELL_DECL_HPP
