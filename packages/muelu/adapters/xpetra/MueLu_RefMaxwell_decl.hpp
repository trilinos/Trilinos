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
#include "MueLu_ThresholdAFilterFactory_fwd.hpp"

#include "MueLu_CoalesceDropFactory_fwd.hpp"
#include "MueLu_CoarseMapFactory_fwd.hpp"
#include "MueLu_CoordinatesTransferFactory_fwd.hpp"
#include "MueLu_TentativePFactory_fwd.hpp"
#include "MueLu_UncoupledAggregationFactory_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"

#ifdef HAVE_MUELU_KOKKOS_REFACTOR
#include "MueLu_CoalesceDropFactory_kokkos_fwd.hpp"
#include "MueLu_CoarseMapFactory_kokkos_fwd.hpp"
#include "MueLu_CoordinatesTransferFactory_kokkos_fwd.hpp"
#include "MueLu_TentativePFactory_kokkos_fwd.hpp"
#include "MueLu_Utilities_kokkos_fwd.hpp"
#include "MueLu_UncoupledAggregationFactory_kokkos_fwd.hpp"
#endif

#include "MueLu_ZoltanInterface_fwd.hpp"
#include "MueLu_Zoltan2Interface_fwd.hpp"
#include "MueLu_RepartitionHeuristicFactory_fwd.hpp"
#include "MueLu_RepartitionFactory_fwd.hpp"
#include "MueLu_RebalanceAcFactory_fwd.hpp"
#include "MueLu_RebalanceTransferFactory_fwd.hpp"

#include "MueLu_SmootherFactory_fwd.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Hierarchy.hpp"

#include "Xpetra_Map_fwd.hpp"
#include "Xpetra_Matrix_fwd.hpp"
#include "Xpetra_MatrixFactory_fwd.hpp"
#include "Xpetra_MultiVectorFactory_fwd.hpp"
#include "Xpetra_VectorFactory_fwd.hpp"
#include "Xpetra_CrsMatrixWrap_fwd.hpp"

#if defined(HAVE_MUELU_IFPACK2) && (!defined(HAVE_MUELU_EPETRA) || defined(HAVE_MUELU_INST_DOUBLE_INT_INT))
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Hiptmair.hpp"
#endif

namespace MueLu {

  /*!
    @brief Preconditioner (wrapped as a Xpetra::Operator) for Maxwell's equations in curl-curl form.

    This uses a 2x2 block reformulation.

    Reference:
    P. Bochev, J. Hu, C. Siefert, and R. Tuminaro. "An algebraic multigrid approach based on
    a compatible gauge reformulation of Maxwell's equations." SIAM Journal on Scientific
    Computing, 31(1), 557-583.

    Parameter list options:
    - <tt>refmaxwell: mode</tt> - a <tt>string</tt> specifying the order of solve of the block system.
                                  Allowed values are: "additive" (default), "121", "212"
    - <tt>refmaxwell: disable addon</tt> - <tt>bool</tt> specifing whether the addon should be built for stabilization.
                                           Default: "true"
    - <tt>refmaxwell: use as preconditioner</tt> - <tt>bool</tt> specifing whether RefMaxwell is used as a preconditioner or as a solver.
    - <tt>refmaxwell: dump matrices</tt> - <tt>bool</tt> specifing whether the matrices should be dumped.
                                           Default: "false"
    - <tt>refmaxwell: prolongator compute algorithm</tt> - a <tt>string</tt> specifying the algorithm to build the prolongator.
                                                           Allowed values are: "mat-mat" and "gustavson"
    - <tt>refmaxwell: 11list</tt> and <tt>refmaxwell: 22list</tt> - parameter list for the multigrid hierarchies on 11 and 22 blocks
    - <tt>refmaxwell: subsolves on subcommunicators</tt> - <tt>bool</tt> redistribute the two subsolves to disjoint sub-communicators (so that the additive solve can occur in parallel)
                                           Default: "false"

    @ingroup MueLuAdapters
  */
  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class RefMaxwell : public VerboseObject, public Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> {

#undef MUELU_REFMAXWELL_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;
    typedef typename Xpetra::MultiVector<magnitudeType,LO,GO,NO> RealValuedMultiVector;

    //! Constructor
    RefMaxwell() :
      HierarchyH_(Teuchos::null),
      Hierarchy22_(Teuchos::null),
      disable_addon_(true),
      mode_("additive")
    {
    }

    //! Constructor with Hierarchies
    RefMaxwell(Teuchos::RCP<Hierarchy> HH, Teuchos::RCP<Hierarchy> H22) :
      HierarchyH_(HH),
      Hierarchy22_(H22),
      disable_addon_(false),
      mode_("additive")
    {
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
    RefMaxwell(const Teuchos::RCP<Matrix> & SM_Matrix,
               const Teuchos::RCP<Matrix> & D0_Matrix,
               const Teuchos::RCP<Matrix> & M0inv_Matrix,
               const Teuchos::RCP<Matrix> & M1_Matrix,
               const Teuchos::RCP<MultiVector> & Nullspace,
               const Teuchos::RCP<RealValuedMultiVector> & Coords,
               Teuchos::ParameterList& List,
               bool ComputePrec = true)
    {
      initialize(D0_Matrix,M0inv_Matrix,M1_Matrix,Nullspace,Coords,List);

      resetMatrix(SM_Matrix);

      // compute preconditioner (optionally)
      if(ComputePrec)
        compute();
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
    RefMaxwell(const Teuchos::RCP<Matrix> & D0_Matrix,
               const Teuchos::RCP<Matrix> & M0inv_Matrix,
               const Teuchos::RCP<Matrix> & M1_Matrix,
               const Teuchos::RCP<MultiVector> & Nullspace,
               const Teuchos::RCP<RealValuedMultiVector> & Coords,
               Teuchos::ParameterList& List) : SM_Matrix_(Teuchos::null)
    {
      initialize(D0_Matrix,M0inv_Matrix,M1_Matrix,Nullspace,Coords,List);
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
    RefMaxwell(const Teuchos::RCP<Matrix> & SM_Matrix,
               const Teuchos::RCP<Matrix> & D0_Matrix,
               const Teuchos::RCP<Matrix> & M1_Matrix,
               const Teuchos::RCP<MultiVector>  & Nullspace,
               const Teuchos::RCP<RealValuedMultiVector>  & Coords,
               Teuchos::ParameterList& List,
               bool ComputePrec = true)
    {
      initialize(D0_Matrix,Teuchos::null,M1_Matrix,Nullspace,Coords,List);

      resetMatrix(SM_Matrix);

      // compute preconditioner (optionally)
      if(ComputePrec)
        compute();
    }

    /** Constructor without Jacobian (no add on)
     *
     * \param[in] D0_Matrix Discrete Gradient
     * \param[in] M1_Matrix Edge mass matrix for the
     * \param[in] Nullspace Null space (needed for periodic)
     * \param[in] Coords Nodal coordinates
     * \param[in] List Parameter list
     */
    RefMaxwell(const Teuchos::RCP<Matrix> & D0_Matrix,
               const Teuchos::RCP<Matrix> & M1_Matrix,
               const Teuchos::RCP<MultiVector>  & Nullspace,
               const Teuchos::RCP<RealValuedMultiVector>  & Coords,
               Teuchos::ParameterList& List) : SM_Matrix_(Teuchos::null)
    {
      initialize(D0_Matrix,Teuchos::null,M1_Matrix,Nullspace,Coords,List);
    }

    /** Constructor with parameter list
     *
     * \param[in] SM_Matrix Jacobian
     * \param[in] List Parameter list
     * \param[in] ComputePrec If true, compute the preconditioner immediately
     */
    RefMaxwell(const Teuchos::RCP<Matrix> & SM_Matrix,
               Teuchos::ParameterList& List,
               bool ComputePrec = true)
    {

      RCP<MultiVector> Nullspace = List.get<RCP<MultiVector> >("Nullspace", Teuchos::null);
      RCP<RealValuedMultiVector> Coords = List.get<RCP<RealValuedMultiVector> >("Coordinates", Teuchos::null);
      RCP<Matrix> D0_Matrix = List.get<RCP<Matrix> >("D0");
      RCP<Matrix> M1_Matrix = List.get<RCP<Matrix> >("M1");
      RCP<Matrix> M0inv_Matrix = List.get<RCP<Matrix> >("M0inv", Teuchos::null);

      initialize(D0_Matrix,M0inv_Matrix,M1_Matrix,Nullspace,Coords,List);

      if (SM_Matrix != Teuchos::null)
        resetMatrix(SM_Matrix);

      // compute preconditioner (optionally)
      if(ComputePrec)
        compute();
    }

    //! Destructor.
    virtual ~RefMaxwell() {}

    //! Returns the Xpetra::Map object associated with the domain of this operator.
    Teuchos::RCP<const Map> getDomainMap() const;

    //! Returns the Xpetra::Map object associated with the range of this operator.
    Teuchos::RCP<const Map> getRangeMap() const;

    //! Returns Jacobian matrix SM
    const Teuchos::RCP<Matrix> & getJacobian() const {
      return SM_Matrix_;
    }

    //! Set parameters
    void setParameters(Teuchos::ParameterList& list);

    //! Setup the preconditioner
    void compute();

    //! Setup the prolongator for the (1,1)-block
    void buildProlongator();

    //! Compute P11^{T}*A*P11 efficiently
    void formCoarseMatrix();

    //! Reset system matrix
    void resetMatrix(Teuchos::RCP<Matrix> SM_Matrix_new);

    //! Returns in Y the result of a Xpetra::Operator applied to a Xpetra::MultiVector X.
    //! \param[in]  X - MultiVector of dimension NumVectors to multiply with matrix.
    //! \param[out] Y - MultiVector of dimension NumVectors containing result.
    void apply (const MultiVector& X, MultiVector& Y,
                Teuchos::ETransp mode = Teuchos::NO_TRANS,
                Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
                Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;

    //! Indicates whether this operator supports applying the adjoint operator.
    bool hasTransposeApply() const;

    template <class NewNode>
    Teuchos::RCP< RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, NewNode> >
    clone (const RCP<NewNode>& new_node) const {
      return Teuchos::rcp (new RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, NewNode>
                           (HierarchyH_->template clone<NewNode> (new_node),
                            Hierarchy22_->template clone<NewNode> (new_node)));
    }

    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_HIGH) const;

  private:

    /** Initialize with matrices except the Jacobian (don't compute the preconditioner)
     *
     * \param[in] D0_Matrix Discrete Gradient
     * \param[in] M0inv_Matrix Inverse of lumped nodal mass matrix (add on only)
     * \param[in] M1_Matrix Edge mass matrix
     * \param[in] Nullspace Null space (needed for periodic)
     * \param[in] Coords Nodal coordinates
     * \param[in] List Parameter list
     */
    void initialize(const Teuchos::RCP<Matrix> & D0_Matrix,
                    const Teuchos::RCP<Matrix> & M0inv_Matrix,
                    const Teuchos::RCP<Matrix> & M1_Matrix,
                    const Teuchos::RCP<MultiVector> & Nullspace,
                    const Teuchos::RCP<RealValuedMultiVector> & Coords,
                    Teuchos::ParameterList& List);

    //! solve coarse (1,1) block
    void solveH() const;

    //! solve (2,2) block
    void solve22() const;

    //! apply additive algorithm for 2x2 solve
    void applyInverseAdditive(const MultiVector& RHS, MultiVector& X) const;

    //! apply 1-2-1 algorithm for 2x2 solve
    void applyInverse121(const MultiVector& RHS, MultiVector& X) const;

    //! apply 2-1-2 algorithm for 2x2 solve
    void applyInverse212(const MultiVector& RHS, MultiVector& X) const;

    //! apply solve to 1-1 block only
    void applyInverse11only(const MultiVector& RHS, MultiVector& X) const;

    //! Two hierarchies: one for the coarse (1,1)-block, another for the (2,2)-block
    Teuchos::RCP<Hierarchy> HierarchyH_, Hierarchy22_;
    Teuchos::RCP<SmootherBase> PreSmoother_, PostSmoother_;
#if defined(HAVE_MUELU_IFPACK2) && (!defined(HAVE_MUELU_EPETRA) || defined(HAVE_MUELU_INST_DOUBLE_INT_INT))
    Teuchos::RCP<Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> > hiptmairPreSmoother_, hiptmairPostSmoother_;
#endif
    bool useHiptmairSmoothing_;
    //! Various matrices
    Teuchos::RCP<Matrix> SM_Matrix_, D0_Matrix_, D0_T_Matrix_, M0inv_Matrix_, M1_Matrix_, Ms_Matrix_;
    Teuchos::RCP<Matrix> A_nodal_Matrix_, P11_, R11_, AH_, A22_;
    //! Vectors for BCs
#ifdef HAVE_MUELU_KOKKOS_REFACTOR
    Kokkos::View<const bool*, typename Node::device_type> BCrowsKokkos_;
    Kokkos::View<const bool*, typename Node::device_type> BCcolsKokkos_;
#endif
    Teuchos::ArrayRCP<const bool> BCrows_;
    Teuchos::ArrayRCP<const bool> BCcols_;
    //! Nullspace
    Teuchos::RCP<MultiVector> Nullspace_;
    //! Coordinates
    Teuchos::RCP<RealValuedMultiVector> Coords_, CoordsH_;
    //! Importer to coarse (1,1) hierarchy
    Teuchos::RCP<const Import> ImporterH_, Importer22_;
    //! Parameter lists
    Teuchos::ParameterList parameterList_, precList11_, precList22_, smootherList_;
    //! Some options
    bool disable_addon_, dump_matrices_,useKokkos_,use_as_preconditioner_;
    std::string mode_;
    //! Temporary memory
    mutable Teuchos::RCP<MultiVector> P11res_, P11x_, D0res_, D0x_, residual_, P11resTmp_, P11xTmp_, D0resTmp_, D0xTmp_;
  };

} // namespace

#define MUELU_REFMAXWELL_SHORT
#endif // MUELU_REFMAXWELL_DECL_HPP
