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
#include "MueLu_Utilities_fwd.hpp"
#include "MueLu_TentativePFactory_fwd.hpp"
#include "MueLu_SaPFactory_fwd.hpp"
#include "MueLu_UncoupledAggregationFactory_fwd.hpp"
#include "MueLu_SmootherFactory_fwd.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Hierarchy.hpp"

#include "Xpetra_Map_fwd.hpp"
#include "Xpetra_Matrix_fwd.hpp"
#include "Xpetra_MatrixFactory_fwd.hpp"
#include "Xpetra_MultiVectorFactory_fwd.hpp"
#include "Xpetra_CrsMatrixWrap_fwd.hpp"
#include "Xpetra_BlockedCrsMatrix_fwd.hpp"
#include "Xpetra_ExportFactory_fwd.hpp"

#ifdef HAVE_MUELU_IFPACK2
#include "Ifpack2_Preconditioner.hpp"
#if defined(HAVE_TPETRA_INST_INT_INT)
#include "Ifpack2_Hiptmair.hpp"
#endif
#endif

namespace MueLu {

  /*!
    @brief Preconditioner (wrapped as a Tpetra::Operator) for Maxwell's equations in curl-curl form.

    This uses a 2x2 block reformulation.

    Reference:
    P. Bochev, J. Hu, C. Siefert, and R. Tuminaro. "An algebraic multigrid approach based on
    a compatible gauge reformulation of Maxwell's equations." SIAM Journal on Scientific
    Computing, 31(1), 557-583.

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
     * \param[in] coords Nodal coordinates
     * \param[in] precList Parameter list
     * \param[in] ComputePrec If true, compute the preconditioner immediately
     */
    RefMaxwell(const Teuchos::RCP<Matrix> & SM_Matrix,
               const Teuchos::RCP<Matrix> & D0_Matrix,
               const Teuchos::RCP<Matrix> & M0inv_Matrix,
               const Teuchos::RCP<Matrix> & M1_Matrix,
               const Teuchos::RCP<MultiVector> & Nullspace,
               const Teuchos::RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > & Coords,
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
     * \param[in] coords Nodal coordinates
     * \param[in] precList Parameter list
     */
    RefMaxwell(const Teuchos::RCP<Matrix> & D0_Matrix,
               const Teuchos::RCP<Matrix> & M0inv_Matrix,
               const Teuchos::RCP<Matrix> & M1_Matrix,
               const Teuchos::RCP<MultiVector> & Nullspace,
               const Teuchos::RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > & Coords,
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
     * \param[in] coords Nodal coordinates
     * \param[in] precList Parameter list
     * \param[in] ComputePrec If true, compute the preconditioner immediately
     */
    RefMaxwell(const Teuchos::RCP<Matrix> & SM_Matrix,
               const Teuchos::RCP<Matrix> & D0_Matrix,
               const Teuchos::RCP<Matrix> & M1_Matrix,
               const Teuchos::RCP<MultiVector>  & Nullspace,
               const Teuchos::RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> >  & Coords,
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
     * \param[in] coords Nodal coordinates
     * \param[in] precList Parameter list
     */
    RefMaxwell(const Teuchos::RCP<Matrix> & D0_Matrix,
               const Teuchos::RCP<Matrix> & M1_Matrix,
               const Teuchos::RCP<MultiVector>  & Nullspace,
               const Teuchos::RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> >  & Coords,
               Teuchos::ParameterList& List) : SM_Matrix_(Teuchos::null)
    {
      initialize(D0_Matrix,Teuchos::null,M1_Matrix,Nullspace,Coords,List);
    }

    /** Constructor with parameter list
     *
     * \param[in] SM_Matrix Jacobian
     * \param[in] precList Parameter list
     * \param[in] ComputePrec If true, compute the preconditioner immediately
     */
    RefMaxwell(const Teuchos::RCP<Matrix> & SM_Matrix,
               Teuchos::ParameterList& List,
               bool ComputePrec = true)
    {

      RCP<MultiVector> Nullspace = List.get<RCP<MultiVector> >("Nullspace", Teuchos::null);
      RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > Coords = List.get<RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > >("Coordinates", Teuchos::null);
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

    //! apply additive algorithm for 2x2 solve
    void applyInverseAdditive(const MultiVector& RHS, MultiVector& X) const;

    //! apply 1-2-1 algorithm for 2x2 solve
    void applyInverse121(const MultiVector& RHS, MultiVector& X) const;

    //! apply 2-1-2 algorithm for 2x2 solve
    void applyInverse212(const MultiVector& RHS, MultiVector& X) const;

    //! apply solve to 1-1 block only
    void applyInverse11only(const MultiVector& RHS, MultiVector& X) const;

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

  private:

    /** Initialize with matrices except the Jacobian (don't compute the preconditioner)
     *
     * \param[in] D0_Matrix Discrete Gradient
     * \param[in] M0inv_Matrix Inverse of lumped nodal mass matrix (add on only)
     * \param[in] M1_Matrix Edge mass matrix
     * \param[in] Nullspace Null space (needed for periodic)
     * \param[in] coords Nodal coordinates
     * \param[in] precList Parameter list
     */
    void initialize(const Teuchos::RCP<Matrix> & D0_Matrix,
                    const Teuchos::RCP<Matrix> & M0inv_Matrix,
                    const Teuchos::RCP<Matrix> & M1_Matrix,
                    const Teuchos::RCP<MultiVector> & Nullspace,
                    const Teuchos::RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > & Coords,
                    Teuchos::ParameterList& List);

    //! Two hierarchies: one for the coarse (1,1)-block, another for the (2,2)-block
    Teuchos::RCP<Hierarchy> HierarchyH_, Hierarchy22_;
    Teuchos::RCP<SmootherBase> Smoother_;
#if defined(HAVE_MUELU_IFPACK2) && defined(HAVE_TPETRA_INST_INT_INT)
    Teuchos::RCP<Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> > hiptmairPreSmoother_, hiptmairPostSmoother_;
#endif
    bool useHiptmairSmoothing_;
    //! Top Level
    Teuchos::RCP<Level> TopLevel_;
    //! Various matrices
    Teuchos::RCP<Matrix> SM_Matrix_, D0_Matrix_, M0inv_Matrix_, M1_Matrix_, Ms_Matrix_;
    Teuchos::RCP<Matrix> A_nodal_Matrix_, P11_, AH_, A22_;
    //! Vectors for BCs
    Teuchos::ArrayRCP<const bool> BCrows_;
    Teuchos::ArrayRCP<const bool> BCcols_;
    //! Nullspace
    Teuchos::RCP<MultiVector> Nullspace_;
    Teuchos::RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > Coords_;
    //! Parameter lists
    Teuchos::ParameterList parameterList_, precList11_, precList22_, smootherList_;
    //! Some options
    bool disable_addon_, dump_matrices_;
    std::string mode_;
    //! Temporary memory
    mutable Teuchos::RCP<MultiVector> P11res_, P11x_, D0res_, D0x_, residual_;
  };

} // namespace

#define MUELU_REFMAXWELL_SHORT
#endif // MUELU_REFMAXWELL_DECL_HPP
