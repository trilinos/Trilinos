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

#include "MueLu.hpp"
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"

#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_IFPACK2)

#include "Tpetra_Operator.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector_decl.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_MatrixFactory.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_BlockedCrsMatrix.hpp"
#include "Xpetra_TpetraMultiVector.hpp"
#include "XpetraExt_MatrixMatrix.hpp"
#include "Xpetra_ExportFactory.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Factory_decl.hpp"
#include "Ifpack2_Factory_def.hpp"

/*
  
  @class RefMaxwell
  
  Preconditioner for Maxwell's equations in curl-curl form using a 2x2 block reformulation,
  wrapped as a Tpetra::Operator.
  
  Reference:
  P. Bochev, J. Hu, C. Siefert, and R. Tuminaro. "An algebraic multigrid approach based on
  a compatible gauge reformulation of Maxwell's equations." SIAM Journal on Scientific
  Computing, 31(1), 557-583.
  
*/

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal,
            class Node = KokkosClassic::DefaultNode::DefaultNodeType,
            class LocalMatOps = typename KokkosClassic::DefaultKernels<Scalar, LocalOrdinal, Node>::SparseOps >
  class RefMaxwell : public Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> {

#undef MUELU_REFMAXWELL_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>                                        TMap;
    typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>                           TCRS;
    typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>                         TMV;
    typedef Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node>                                        XMap;
    typedef Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>                         XMV;
    typedef Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>                   XTMV;
    typedef Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>                           XCRS;
    typedef Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>                     XTCRS;
    typedef Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>                              XMat;
    typedef Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>                       XCrsWrap;
    
    //! Constructor
    RefMaxwell() :
      Hierarchy11_(Teuchos::null),
      Hierarchy22_(Teuchos::null),
      disable_addon_(true),
      MaxCoarseSize_(1000),
      MaxLevels_(5),
      Cycles_(1),
      precType11_("CHEBYSHEV"),
      precType22_("CHEBYSHEV"),
      mode_("additive")
    {
    }
    
    //! Constructor with Hierarchies
    RefMaxwell(Teuchos::RCP<Hierarchy> H11, Teuchos::RCP<Hierarchy> H22) :
      Hierarchy11_(H11),
      Hierarchy22_(H22),
      disable_addon_(false),
      MaxCoarseSize_(1000),
      MaxLevels_(5),
      Cycles_(1),
      precType11_("CHEBYSHEV"),
      precType22_("CHEBYSHEV"),
      mode_("additive")
    {
    }

    //! Constructor with matrices
    RefMaxwell(Teuchos::RCP<TCRS> SM_Matrix,
	       Teuchos::RCP<TCRS> D0_Matrix,
	       Teuchos::RCP<TCRS> M0inv_Matrix,
	       Teuchos::RCP<TCRS> M1_Matrix,
	       Teuchos::RCP<TMV>  Nullspace,
	       Teuchos::RCP<TMV>  Coords,
	       Teuchos::ParameterList& List,
	       bool ComputePrec = true) :
      Hierarchy11_(Teuchos::null),
      Hierarchy22_(Teuchos::null),
      parameterList_(List),
      disable_addon_(false),
      MaxCoarseSize_(1000),
      MaxLevels_(5),
      Cycles_(1),
      precType11_("CHEBYSHEV"),
      precType22_("CHEBYSHEV"),
      mode_("additive")
    {
      // set parameters
      setParameters(List);
      // convert Tpetra matrices to Xpetra
      Teuchos::RCP<XCRS> SM_tmp = Teuchos::rcp( new XTCRS(SM_Matrix) );
      SM_Matrix_ = Teuchos::rcp( new XCrsWrap(SM_tmp) );
      Teuchos::RCP<XCRS> D0_tmp = Teuchos::rcp( new XTCRS(D0_Matrix) );
      D0_Matrix_ = Teuchos::rcp( new XCrsWrap(D0_tmp) );
      Teuchos::RCP<XCRS> M0inv_tmp = Teuchos::rcp( new XTCRS(M0inv_Matrix) );
      M0inv_Matrix_ = Teuchos::rcp( new XCrsWrap(M0inv_tmp) );
      Teuchos::RCP<XCRS> M1_tmp = Teuchos::rcp( new XTCRS(M1_Matrix) );
      M1_Matrix_ = Teuchos::rcp( new XCrsWrap(M1_tmp) );
      // convert Tpetra MultiVector to Xpetra
      if(Coords != Teuchos::null)
	Coords_ = Xpetra::toXpetra(Coords);
      if(Nullspace != Teuchos::null)
	Nullspace_ = Xpetra::toXpetra(Nullspace);
      // compute preconditioner
      compute();
    }

    //! Destructor.
    virtual ~RefMaxwell() {
      // clean up
      Hierarchy11_=Teuchos::null;
      Hierarchy22_=Teuchos::null;
      SM_Matrix_=Teuchos::null;
      D0_Matrix_=Teuchos::null;
      M0inv_Matrix_=Teuchos::null;
      M1_Matrix_=Teuchos::null;
      Ms_Matrix_=Teuchos::null;
      Nullspace_=Teuchos::null;
      Coords_=Teuchos::null;
      TMT_Matrix_=Teuchos::null;
      TMT_Agg_Matrix_=Teuchos::null;
    }

    //! Returns the Tpetra::Map object associated with the domain of this operator.
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getDomainMap() const;

    //! Returns the Tpetra::Map object associated with the range of this operator.
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getRangeMap() const;

    //! Set parameters
    void setParameters(Teuchos::ParameterList& list);
    
    //! Setup the preconditioner
    void compute();

    //! find rows associated with Dirichlet BCs
    void findDirichletRows(Teuchos::RCP<XMat> A,
			   std::vector<LocalOrdinal>& dirichletRows);

    //! find cols associated with Dirichlet BCs
    void findDirichletCols(Teuchos::RCP<XMat> A,
			   std::vector<LocalOrdinal>& dirichletRows,
			   std::vector<LocalOrdinal>& dirichletCols);

    //! apply BCs to rows
    void Apply_BCsToMatrixRows(Teuchos::RCP<XMat>& A, std::vector<LocalOrdinal>& dirichletRows);

    //! apply BCs to cols
    void Apply_BCsToMatrixCols(Teuchos::RCP<XMat>& A, std::vector<LocalOrdinal>& dirichletCols);

    //! add 1's to the diagonal for zeroed out rows
    void Remove_Zeroed_Rows(Teuchos::RCP<XMat>& A, double tol=1.0e-14);

    //! Setup the prolongator for the (1,1)-block
    void buildProlongator();

    //! Compute P11^{T}*A*P11 efficiently
    void formCoarseMatrix();
    
    //! Reset system matrix
    void resetMatrix(Teuchos::RCP<TCRS> SM_Matrix_new);

    //! apply Hiptmair smoothing
    void applyHiptmairSmoother(const XTMV& RHS, XTMV& X) const;

    //! apply additive algorithm for 2x2 solve
    void applyInverseAdditive(const XTMV& RHS, XTMV& X) const;

    //! apply 1-2-1 algorithm for 2x2 solve
    void applyInverse121(const XTMV& RHS, XTMV& X) const;

    //! apply 2-1-2 algorithm for 2x2 solve
    void applyInverse212(const XTMV& RHS, XTMV& X) const;

    //! Returns in Y the result of a Tpetra::Operator applied to a Tpetra::MultiVector X.
    //! \param[in]  X - Tpetra::MultiVector of dimension NumVectors to multiply with matrix.
    //! \param[out] Y - Tpetra::MultiVector of dimension NumVectors containing result.
    void apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
	       Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
	       Teuchos::ETransp mode = Teuchos::NO_TRANS,
	       Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
	       Scalar beta  = Teuchos::ScalarTraits<Scalar>::one()) const;
    
    //! Indicates whether this operator supports applying the adjoint operator.
    bool hasTransposeApply() const;

    template <class NewNode, class NewLocalMatOps>
    Teuchos::RCP< RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, NewNode, NewLocalMatOps> >
    clone(const RCP<NewNode>& new_node) const {
      return Teuchos::rcp(new RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, NewNode, NewLocalMatOps>
			  (Hierarchy11_->template clone<NewNode,NewLocalMatOps>(new_node),
			   Hierarchy22_->template clone<NewNode,NewLocalMatOps>(new_node)));
    }

  private:
    
    //! Two hierarchies: one for the (1,1)-block, another for the (2,2)-block
    Teuchos::RCP<Hierarchy> Hierarchy11_, Hierarchy22_;
    //! Various matrices
    Teuchos::RCP<XMat> SM_Matrix_, D0_Matrix_, M0inv_Matrix_, M1_Matrix_, Ms_Matrix_;
    Teuchos::RCP<XMat> TMT_Matrix_, TMT_Agg_Matrix_, P11_, A11_, A22_;
    //! Vectors for BCs
    std::vector<LocalOrdinal> BCrows_, BCcols_;
    //! Nullspace
    Teuchos::RCP<XMV>  Nullspace_, Coords_;
    //! Parameter lists
    Teuchos::ParameterList parameterList_, precList11_, precList22_;
    //! Ifpack preconditioners for Hiptmair smoothing
    Teuchos::RCP< Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> > nodePrec_, edgePrec_;
    //! Some options
    bool disable_addon_;
    int MaxCoarseSize_, MaxLevels_, Cycles_;
    std::string precType11_, precType22_, mode_;

  };

} // namespace

#endif //ifdef HAVE_MUELU_TPETRA

#define MUELU_REFMAXWELL_SHORT
#endif // MUELU_REFMAXWELL_DECL_HPP
