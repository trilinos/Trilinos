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

// TODO move this file to xpetra subfolder

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
  class RefMaxwell : public Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> {

#undef MUELU_REFMAXWELL_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! Constructor
    RefMaxwell() :
      Hierarchy11_(Teuchos::null),
      Hierarchy22_(Teuchos::null),
      disable_addon_(true),
      mode_("additive")
    {
    }

    //! Constructor with Hierarchies
    RefMaxwell(Teuchos::RCP<Hierarchy> H11, Teuchos::RCP<Hierarchy> H22) :
      Hierarchy11_(H11),
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
               const Teuchos::RCP<MultiVector> & Coords,
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
               const Teuchos::RCP<MultiVector> & Coords,
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
               const Teuchos::RCP<MultiVector>  & Coords,
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
               const Teuchos::RCP<MultiVector>  & Coords,
               Teuchos::ParameterList& List) : SM_Matrix_(Teuchos::null)
    {
      initialize(D0_Matrix,Teuchos::null,M1_Matrix,Nullspace,Coords,List);
    }

    //! Destructor.
    virtual ~RefMaxwell() {}

    //! Returns the Xpetra::Map object associated with the domain of this operator.
    Teuchos::RCP<const Map> getDomainMap() const;

    //! Returns the Xpetra::Map object associated with the range of this operator.
    Teuchos::RCP<const Map> getRangeMap() const;

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
                           (Hierarchy11_->template clone<NewNode> (new_node),
                            Hierarchy22_->template clone<NewNode> (new_node)));
    }

  private:

    void findDirichletRows(Teuchos::RCP<Matrix> A,
                                  std::vector<LocalOrdinal>& dirichletRows) {
      dirichletRows.resize(0);
      for(size_t i=0; i<A->getNodeNumRows(); i++) {
        Teuchos::ArrayView<const LocalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> values;
        A->getLocalRowView(i,indices,values);
        int nnz=0;
        for (int j=0; j<indices.size(); j++) {
          // FIXME (mfh 12 Sep 2015) I just replaced abs with the
          // appropriate ScalarTraits call.  However, this is NOT
          // correct for arbitrary scalar types!!!  I'm guessing you
          // should use the equivalent of LAPACK's SFMIN or machine
          // epsilon here.
          if (Teuchos::ScalarTraits<Scalar>::magnitude(values[j]) > 1.0e-16) {
            nnz++;
          }
        }
        if (nnz == 1 || nnz == 2) {
          dirichletRows.push_back(i);
        }
      }
    }

    void findDirichletCols(Teuchos::RCP<Matrix> A,
                                  std::vector<LocalOrdinal>& dirichletRows,
                                  std::vector<LocalOrdinal>& dirichletCols) {
      Teuchos::RCP<const Map> domMap = A->getDomainMap();
      Teuchos::RCP<const Map> colMap = A->getColMap();
      Teuchos::RCP<Export> exporter = ExportFactory::Build(colMap,domMap);
      Teuchos::RCP<MultiVector> myColsToZero = MultiVectorFactory::Build(colMap,1);
      Teuchos::RCP<MultiVector> globalColsToZero = MultiVectorFactory::Build(domMap,1);
      myColsToZero->putScalar((Scalar)0.0);
      globalColsToZero->putScalar((Scalar)0.0);
      for(size_t i=0; i<dirichletRows.size(); i++) {
        Teuchos::ArrayView<const LocalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> values;
        A->getLocalRowView(dirichletRows[i],indices,values);
        for(int j=0; j<indices.size(); j++)
          myColsToZero->replaceLocalValue(indices[j],0,(Scalar)1.0);
      }
      globalColsToZero->doExport(*myColsToZero,*exporter,Xpetra::ADD);
      myColsToZero->doImport(*globalColsToZero,*exporter,Xpetra::INSERT);
      Teuchos::ArrayRCP<const Scalar> myCols = myColsToZero->getData(0);
      dirichletCols.resize(colMap->getNodeNumElements());
      for(size_t i=0; i<colMap->getNodeNumElements(); i++) {
        if(Teuchos::ScalarTraits<Scalar>::magnitude(myCols[i])>0.0)
          dirichletCols[i]=1;
        else
          dirichletCols[i]=0;
      }
    }

    void Apply_BCsToMatrixRows(Teuchos::RCP<Matrix>& A,
                                      std::vector<LocalOrdinal>& dirichletRows) {
      for(size_t i=0; i<dirichletRows.size(); i++) {
        Teuchos::ArrayView<const LocalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> values;
        A->getLocalRowView(dirichletRows[i],indices,values);
        std::vector<Scalar> vec;
        vec.resize(indices.size());
        Teuchos::ArrayView<Scalar> zerovalues(vec);
        for(int j=0; j<indices.size(); j++)
          zerovalues[j]=(Scalar)1.0e-32;
        A->replaceLocalValues(dirichletRows[i],indices,zerovalues);
      }
    }

    void Apply_BCsToMatrixCols(Teuchos::RCP<Matrix>& A,
                                      std::vector<LocalOrdinal>& dirichletCols) {
      for(size_t i=0; i<A->getNodeNumRows(); i++) {
        Teuchos::ArrayView<const LocalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> values;
        A->getLocalRowView(i,indices,values);
        std::vector<Scalar> vec;
        vec.resize(indices.size());
        Teuchos::ArrayView<Scalar> zerovalues(vec);
        for(int j=0; j<indices.size(); j++) {
          if(dirichletCols[indices[j]]==1)
            zerovalues[j]=(Scalar)1.0e-32;
          else
            zerovalues[j]=values[j];
        }
        A->replaceLocalValues(i,indices,zerovalues);
      }
    }

    void Remove_Zeroed_Rows(Teuchos::RCP<Matrix>& A, double tol=1.0e-14) {
      Teuchos::RCP<const Map> rowMap = A->getRowMap();
      RCP<Matrix> DiagMatrix = MatrixFactory::Build(rowMap,1);
      RCP<Matrix> NewMatrix  = MatrixFactory::Build(rowMap,1);
      for(size_t i=0; i<A->getNodeNumRows(); i++) {
        Teuchos::ArrayView<const LocalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> values;
        A->getLocalRowView(i,indices,values);
        int nnz=0;
        for (int j=0; j<indices.size(); j++) {
          if (Teuchos::ScalarTraits<Scalar>::magnitude(values[j]) > tol) {
            nnz++;
          }
        }
        Scalar one = (Scalar)1.0;
        Scalar zero = (Scalar)0.0;
        GlobalOrdinal row = rowMap->getGlobalElement(i);
        if (nnz == 0) {
          DiagMatrix->insertGlobalValues(row,
                                         Teuchos::ArrayView<GlobalOrdinal>(&row,1),
                                         Teuchos::ArrayView<Scalar>(&one,1));
        }
        else {
          DiagMatrix->insertGlobalValues(row,
                                         Teuchos::ArrayView<GlobalOrdinal>(&row,1),
                                         Teuchos::ArrayView<Scalar>(&zero,1));
        }
      }
      DiagMatrix->fillComplete();
      A->fillComplete();
      // add matrices together
      RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TwoMatrixAdd(*DiagMatrix,false,(Scalar)1.0,*A,false,(Scalar)1.0,NewMatrix,*out);
      NewMatrix->fillComplete();
      A=NewMatrix;
    }


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
                    const Teuchos::RCP<MultiVector> & Coords,
                    Teuchos::ParameterList& List);

    //! Two hierarchies: one for the (1,1)-block, another for the (2,2)-block
    Teuchos::RCP<Hierarchy> Hierarchy11_, Hierarchy22_, HierarchySmoother_;
    //! Top Level
    Teuchos::RCP<Level> TopLevel_;
    //! Various matrices
    Teuchos::RCP<Matrix> SM_Matrix_, D0_Matrix_, M0inv_Matrix_, M1_Matrix_, Ms_Matrix_;
    Teuchos::RCP<Matrix> TMT_Matrix_, TMT_Agg_Matrix_, P11_, A11_, A22_;
    //! Vectors for BCs
    std::vector<LocalOrdinal> BCrows_, BCcols_;
    //! Nullspace
    Teuchos::RCP<MultiVector>  Nullspace_, Coords_;
    //! Parameter lists
    Teuchos::ParameterList parameterList_, precList11_, precList22_, smootherList_;
    //! Some options
    bool disable_addon_;
    std::string mode_;

  };

} // namespace

#define MUELU_REFMAXWELL_SHORT
#endif // MUELU_REFMAXWELL_DECL_HPP
