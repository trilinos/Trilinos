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
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Tpetra_RowMatrix.hpp>

#include <Xpetra_IO.hpp>

#include <MueLu.hpp>
#include "MueLu_Exceptions.hpp"
#include <MueLu_CreateXpetraPreconditioner.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp
#endif

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;

namespace Tpetra {

  template <class Scalar = Tpetra::Operator<>::scalar_type,
            class LocalOrdinal = typename Tpetra::Operator<Scalar>::local_ordinal_type,
            class GlobalOrdinal = typename Tpetra::Operator<Scalar, LocalOrdinal>::global_ordinal_type,
            class Node = typename Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
  class HierarchicalOperator : public Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> {

  public:
    using matrix_type = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using vec_type = Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
    using map_type = Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>;

    //! The RowMatrix representing the base class of CrsMatrix
    using row_matrix_type = RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

    using impl_scalar_type = typename row_matrix_type::impl_scalar_type;
    using mag_type = typename Kokkos::ArithTraits<impl_scalar_type>::mag_type;

    using local_inds_device_view_type =
          typename row_matrix_type::local_inds_device_view_type;
    using local_inds_host_view_type =
          typename row_matrix_type::local_inds_host_view_type;
    using nonconst_local_inds_host_view_type =
          typename row_matrix_type::nonconst_local_inds_host_view_type;

    using global_inds_device_view_type =
          typename row_matrix_type::global_inds_device_view_type;
    using global_inds_host_view_type =
          typename row_matrix_type::global_inds_host_view_type;
    using nonconst_global_inds_host_view_type =
          typename row_matrix_type::nonconst_global_inds_host_view_type;

    using values_device_view_type =
          typename row_matrix_type::values_device_view_type;
    using values_host_view_type =
          typename row_matrix_type::values_host_view_type;
    using nonconst_values_host_view_type =
          typename row_matrix_type::nonconst_values_host_view_type;

    //! @name Constructor/Destructor
    //@{

    //! Constructor
    HierarchicalOperator(const RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& nearField,
                         const RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& kernelApproximations,
                         const RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& basisMatrix,
                         std::vector<RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& transferMatrices)
      :
      nearField_(nearField),
      kernelApproximations_(kernelApproximations),
      basisMatrix_(basisMatrix),
      transferMatrices_(transferMatrices)
    {
      auto map = nearField_->getDomainMap();
      clusterCoeffMap_ = basisMatrix_->getDomainMap();
      TEUCHOS_ASSERT(map->isSameAs(*nearField_->getRangeMap()));
      TEUCHOS_ASSERT(map->isSameAs(*basisMatrix->getRangeMap()));
      // TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*basisMatrix->getDomainMap()));
      TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*kernelApproximations_->getDomainMap()));
      TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*kernelApproximations_->getRangeMap()));
      TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*kernelApproximations_->getRowMap()));

      for (size_t i = 0; i<transferMatrices_.size(); i++) {
        TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*transferMatrices_[i]->getDomainMap()));
        TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*transferMatrices_[i]->getRangeMap()));
      }

      allocateMemory(1);
    }

    //! Returns the Tpetra::Map object associated with the domain of this operator.
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getDomainMap() const {
      return nearField_->getDomainMap();
    }

    //! Returns the Tpetra::Map object associated with the range of this operator.
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getRangeMap() const {
      return nearField_->getRangeMap();
    }

    //! Returns in Y the result of a Tpetra::Operator applied to a Tpetra::MultiVector X.
    /*!
      \param[in]  X - Tpetra::MultiVector of dimension NumVectors to multiply with matrix.
      \param[out] Y -Tpetra::MultiVector of dimension NumVectors containing result.
    */
    void apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
               Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
               Teuchos::ETransp mode = Teuchos::NO_TRANS,
               Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
               Scalar beta  = Teuchos::ScalarTraits<Scalar>::zero()) const {
      const Scalar one = Teuchos::ScalarTraits<Scalar>::one();

      allocateMemory(X.getNumVectors());

      // near field
      nearField_->apply(X, Y, mode, alpha, beta);

      // auto out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr));
      // Y.describe(*out, Teuchos::VERB_EXTREME);

      // upward pass
      basisMatrix_->apply(X, *coefficients_, Teuchos::TRANS);

      bool flip = true;
      for (int i = Teuchos::as<int>(transferMatrices_.size())-1; i>=0; i--)
        if (flip) {
          coefficients2_->assign(*coefficients_);
          transferMatrices_[i]->apply(*coefficients_, *coefficients2_, Teuchos::NO_TRANS, one, one);
          flip = false;
        } else {
          coefficients_->assign(*coefficients2_);
          transferMatrices_[i]->apply(*coefficients2_, *coefficients_, Teuchos::NO_TRANS, one, one);
          flip = true;
        }

      if (flip)
        kernelApproximations_->apply(*coefficients_, *coefficients2_, mode, alpha);
      else
        kernelApproximations_->apply(*coefficients2_, *coefficients_, mode, alpha);
      // coefficients2_->describe(*out, Teuchos::VERB_EXTREME);

      // downward pass
      for (size_t i = 0; i<transferMatrices_.size(); i++)
        if (flip) {
          coefficients_->assign(*coefficients2_);
          transferMatrices_[i]->apply(*coefficients2_, *coefficients_, Teuchos::TRANS, one, one);
          flip = false;
        } else {
          coefficients2_->assign(*coefficients_);
          transferMatrices_[i]->apply(*coefficients_, *coefficients2_, Teuchos::TRANS, one, one);
          flip = true;
        }
      if (flip)
        basisMatrix_->apply(*coefficients2_, Y, Teuchos::NO_TRANS, one, one);
      else
        basisMatrix_->apply(*coefficients_, Y, Teuchos::NO_TRANS, one, one);
    }


    RCP<HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > restrict(const RCP<matrix_type>& P) {

      RCP<matrix_type> temp = rcp(new matrix_type(nearField_->getRowMap(), 0));
      MatrixMatrix::Multiply(*nearField_, false, *P, false, *temp);
      RCP<matrix_type> newNearField = rcp(new matrix_type(P->getDomainMap(), 0));
      MatrixMatrix::Multiply(*P, true, *temp, false, *newNearField);

      RCP<matrix_type> newBasisMatrix = rcp(new matrix_type(P->getDomainMap(), clusterCoeffMap_, 0));
      MatrixMatrix::Multiply(*P, true, *basisMatrix_, false, *newBasisMatrix);

      return rcp(new HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>(newNearField, kernelApproximations_, newBasisMatrix, transferMatrices_));
    }

    RCP<matrix_type> toMatrix() {
      const Scalar one = Teuchos::ScalarTraits<Scalar>::one();

      // construct identity on clusterCoeffMap_
      RCP<matrix_type> identity = rcp(new matrix_type(clusterCoeffMap_, 1));
      Teuchos::ArrayView<const GlobalOrdinal> gblRows = clusterCoeffMap_->getNodeElementList ();
      for (auto it = gblRows.begin (); it != gblRows.end (); ++it) {
        Teuchos::Array<GlobalOrdinal> col (1, *it);
        Teuchos::Array<Scalar> val (1, one);
        identity->insertGlobalValues (*it, col (), val ());
      }
      identity->fillComplete ();

      // transfer = basisMatrix_ * (identity + transferMatrices_[0]) * ... * (identity + transferMatrices_[n-1])
      RCP<matrix_type> transfer = rcp(new matrix_type(*basisMatrix_));
      for (size_t i = 0; i<transferMatrices_.size(); i++) {
        RCP<matrix_type> temp = MatrixMatrix::add(one, false, *identity, one, false, *transferMatrices_[i]);
        RCP<matrix_type> temp2 = rcp(new matrix_type(basisMatrix_->getRowMap(), 0));
        MatrixMatrix::Multiply(*transfer, false, *temp, true, *temp2);
        transfer = temp2;
      }

      // farField = transfer * kernelApproximations_ * transfer^T
      RCP<matrix_type> temp = rcp(new matrix_type(basisMatrix_->getRowMap(), 0));
      MatrixMatrix::Multiply(*transfer, false, *kernelApproximations_, false, *temp);
      RCP<matrix_type> farField = rcp(new matrix_type(basisMatrix_->getRowMap(), 0));
      MatrixMatrix::Multiply(*temp, false, *transfer, true, *farField);

      // nearField_ + farField
      return MatrixMatrix::add(one, false, *nearField_, one, false, *farField);

    }

    double getCompression() {
      size_t nnz = (nearField_->getGlobalNumEntries() +
                    kernelApproximations_->getGlobalNumEntries() +
                    basisMatrix_->getGlobalNumEntries());
      for (size_t i = 0; i < transferMatrices_.size(); i++)
        nnz += transferMatrices_[i]->getGlobalNumEntries();
      return Teuchos::as<double>(nnz) / (getDomainMap()->getGlobalNumElements()*getDomainMap()->getGlobalNumElements());
    }

    // Fake RowMatrix interface
    Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getRowMap() const {
      return nearField_->getRowMap();
    }

    Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getColMap() const {
      return nearField_->getColMap();
    }

    Teuchos::RCP<const Teuchos::Comm<int> > getComm() const {
      return nearField_->getDomainMap()->getComm();
    }

    Teuchos::RCP<const RowGraph<LocalOrdinal,GlobalOrdinal,Node> > getGraph() const {
      return nearField_->getCrsGraph();
    }

    global_size_t getGlobalNumRows() const {
      return nearField_->getGlobalNumRows();
    }

    global_size_t getGlobalNumCols() const {
      return nearField_->getGlobalNumCols();
    }

    size_t getNodeNumRows() const {
      return nearField_->getNodeNumRows();
    }

    size_t getNodeNumCols() const {
      return nearField_->getNodeNumCols();
    }

    GlobalOrdinal getIndexBase() const {
      return nearField_->getIndexBase();
    }

    global_size_t getGlobalNumEntries() const {
      return nearField_->getGlobalNumEntries();
    }

    size_t getNodeNumEntries() const {
      return nearField_->getNodeNumEntries();
    }

    size_t getNumEntriesInGlobalRow (GlobalOrdinal globalRow) const {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }

    size_t getNumEntriesInLocalRow (LocalOrdinal localRow) const {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }

    size_t getGlobalMaxNumRowEntries () const {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }

    size_t getNodeMaxNumRowEntries () const {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }

    bool hasColMap () const {
      return false;
    }

    bool isLocallyIndexed() const {
      return true;
    }

    bool isGloballyIndexed() const {
      return true;
    }

    bool isFillComplete() const {
      return true;
    }

    bool supportsRowViews() const {
      return false;
    }

    void
    getGlobalRowCopy (GlobalOrdinal GlobalRow,
                      nonconst_global_inds_host_view_type &Indices,
                      nonconst_values_host_view_type &Values,
                      size_t& NumEntries) const {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }

#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    void
    getGlobalRowCopy (GlobalOrdinal GlobalRow,
                      const Teuchos::ArrayView<GlobalOrdinal> &Indices,
                      const Teuchos::ArrayView<Scalar> &Values,
                      size_t &NumEntries) const {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }
#endif

    void
    getLocalRowCopy (LocalOrdinal LocalRow,
                     nonconst_local_inds_host_view_type &Indices,
                     nonconst_values_host_view_type &Values,
                     size_t& NumEntries) const {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }

#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    void
    getLocalRowCopy (LocalOrdinal LocalRow,
                     const Teuchos::ArrayView<LocalOrdinal> &Indices,
                     const Teuchos::ArrayView<Scalar> &Values,
                     size_t &NumEntries) const {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }
#endif

    void
    getGlobalRowView (GlobalOrdinal GlobalRow,
                      global_inds_host_view_type &indices,
                      values_host_view_type &values) const {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }

#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    void
    getGlobalRowView (GlobalOrdinal GlobalRow,
                      Teuchos::ArrayView<const GlobalOrdinal> &indices,
                      Teuchos::ArrayView<const Scalar> &values) const {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }
#endif

    void
    getLocalRowView (LocalOrdinal LocalRow,
                     local_inds_host_view_type & indices,
                     values_host_view_type & values) const {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }

#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    void
    getLocalRowView (LocalOrdinal LocalRow,
                     Teuchos::ArrayView<const LocalOrdinal>& indices,
                     Teuchos::ArrayView<const Scalar>& values) const {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }
#endif

    void getLocalDiagCopy (Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &diag) const {
      nearField_->getLocalDiagCopy(diag);
    }

    void leftScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }

    void rightScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }

    mag_type getFrobeniusNorm() const {
      return 0.;
    }

  private:

    void allocateMemory(size_t numVectors) const {
      if (coefficients_.is_null() || coefficients_->getNumVectors() != numVectors) {
        coefficients_  = rcp(new vec_type(clusterCoeffMap_, numVectors));
        coefficients2_ = rcp(new vec_type(clusterCoeffMap_, numVectors));
      }
    }

    RCP<matrix_type> nearField_;
    RCP<matrix_type> kernelApproximations_;
    RCP<matrix_type> basisMatrix_;
    std::vector<RCP<matrix_type> > transferMatrices_;
    RCP<const map_type> clusterCoeffMap_;
    mutable RCP<vec_type> coefficients_, coefficients2_;
  };

}

namespace Xpetra {

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class HierarchicalOperator : public TpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> {

  public:
    using tHOp = Tpetra::HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
    using map_type = Map<LocalOrdinal,GlobalOrdinal,Node>;
    using vec_type = MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
    using matrix_type = Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>;

    //! @name Constructor/Destructor
    //@{

    //! Constructor
    HierarchicalOperator(const RCP<tHOp>& op) : op_(op) { }

    HierarchicalOperator(const RCP<matrix_type>& nearField,
                         const RCP<matrix_type>& kernelApproximations,
                         const RCP<matrix_type>& basisMatrix,
                         std::vector<RCP<matrix_type> >& transferMatrices) {
      using TpCrs = TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
      using CrsWrap = CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>;

      std::vector<RCP<typename tHOp::matrix_type> > tTransferMatrices;
      for (size_t i = 0; i<transferMatrices.size(); i++) {
        auto transferT = rcp_dynamic_cast<TpCrs>(rcp_dynamic_cast<CrsWrap>(transferMatrices[i])->getCrsMatrix(), true)->getTpetra_CrsMatrixNonConst();
        tTransferMatrices.push_back(transferT);
      }

      op_ = rcp(new tHOp(rcp_dynamic_cast<TpCrs>(rcp_dynamic_cast<CrsWrap>(nearField)->getCrsMatrix(), true)->getTpetra_CrsMatrixNonConst(),
                         rcp_dynamic_cast<TpCrs>(rcp_dynamic_cast<CrsWrap>(kernelApproximations)->getCrsMatrix(), true)->getTpetra_CrsMatrixNonConst(),
                         rcp_dynamic_cast<TpCrs>(rcp_dynamic_cast<CrsWrap>(basisMatrix)->getCrsMatrix(), true)->getTpetra_CrsMatrixNonConst(),
                         tTransferMatrices));
    }

    //! Returns the Tpetra::Map object associated with the domain of this operator.
    Teuchos::RCP<const map_type> getDomainMap() const {
      return toXpetra(op_->getDomainMap());
    }

    //! Returns the Tpetra::Map object associated with the range of this operator.
    Teuchos::RCP<const map_type> getRangeMap() const {
      return toXpetra(op_->getRangeMap());
    }

    //! \brief Computes the operator-multivector application.
    /*! Loosely, performs \f$Y = \alpha \cdot A^{\textrm{mode}} \cdot X + \beta \cdot Y\f$. However, the details of operation
        vary according to the values of \c alpha and \c beta. Specifically
        - if <tt>beta == 0</tt>, apply() <b>must</b> overwrite \c Y, so that any values in \c Y (including NaNs) are ignored.
        - if <tt>alpha == 0</tt>, apply() <b>may</b> short-circuit the operator, so that any values in \c X (including NaNs) are ignored.
     */
    void apply (const vec_type& X, vec_type& Y,
                Teuchos::ETransp mode = Teuchos::NO_TRANS,
                Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
                Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const {
      op_->apply(Xpetra::toTpetra(X), Xpetra::toTpetra(Y), mode, alpha, beta);
    }

    //! Compute a residual R = B - (*this) * X
    void residual(const vec_type & X,
                  const vec_type & B,
                  vec_type& R) const {
      Tpetra::Details::residual(*op_, toTpetra(X), toTpetra(B), toTpetra(R));
    }

    RCP<HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > restrict(const RCP<matrix_type>& P) {
      return rcp(new HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>(op_->restrict(rcp_dynamic_cast<TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(rcp_dynamic_cast<CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(P)->getCrsMatrix(), true)->getTpetra_CrsMatrixNonConst())));
    }

    RCP<matrix_type> toMatrix() {
      auto tpMat = rcp(new TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(op_->toMatrix()));
      return rcp(new CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rcp_dynamic_cast<CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(tpMat)));
    }

    double getCompression() {
      return op_->getCompression();
    }

    //! Gets the operator out
    RCP<Tpetra::Operator< Scalar, LocalOrdinal, GlobalOrdinal, Node> > getOperator() { return op_; }

    RCP<const Tpetra::Operator< Scalar, LocalOrdinal, GlobalOrdinal, Node> > getOperatorConst() const { return op_; }

  private:
    RCP<tHOp> op_;
  };
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
  #include "MueLu_UseShortNames.hpp"

  std::string xmlHierachical  = "hierarchical.xml";
  std::string xmlMueLu        = "muelu.xml";
  std::string xmlAuxHierarchy = "aux.xml";
  clp.setOption("xml",    &xmlHierachical);
  clp.setOption("xmlMueLu", &xmlMueLu);
  clp.setOption("xmlAux", &xmlAuxHierarchy);

  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  using HOp = Xpetra::HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>;

  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::FancyOStream& out = *fancy;
  out.setOutputToRootOnly(0);
  bool success = true;
  const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  Teuchos::ParameterList hierachicalParams;
  Teuchos::updateParametersFromXmlFileAndBroadcast(xmlHierachical, Teuchos::Ptr<Teuchos::ParameterList>(&hierachicalParams), *comm);

  // row, domain and range map of the operator
  auto map = Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ReadMap(hierachicalParams.get<std::string>("map"), lib, comm);
  // 1-to-1 map for the cluster coefficients
  auto clusterCoeffMap = Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ReadMap(hierachicalParams.get<std::string>("coefficient map"), lib, comm);
  // overlapping map for the cluster coefficients
  auto ghosted_clusterCoeffMap = Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ReadMap(hierachicalParams.get<std::string>("ghosted coefficient map"), lib, comm);

  // near field interactions
  auto nearField = Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Read(hierachicalParams.get<std::string>("near field matrix"), map);

  // far field basis expansion coefficients
  auto basisMatrix = Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Read(hierachicalParams.get<std::string>("basis expansion coefficient matrix"), map, clusterCoeffMap, clusterCoeffMap, map);
  // far field interactions
  auto kernelApproximations = Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Read(hierachicalParams.get<std::string>("far field interaction matrix"), clusterCoeffMap, ghosted_clusterCoeffMap, clusterCoeffMap, clusterCoeffMap);

  auto transfersList = hierachicalParams.sublist("shift coefficient matrices");
  std::vector<RCP<typename HOp::matrix_type> > transferMatrices;
  for (int i = 0; i < transfersList.numParams(); i++) {
    std::string filename = transfersList.get<std::string>(std::to_string(i));
    auto transfer = Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Read(filename, clusterCoeffMap, clusterCoeffMap, clusterCoeffMap, clusterCoeffMap);
    transferMatrices.push_back(transfer);
  }

  auto op = rcp(new HOp(nearField, kernelApproximations, basisMatrix, transferMatrices));

  out << "Compression: " << op->getCompression() << " of dense matrix."<< std::endl;

  auto X_ex = Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ReadMultiVector(hierachicalParams.get<std::string>("exact solution"), map);
  auto RHS  = Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ReadMultiVector(hierachicalParams.get<std::string>("right-hand side"), map);
  auto X    = MultiVectorFactory::Build(map, 1);

  {
    op->apply(*X_ex, *X);

    // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("X.mtx", *X);
    X->update(one, *RHS, -one);
    out << "|op*X_ex - RHS| = " << X->getVector(0)->norm2() << std::endl;
    // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("diff.mtx", *X);
  }

  {
    op->apply(*X_ex, *X, Teuchos::NO_TRANS, -one);

    // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("X2.mtx", *X);
    X->update(one, *RHS, one);
    out << "|(-op)*X_ex + RHS| = " << X->getVector(0)->norm2() << std::endl;
    // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("diff2.mtx", *X);
  }

  {
    op->apply(*X_ex, *X, Teuchos::TRANS, -one);

    // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("X2.mtx", *X);
    X->update(one, *RHS, one);
    out << "|(-op^T)*X_ex + RHS| = " << X->getVector(0)->norm2() << std::endl;
    // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("diff2.mtx", *X);
  }

#ifdef HAVE_MUELU_BELOS
  {
    // Solve linear system using unpreconditioned Krylov method
    out << "\n*********************************************************\n";
    out << "Unpreconditioned Krylov method\n";
    out << "*********************************************************\n\n";

    using MV = typename HOp::vec_type;
    using OP = Belos::OperatorT<MV>;

    X->putScalar(zero);
    RCP<OP> belosOp = rcp(new Belos::XpetraOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(op));
    RCP<Belos::LinearProblem<Scalar, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<Scalar, MV, OP>(belosOp, X, RHS));

    std::string belosType = "Pseudoblock CG";
    RCP<Teuchos::ParameterList> belosList = Teuchos::parameterList();
    belosList->set("Maximum Iterations",    1000); // Maximum number of iterations allowed
    belosList->set("Convergence Tolerance", 1e-5);    // Relative convergence tolerance requested
    belosList->set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
    belosList->set("Output Frequency",      1);
    belosList->set("Output Style",          Belos::Brief);

    bool set = belosProblem->setProblem();
    if (set == false) {
      throw MueLu::Exceptions::RuntimeError("ERROR:  Belos::LinearProblem failed to set up correctly!");
    }

    // Create an iterative solver manager
    Belos::SolverFactory<Scalar, MV, OP> solverFactory;
    RCP< Belos::SolverManager<Scalar, MV, OP> > solver = solverFactory.create(belosType, belosList);
    solver->setProblem(belosProblem);

    // Perform solve
    Belos::ReturnType ret = solver->solve();
    int numIts = solver->getNumIters();

    // Get the number of iterations for this solve.
    out << "Number of iterations performed for this solve: " << numIts << std::endl;

    // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("X.mtx", *X);
    X->update(one, *X_ex, -one);
    out << "|X-X_ex| = " << X->getVector(0)->norm2() << std::endl << std::endl;

    success &= (ret == Belos::Converged);

  }
#endif // HAVE_MUELU_BELOS

  {
    // Solve linear system using a AMG preconditioned Krylov method

    ////////////////////////////////////////////////////////////////
    // Build the auxiliary hierachy
    out << "\n*********************************************************\n";
    out << "Building the auxiliary hierachy\n";
    out << "*********************************************************\n\n";

    auto auxOp  = Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Read(hierachicalParams.get<std::string>("auxiliary operator"), map);
    auto coords = Xpetra::IO<typename Teuchos::ScalarTraits<Scalar>::coordinateType,LocalOrdinal,GlobalOrdinal,Node>::ReadMultiVector(hierachicalParams.get<std::string>("coordinates"), map);

    Teuchos::ParameterList auxParams;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlAuxHierarchy, Teuchos::Ptr<Teuchos::ParameterList>(&auxParams), *comm);
    auxParams.sublist("user data").set("Coordinates", coords);

    auto auxH = MueLu::CreateXpetraPreconditioner(auxOp, auxParams);

    ////////////////////////////////////////////////////////////////
    // Construct the main hierarchy
    out << "\n*********************************************************\n";
    out << "Building the main hierachy\n";
    out << "*********************************************************\n\n";

    Teuchos::ParameterList params;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlMueLu, Teuchos::Ptr<Teuchos::ParameterList>(&params), *comm);
    params.set("coarse: max size", 1);
    params.set("max levels", auxH->GetNumLevels());

    auto H = rcp(new Hierarchy());
    RCP<Level> lvl = H->GetLevel(0);
    lvl->Set("A", rcp_dynamic_cast<Operator>(op));
    lvl->Set("Coordinates", coords);
    for(int lvlNo = 1; lvlNo<auxH->GetNumLevels(); lvlNo++) {
      H->AddNewLevel();
      RCP<Level> auxLvl = auxH->GetLevel(lvlNo);
      // auto mgr = auxLvl->GetFactoryManager();
      // auxLvl->print(std::cout, MueLu::Debug);
      RCP<Level> fineLvl = H->GetLevel(lvlNo-1);
      lvl = H->GetLevel(lvlNo);
      auto P = auxLvl->Get<RCP<Matrix> >("P");
      lvl->Set("P", P);
      params.sublist("level "+std::to_string(lvlNo)).set("P", P);

      auto fineA = rcp_dynamic_cast<HOp>(fineLvl->Get<RCP<Operator> >("A"));
      auto coarseA = fineA->restrict(P);
      if (lvlNo+1 == auxH->GetNumLevels())
        lvl->Set("A", coarseA->toMatrix());
      else
        lvl->Set("A", rcp_dynamic_cast<Operator>(coarseA));
    }

    RCP<HierarchyManager> mueLuFactory = rcp(new ParameterListInterpreter(params,op->getDomainMap()->getComm()));
    H->setlib(op->getDomainMap()->lib());
    H->SetProcRankVerbose(op->getDomainMap()->getComm()->getRank());
    mueLuFactory->SetupHierarchy(*H);
    H->IsPreconditioner(true);


#ifdef HAVE_MUELU_BELOS
    ////////////////////////////////////////////////////////////////
    // Set up the Krylov solver

    using MV = typename HOp::vec_type;
    using OP = Belos::OperatorT<MV>;

    X->putScalar(zero);
    RCP<OP> belosOp = rcp(new Belos::XpetraOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(op));
    RCP<OP> belosPrec = rcp(new Belos::MueLuOp <SC, LO, GO, NO>(H));
    RCP<Belos::LinearProblem<Scalar, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<Scalar, MV, OP>(belosOp, X, RHS));

    std::string belosType = "Pseudoblock CG";
    RCP<Teuchos::ParameterList> belosList = Teuchos::parameterList();
    belosList->set("Maximum Iterations",    1000); // Maximum number of iterations allowed
    belosList->set("Convergence Tolerance", 1e-5);    // Relative convergence tolerance requested
    belosList->set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
    belosList->set("Output Frequency",      1);
    belosList->set("Output Style",          Belos::Brief);

    belosProblem->setRightPrec(belosPrec);

    bool set = belosProblem->setProblem();
    if (set == false) {
      throw MueLu::Exceptions::RuntimeError("ERROR:  Belos::LinearProblem failed to set up correctly!");
    }

    // Create an iterative solver manager
    Belos::SolverFactory<Scalar, MV, OP> solverFactory;
    RCP< Belos::SolverManager<Scalar, MV, OP> > solver = solverFactory.create(belosType, belosList);
    solver->setProblem(belosProblem);

    // Perform solve
    Belos::ReturnType ret = solver->solve();
    int numIts = solver->getNumIters();

    // Get the number of iterations for this solve.
    out << "Number of iterations performed for this solve: " << numIts << std::endl;

    // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("X.mtx", *X);
    X->update(one, *X_ex, -one);
    out << "|X-X_ex| = " << X->getVector(0)->norm2() << std::endl;

    success &= (ret == Belos::Converged);

#endif // HAVE_MUELU_BELOS
  }

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} //main

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc,argv);
}
