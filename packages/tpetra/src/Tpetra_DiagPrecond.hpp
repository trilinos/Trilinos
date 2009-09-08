#include "Tpetra_Operator.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include <Kokkos_DefaultArithmetic.hpp>

namespace Tpetra {

  template <class Scalar, class LocalOrdinal=int, class GlobalOrdinal=LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class DiagPrecond : public Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
    public:
      DiagPrecond(const Teuchos::RCP<const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &diags);
      DiagPrecond(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap, const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, const Teuchos::RCP<const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &diags);
      virtual ~DiagPrecond();
      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getOperatorDomainMap() const;
      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getOperatorRangeMap() const;
      void apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, Teuchos::ETransp trans) const;
    private:
      typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
      void init();
      Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > rcMap_, rangeMap_, domainMap_;
      Teuchos::RCP<Import<LocalOrdinal,GlobalOrdinal,Node> > importer_;
      Teuchos::RCP<Export<LocalOrdinal,GlobalOrdinal,Node> > exporter_;
      mutable Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > importMV_, exportMV_;
  };

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void DiagPrecond<Scalar,LocalOrdinal,GlobalOrdinal,Node>::init() {
    // construct importer/exporter
    if (!domainMap_->isSameAs(*rcMap_)) {
      importer_ = Teuchos::rcp( new Import<LocalOrdinal,GlobalOrdinal,Node>(domainMap_,rcMap_) );
    }
    else {
      importer_ = Teuchos::null;
    }
    if (!rangeMap_->isSameAs(*rcMap_)) {
      exporter_ = Teuchos::rcp( new Export<LocalOrdinal,GlobalOrdinal,Node>(rcMap_,rangeMap_) );
    }
    else {
      exporter_ = Teuchos::null;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  DiagPrecond<Scalar,LocalOrdinal,GlobalOrdinal,Node>::DiagPrecond(
                                      const Teuchos::RCP<const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &diags) 
  : rcMap_(diags->getMap())
  , rangeMap_(diags->getMap())
  , domainMap_(diags->getMap()) {
    init();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  DiagPrecond<Scalar,LocalOrdinal,GlobalOrdinal,Node>::DiagPrecond(
                                      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap, 
                                      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, 
                                      const Teuchos::RCP<const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &diags)  
  : rcMap_(diags->getMap())
  , rangeMap_(rangeMap)
  , domainMap_(domainMap) { 
    init();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  DiagPrecond<Scalar,LocalOrdinal,GlobalOrdinal,Node>::~DiagPrecond() 
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  DiagPrecond<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getOperatorRangeMap() const
  { return rangeMap_; }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  DiagPrecond<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getOperatorDomainMap() const
  { return domainMap_; }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void DiagPrecond<Scalar,LocalOrdinal,GlobalOrdinal,Node>::apply(
                                      const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
                                            MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
                                      Teuchos::ETransp trans) const {
    using Teuchos::null;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
        Teuchos::typeName(*this) << "::apply(X,Y): X and Y must have the same number of vectors.");
    TEST_FOR_EXCEPTION(trans != Teuchos::NO_TRANS, std::logic_error,
        Teuchos::typeName(*this) << "::apply() does not currently support transposed multiplications.");
    const size_t numVectors = X.getNumVectors();
    Teuchos::RCP<const MV> Xptr = rcpFromRef(X);
    Teuchos::RCP<      MV> Yptr = rcpFromRef(Y);
    if (importer_ != null) {
      if (importMV_ != null && importMV_->getNumVectors() != numVectors) importMV_ = null;
      if (importMV_ == null) {
        importMV_ = Teuchos::rcp( new MV(rcMap_,numVectors) );
      }
    }
    if (exporter_ != null) {
      if (exportMV_ != null && exportMV_->getNumVectors() != numVectors) exportMV_ = null;
      if (exportMV_ == null) {
        exportMV_ = Teuchos::rcp( new MV(rcMap_,numVectors) );
      }
    }
    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if (importer_ != null) {
      importMV_->doImport(X, *importer_, INSERT);
      Xptr = importMV_;
    }
    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    // We will compute solution into the to-be-exported MV; get a view
    if (exporter_ != null) {
      Yptr = exportMV_;
    }
    // compute the solution
    Yptr->multiply(ST::one(),*Yptr,*Xptr,ST::zero());
    // do the export
    if (exporter_ != null) {
      Y.putScalar(0.0);  // Make sure target is zero: necessary because we are adding. may need adjusting for alpha,beta apply()
      Y.doExport(*exportMV_, *exporter_, ADD); // Fill Y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map: in this case, sum contributions from each processor
    if (Y.isDistributed() == false) {
      Y.reduce();
    }
  }

}
