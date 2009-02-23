#include "Tpetra_Operator.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_MultiVector.hpp"

namespace Tpetra {

  template <class Scalar, class LocalOrdinal=int, class GlobalOrdinal=LocalOrdinal> 
  class DiagPrecond : public Operator<Scalar,LocalOrdinal,GlobalOrdinal> {
    public:
      DiagPrecond(const Vector<Scalar,LocalOrdinal,GlobalOrdinal> &diags);
      DiagPrecond(const Map<LocalOrdinal,GlobalOrdinal> &rangeMap, const Map<LocalOrdinal,GlobalOrdinal> &domainMap, const Vector<Scalar,LocalOrdinal,GlobalOrdinal> &diags);
      virtual ~DiagPrecond();
      const Map<LocalOrdinal,GlobalOrdinal> & getDomainMap() const;
      const Map<LocalOrdinal,GlobalOrdinal> & getRangeMap() const;
      void apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>& X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS) const;
    private:
      typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> MV;
      void init(const Vector<Scalar,LocalOrdinal,GlobalOrdinal> &diags);
      Map<LocalOrdinal,GlobalOrdinal> rcMap_, rangeMap_, domainMap_;
      typename Teuchos::ArrayRCP<Scalar> values_;
      Teuchos::RCP<Import<LocalOrdinal,GlobalOrdinal> > importer_;
      Teuchos::RCP<Export<LocalOrdinal,GlobalOrdinal> > exporter_;
      mutable Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> > importMV_, exportMV_;
  };

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void DiagPrecond<Scalar,LocalOrdinal,GlobalOrdinal>::init(const Vector<Scalar,LocalOrdinal,GlobalOrdinal> &diags) 
  {
    Teuchos_Ordinal numLocalRows = rcMap_.getNumMyEntries();
    if (numLocalRows) {
      values_ = Teuchos::arcp<Scalar>(numLocalRows); 
    }
    diags.extractCopy1D(values_());
    // construct importer/exporter
    if (!domainMap_.isSameAs(rcMap_)) {
      importer_ = Teuchos::rcp( new Import<LocalOrdinal,GlobalOrdinal>(domainMap_,rcMap_) );
    }
    else {
      importer_ = Teuchos::null;
    }
    if (!rangeMap_.isSameAs(rcMap_)) {
      exporter_ = Teuchos::rcp( new Export<LocalOrdinal,GlobalOrdinal>(rcMap_,rangeMap_) );
    }
    else {
      exporter_ = Teuchos::null;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  DiagPrecond<Scalar,LocalOrdinal,GlobalOrdinal>::DiagPrecond(const Vector<Scalar,LocalOrdinal,GlobalOrdinal> &diags) 
  : rcMap_(diags.getMap())
  , rangeMap_(diags.getMap())
  , domainMap_(diags.getMap())
  {
    init(diags);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  DiagPrecond<Scalar,LocalOrdinal,GlobalOrdinal>::DiagPrecond(const Map<LocalOrdinal,GlobalOrdinal> &rangeMap, const Map<LocalOrdinal,GlobalOrdinal> &domainMap, const Vector<Scalar,LocalOrdinal,GlobalOrdinal> &diags)  
  : rcMap_(diags.getMap())
  , rangeMap_(rangeMap)
  , domainMap_(domainMap)
  { 
    init(diags);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  DiagPrecond<Scalar,LocalOrdinal,GlobalOrdinal>::~DiagPrecond() 
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  const Map<LocalOrdinal,GlobalOrdinal> & 
  DiagPrecond<Scalar,LocalOrdinal,GlobalOrdinal>::getRangeMap() const
  { return rangeMap_; }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  const Map<LocalOrdinal,GlobalOrdinal> & 
  DiagPrecond<Scalar,LocalOrdinal,GlobalOrdinal>::getDomainMap() const
  { return domainMap_; }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void DiagPrecond<Scalar,LocalOrdinal,GlobalOrdinal>::apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &X, 
                                                                 MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &Y, 
                                                                 Teuchos::ETransp mode) const
  {
    using Teuchos::null;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    TEST_FOR_EXCEPTION(X.numVectors() != Y.numVectors(), std::runtime_error,
        Teuchos::typeName(*this) << "::apply(X,Y): X and Y must have the same number of vectors.");
    Teuchos_Ordinal numVectors = X.numVectors();
    Teuchos_Ordinal numLocalRows = rcMap_.getNumMyEntries();
    Y = X;
    Teuchos::Array<typename MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>::const_pointer> xptrs(numVectors);
    Teuchos::Array<typename MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>::pointer      > yptrs(numVectors);
    for (Teuchos_Ordinal j=0; j<numVectors; ++j) {
      xptrs[j] = X[j];
      yptrs[j] = Y[j];
    }
    if (importer_ != null) {
      if (importMV_ != null && importMV_->numVectors() != numVectors) importMV_ = null;
      if (importMV_ == null) {
        importMV_ = Teuchos::rcp( new MV(rcMap_,numVectors) );
      }
    }
    if (exporter_ != null) {
      if (exportMV_ != null && exportMV_->numVectors() != numVectors) exportMV_ = null;
      if (exportMV_ == null) {
        exportMV_ = Teuchos::rcp( new MV(rcMap_,numVectors) );
      }
    }
    if (mode == Teuchos::NO_TRANS) {
      // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
      if (importer_ != null) {
        importMV_->doImport(X, *importer_, INSERT);
        for (Teuchos_Ordinal j=0; j<numVectors; ++j) {
          xptrs[j] = (*importMV_)[j];
        }
      }
      // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
      // We will compute solution into the to-be-exported MV; get a view
      if (exporter_ != null) {
        for (Teuchos_Ordinal j=0; j<numVectors; ++j) {
          yptrs[j] = (*exportMV_)[j];
        }
      }
      for (Teuchos_Ordinal j=0; j<numVectors; ++j) {
        typename MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>::pointer       yptr = yptrs[j];
        typename MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>::const_pointer xptr = xptrs[j];
        for (Teuchos_Ordinal i=0; i<numLocalRows; ++i) {
          yptr[i] = values_[i] * xptr[i];
        }
      }
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
    else {
      // mode == CONJ_TRANS or TRANS
      TEST_FOR_EXCEPTION(Teuchos::ScalarTraits<Scalar>::isComplex && mode == Teuchos::TRANS, std::logic_error,
          Teuchos::typeName(*this) << "::apply() does not currently support transposed multiplications for complex scalar types.");
      // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors
      if (exporter_ != null) {
        exportMV_->doImport(X,*exporter_,INSERT);
        for (Teuchos_Ordinal j=0; j<numVectors; ++j) {
          xptrs[j] = (*exportMV_)[j];
        }
      }
      // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
      // We will compute colutioni into the to-be-exported MV; get a view
      if (importer_ != null) {
        for (Teuchos_Ordinal j=0; j<numVectors; ++j) {
          yptrs[j] = (*importMV_)[j];
        }
      }
      for (Teuchos_Ordinal j=0; j<numVectors; ++j) {
        typename MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>::pointer       yptr = yptrs[j];
        typename MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>::const_pointer xptr = xptrs[j];
        for (Teuchos_Ordinal i=0; i<numLocalRows; ++i) {
          yptr[i] = ST::conjugate(values_[i]) * xptr[i];
        }
      }
      if (importer_ != null) {
        Y.putScalar(0.0); // Make sure target is zero: necessary because we are adding. may need adjusting for alpha,beta apply()
        Y.doExport(*importMV_,*importer_,ADD);
      }
      // Handle case of rangemap being a local replicated map: in this case, sum contributions from each processor
      if (Y.isDistributed() == false) {
        Y.reduce();
      }
    }
  }

};
