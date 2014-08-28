// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
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
#include <Teuchos_Array.hpp>
#include "Xpetra_EpetraCrsMatrix.hpp"

namespace Xpetra {

  template<class EpetraGlobalOrdinal>
  EpetraCrsMatrixT<EpetraGlobalOrdinal>::EpetraCrsMatrixT(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : mtx_(Teuchos::rcp(new Epetra_CrsMatrix(Copy, toEpetra(rowMap), maxNumEntriesPerRow, toEpetra(pftype)))), isFillResumed_(false) { }

  template<class EpetraGlobalOrdinal>
  EpetraCrsMatrixT<EpetraGlobalOrdinal>::EpetraCrsMatrixT(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : isFillResumed_(false)
  {
    Teuchos::Array<int> numEntriesPerRowToAlloc(NumEntriesPerRowToAlloc.begin(), NumEntriesPerRowToAlloc.end()); // convert array of "size_t" to array of "int"
    mtx_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, toEpetra(rowMap), numEntriesPerRowToAlloc.getRawPtr(), toEpetra(pftype)));
  }

  template<class EpetraGlobalOrdinal>
  EpetraCrsMatrixT<EpetraGlobalOrdinal>::EpetraCrsMatrixT(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : mtx_(Teuchos::rcp(new Epetra_CrsMatrix(Copy, toEpetra(rowMap), toEpetra(colMap), maxNumEntriesPerRow, toEpetra(pftype)))), isFillResumed_(false) { }

  template<class EpetraGlobalOrdinal>
  EpetraCrsMatrixT<EpetraGlobalOrdinal>::EpetraCrsMatrixT(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : isFillResumed_(false)
  {
    Teuchos::Array<int> numEntriesPerRowToAlloc(NumEntriesPerRowToAlloc.begin(), NumEntriesPerRowToAlloc.end()); // convert array of "size_t" to array of "int"
    mtx_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, toEpetra(rowMap), toEpetra(colMap), numEntriesPerRowToAlloc.getRawPtr(), toEpetra(pftype)));
  }

  template<class EpetraGlobalOrdinal>
  EpetraCrsMatrixT<EpetraGlobalOrdinal>::EpetraCrsMatrixT(const Teuchos::RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node > > &graph, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : mtx_(Teuchos::rcp(new Epetra_CrsMatrix(Copy, toEpetra(graph)))), isFillResumed_(false) { }

  template<class EpetraGlobalOrdinal>
  EpetraCrsMatrixT<EpetraGlobalOrdinal>::EpetraCrsMatrixT(const EpetraCrsMatrixT& matrix)
    : mtx_(Teuchos::rcp(new Epetra_CrsMatrix(*(matrix.mtx_)))), isFillResumed_(false) { }


  template<class EpetraGlobalOrdinal>
  EpetraCrsMatrixT<EpetraGlobalOrdinal>::EpetraCrsMatrixT(const Teuchos::RCP<const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& sourceMatrix,
                                   const Import<LocalOrdinal,GlobalOrdinal,Node> &importer,
                                   const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
                                   const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
                                   const Teuchos::RCP<Teuchos::ParameterList>& params):
    isFillResumed_(false)
  {
    XPETRA_DYNAMIC_CAST(const EpetraCrsMatrixT<GlobalOrdinal>, *sourceMatrix, tSourceMatrix, "Xpetra::EpetraCrsMatrixT constructor only accepts Xpetra::EpetraCrsMatrixT as an input argument.");
    XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal>, importer, tImporter, "Xpetra::EpetraCrsMatrixT constructor only accepts Xpetra::EpetraImportT as an input argument.");

    const Epetra_Map* myDomainMap = (domainMap!=Teuchos::null)? &toEpetra(domainMap): 0;
    const Epetra_Map* myRangeMap  = (rangeMap !=Teuchos::null)? &toEpetra(rangeMap) : 0;

    // Follows the Tpetra parameters
    bool restrictComm=false;
    if(!params.is_null()) restrictComm = params->get("Restrict Communicator",restrictComm);
    mtx_ = Teuchos::rcp(new Epetra_CrsMatrix(*tSourceMatrix.getEpetra_CrsMatrix(),*tImporter.getEpetra_Import(),myDomainMap,myRangeMap,restrictComm));
    if(restrictComm && mtx_->NumMyRows()==0)
      mtx_=Teuchos::null;
  }

  template<class EpetraGlobalOrdinal>
  EpetraCrsMatrixT<EpetraGlobalOrdinal>::EpetraCrsMatrixT(const Teuchos::RCP<const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& sourceMatrix,
                                   const Export<LocalOrdinal,GlobalOrdinal,Node> &exporter,
                                   const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
                                   const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
                                   const Teuchos::RCP<Teuchos::ParameterList>& params):
    isFillResumed_(false)
  {
    XPETRA_DYNAMIC_CAST(const EpetraCrsMatrixT<GlobalOrdinal>, *sourceMatrix, tSourceMatrix, "Xpetra::EpetraCrsMatrixT constructor only accepts Xpetra::EpetraCrsMatrixT as an input argument.");
    XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal>, exporter, tExporter, "Xpetra::EpetraCrsMatrixT constructor only accepts Xpetra::EpetraExportT as an input argument.");

    const Epetra_Map* myDomainMap = (domainMap!=Teuchos::null)? &toEpetra(domainMap): 0;
    const Epetra_Map* myRangeMap  = (rangeMap !=Teuchos::null)? &toEpetra(rangeMap) : 0;

    // Follows the Tpetra parameters
    bool restrictComm=false;
    if(!params.is_null()) restrictComm = params->get("Restrict Communicator",restrictComm);

    mtx_ = Teuchos::rcp(new Epetra_CrsMatrix(*tSourceMatrix.getEpetra_CrsMatrix(),*tExporter.getEpetra_Export(),myDomainMap,myRangeMap,restrictComm));
  }



  template<class EpetraGlobalOrdinal>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal>::insertGlobalValues(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &cols, const ArrayView<const double> &vals) {
    XPETRA_MONITOR("EpetraCrsMatrixT::insertGlobalValues");
    XPETRA_ERR_CHECK(mtx_->InsertGlobalValues(globalRow, vals.size(), vals.getRawPtr(), cols.getRawPtr()));
  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal>::insertLocalValues(int localRow, const ArrayView<const int> &cols, const ArrayView<const double> &vals) {
    XPETRA_MONITOR("EpetraCrsMatrixT::insertLocalValues");
    XPETRA_ERR_CHECK(mtx_->InsertMyValues(localRow, vals.size(), vals.getRawPtr(), cols.getRawPtr()));
  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal>::replaceGlobalValues(GlobalOrdinal globalRow, const ArrayView< const GlobalOrdinal > &indices, const ArrayView< const Scalar > &values) {
    XPETRA_MONITOR("EpetraCrsMatrixT::replaceGlobalValues");

    {
      const std::string tfecfFuncName("replaceGlobalValues");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(! isFillActive(), std::runtime_error,
                                            ": Fill must be active in order to call this method.  If you have already "
                                            "called fillComplete(), you need to call resumeFill() before you can "
                                            "replace values.");

      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(values.size() != indices.size(),
                                            std::runtime_error, ": values.size() must equal indices.size().");
    }

    XPETRA_ERR_CHECK(mtx_->ReplaceGlobalValues(globalRow, indices.size(), values.getRawPtr(), indices.getRawPtr()));

  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal>::replaceLocalValues(LocalOrdinal localRow, const ArrayView< const LocalOrdinal > &indices, const ArrayView< const Scalar > &values) {
    XPETRA_MONITOR("EpetraCrsMatrixT::replaceLocalValues");

    {
      const std::string tfecfFuncName("replaceLocalValues");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(! isFillActive(), std::runtime_error,
                                            ": Fill must be active in order to call this method.  If you have already "
                                            "called fillComplete(), you need to call resumeFill() before you can "
                                            "replace values.");

      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(values.size() != indices.size(),
                                            std::runtime_error, ": values.size() must equal indices.size().");
    }

    XPETRA_ERR_CHECK(mtx_->ReplaceMyValues(localRow, indices.size(), values.getRawPtr(), indices.getRawPtr()));

  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal>::allocateAllValues(size_t numNonZeros, ArrayRCP<size_t>& rowptr, ArrayRCP<int>& colind, ArrayRCP<double>& values) {
     XPETRA_MONITOR("EpetraCrsMatrixT::allocateAllValues");

    // Row offsets
    // Unfortunately, we cannot do this in the same manner as column indices
    // and values (see below).  The problem is that Tpetra insists on using
    // size_t, and Epetra uses int internally.  So we only resize here, and
    // will need to copy in setAllValues
    rowptr.resize(getNodeNumRows()+1);

    int  lowerOffset = 0;
    bool ownMemory   = false;

    // Column indices
    // Extract, resize, set colind
    Epetra_IntSerialDenseVector& myColind = mtx_->ExpertExtractIndices();
    myColind.Resize(numNonZeros);
    colind = Teuchos::arcp(myColind.Values(), lowerOffset, numNonZeros, ownMemory);

    // Values
    // Extract, reallocate, set values
    double *& myValues = mtx_->ExpertExtractValues();
    delete [] myValues;
    myValues = new double[numNonZeros];
    values = Teuchos::arcp(myValues,lowerOffset,numNonZeros,ownMemory);
  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal>::setAllValues(const ArrayRCP<size_t>& rowptr, const ArrayRCP<int>& colind, const ArrayRCP<double>& values) {
    XPETRA_MONITOR("EpetraCrsMatrixT::setAllValues");

    // Check sizes
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(rowptr.size()) != getNodeNumRows()+1, Xpetra::Exceptions::RuntimeError,
                               "An exception is thrown to let you know that the size of your rowptr array is incorrect.");
    TEUCHOS_TEST_FOR_EXCEPTION(values.size() != colind.size(), Xpetra::Exceptions::RuntimeError,
                               "An exception is thrown to let you know that you mismatched your pointers.");

    // Check pointers
    if (values.size() > 0) {
      TEUCHOS_TEST_FOR_EXCEPTION(colind.getRawPtr() != mtx_->ExpertExtractIndices().Values(), Xpetra::Exceptions::RuntimeError,
                                 "An exception is thrown to let you know that you mismatched your pointers.");
      TEUCHOS_TEST_FOR_EXCEPTION(values.getRawPtr() != mtx_->ExpertExtractValues(), Xpetra::Exceptions::RuntimeError,
                                 "An exception is thrown to let you know that you mismatched your pointers.");
    }

    // We have to make a copy here, it is unavoidable
    // See comments in allocateAllValues
    const size_t N = getNodeNumRows();

    Epetra_IntSerialDenseVector& myRowptr = mtx_->ExpertExtractIndexOffset();
    myRowptr.Resize(N+1);
    for (size_t i = 0; i < N+1; i++)
      myRowptr[i] = Teuchos::as<int>(rowptr[i]);
  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal>::getAllValues(ArrayRCP<const size_t>& rowptr, ArrayRCP<const LocalOrdinal>& colind, ArrayRCP<const Scalar>& values) const {
    XPETRA_MONITOR("EpetraCrsMatrixT::getAllValues");

    int  lowerOffset = 0;
    bool ownMemory   = false;

    const size_t n   = getNodeNumRows();
    const size_t nnz = getNodeNumEntries();

    // Row offsets
    // We have to make a copy here, it is unavoidable (see comments in allocateAllValues)
    Epetra_IntSerialDenseVector& myRowptr = mtx_->ExpertExtractIndexOffset();
    rowptr.resize(n+1);
    for (size_t i = 0; i < n+1; i++)
      (*const_cast<size_t*>(&rowptr[i])) = Teuchos::as<size_t>(myRowptr[i]);

    // Column indices
    colind = Teuchos::arcp(mtx_->ExpertExtractIndices().Values(), lowerOffset, nnz, ownMemory);

    // Values
    values = Teuchos::arcp(mtx_->ExpertExtractValues(), lowerOffset, nnz, ownMemory);
  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal>::resumeFill(const RCP< ParameterList > &params) {
    XPETRA_MONITOR("EpetraCrsMatrixT::resumeFill");

    // According to Tpetra documentation, resumeFill() may be called repeatedly.
    isFillResumed_ = true;
  }

  template<class EpetraGlobalOrdinal>
  bool EpetraCrsMatrixT<EpetraGlobalOrdinal>::isFillComplete() const { XPETRA_MONITOR("EpetraCrsMatrixT::isFillComplete"); if (isFillResumed_) return false; else return mtx_->Filled(); }

  template<class EpetraGlobalOrdinal>
  bool EpetraCrsMatrixT<EpetraGlobalOrdinal>::isFillActive() const { XPETRA_MONITOR("EpetraCrsMatrixT::isFillActive"); return !isFillComplete(); }

  template<class EpetraGlobalOrdinal>
  bool EpetraCrsMatrixT<EpetraGlobalOrdinal>::supportsRowViews() const { XPETRA_MONITOR("EpetraCrsMatrixT::supportsRowViews"); return true; }

  //TODO: throw same exception as Tpetra
  template<class EpetraGlobalOrdinal>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal>::getLocalRowCopy(int LocalRow, const ArrayView<int> &Indices, const ArrayView<double> &Values, size_t &NumEntries) const {
    XPETRA_MONITOR("EpetraCrsMatrixT::getLocalRowCopy");

    int numEntries = -1;
    XPETRA_ERR_CHECK(mtx_->ExtractMyRowCopy(LocalRow, Indices.size(), numEntries, Values.getRawPtr(), Indices.getRawPtr()));
    NumEntries = numEntries;
  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal>::getGlobalRowCopy(GlobalOrdinal GlobalRow, const ArrayView<GlobalOrdinal> &Indices, const ArrayView<double> &Values, size_t &NumEntries) const {
    XPETRA_MONITOR("EpetraCrsMatrixT::getGlobalRowCopy");

    int numEntries = -1;
    XPETRA_ERR_CHECK(mtx_->ExtractGlobalRowCopy(GlobalRow, Indices.size(), numEntries, Values.getRawPtr(), Indices.getRawPtr()));
    NumEntries = numEntries;
  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal>::getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &indices, ArrayView<const double> &values) const {
    XPETRA_MONITOR("EpetraCrsMatrixT::getGlobalRowView");

    int      numEntries;
    double * eValues;
    GlobalOrdinal * eIndices;

    XPETRA_ERR_CHECK(mtx_->ExtractGlobalRowView(GlobalRow, numEntries, eValues, eIndices));
    if (numEntries == 0) { eValues = NULL; eIndices = NULL; } // Cf. TEUCHOS_TEST_FOR_EXCEPT( p == 0 && size_in != 0 ) in Teuchos ArrayView constructor.

    indices = ArrayView<const GlobalOrdinal>(eIndices, numEntries);
    values  = ArrayView<const double>(eValues, numEntries);
  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal>::getLocalRowView(int LocalRow, ArrayView<const int> &indices, ArrayView<const double> &values) const {
    XPETRA_MONITOR("EpetraCrsMatrixT::getLocalRowView");

    int      numEntries;
    double * eValues;
    int    * eIndices;

    XPETRA_ERR_CHECK(mtx_->ExtractMyRowView(LocalRow, numEntries, eValues, eIndices));
    if (numEntries == 0) { eValues = NULL; eIndices = NULL; } // Cf. TEUCHOS_TEST_FOR_EXCEPT( p == 0 && size_in != 0 ) in Teuchos ArrayView constructor.

    indices = ArrayView<const int>(eIndices, numEntries);
    values  = ArrayView<const double>(eValues, numEntries);
  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal>::apply(const MultiVector<double, int, GlobalOrdinal> & X, MultiVector<double, int, GlobalOrdinal> &Y, Teuchos::ETransp mode, double alpha, double beta) const {
    XPETRA_MONITOR("EpetraCrsMatrixT::apply");

    //TEUCHOS_TEST_FOR_EXCEPTION((alpha != 1) || (beta != 0), Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraCrsMatrixT.multiply() only accept alpha==1 and beta==0");

    XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT<GlobalOrdinal>, X, eX, "Xpetra::EpetraCrsMatrixT->apply() only accept Xpetra::EpetraMultiVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(      EpetraMultiVectorT<GlobalOrdinal>, Y, eY, "Xpetra::EpetraCrsMatrixT->apply() only accept Xpetra::EpetraMultiVectorT as input arguments.");

    TEUCHOS_TEST_FOR_EXCEPTION((mode != Teuchos::NO_TRANS) && (mode != Teuchos::TRANS), Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraCrsMatrixT->apply() only accept mode == NO_TRANS or mode == TRANS");
    bool eTrans = toEpetra(mode);

    // /!\ UseTranspose value
    TEUCHOS_TEST_FOR_EXCEPTION(mtx_->UseTranspose(), Xpetra::Exceptions::NotImplemented, "An exception is throw to let you know that Xpetra::EpetraCrsMatrixT->apply() do not take into account the UseTranspose() parameter of Epetra_CrsMatrix.");

    RCP<Epetra_MultiVector> epY = eY.getEpetra_MultiVector();

    // helper vector: tmp = A*x
    RCP<Epetra_MultiVector> tmp = Teuchos::rcp(new Epetra_MultiVector(*epY));
    tmp->PutScalar(0.0);
    XPETRA_ERR_CHECK(mtx_->Multiply(eTrans, *eX.getEpetra_MultiVector(), *tmp));

    // calculate alpha * A * x + beta * y
    XPETRA_ERR_CHECK(eY.getEpetra_MultiVector()->Update(alpha,*tmp,beta));
  }

  template<class EpetraGlobalOrdinal>
  std::string EpetraCrsMatrixT<EpetraGlobalOrdinal>::description() const {
    XPETRA_MONITOR("EpetraCrsMatrixT::description");

    // This implementation come from Tpetra_CrsMatrix_def.hpp (without modification)
    std::ostringstream oss;
    //TODO: oss << DistObject<char, LocalOrdinal,GlobalOrdinal>::description();
    if (isFillComplete()) {
      oss << "{status = fill complete"
          << ", global rows = " << getGlobalNumRows()
          << ", global cols = " << getGlobalNumCols()
          << ", global num entries = " << getGlobalNumEntries()
          << "}";
    }
    else {
      oss << "{status = fill not complete"
          << ", global rows = " << getGlobalNumRows()
          << "}";
    }
    return oss.str();

  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
    XPETRA_MONITOR("EpetraCrsMatrixT::describe");

    // This implementation come from Tpetra_CrsMatrix_def.hpp (without modification)
    using std::endl;
    using std::setw;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    Teuchos::EVerbosityLevel vl = verbLevel;
    if (vl == VERB_DEFAULT) vl = VERB_LOW;
    RCP<const Comm<int> > comm = this->getComm();
    const int myImageID = comm->getRank(),
      numImages = comm->getSize();
    size_t width = 1;
    for (size_t dec=10; dec<getGlobalNumRows(); dec *= 10) {
      ++width;
    }
    width = std::max<size_t>(width,11) + 2;
    Teuchos::OSTab tab(out);
    //    none: print nothing
    //     low: print O(1) info from node 0
    //  medium: print O(P) info, num entries per node
    //    high: print O(N) info, num entries per row
    // extreme: print O(NNZ) info: print indices and values
    //
    // for medium and higher, print constituent objects at specified verbLevel
    if (vl != VERB_NONE) {
      if (myImageID == 0) out << this->description() << std::endl;
      // O(1) globals, minus what was already printed by description()
      if (isFillComplete() && myImageID == 0) {
        out << "Global number of diagonals = " << getGlobalNumDiags() << std::endl;
        out << "Global max number of entries = " << getGlobalMaxNumRowEntries() << std::endl;
      }
      // constituent objects
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        if (myImageID == 0) out << "\nRow map: " << std::endl;
        getRowMap()->describe(out,vl);
        //
        if (getColMap() != null) {
          if (getColMap() == getRowMap()) {
            if (myImageID == 0) out << "\nColumn map is row map.";
          }
          else {
            if (myImageID == 0) out << "\nColumn map: " << std::endl;
            getColMap()->describe(out,vl);
          }
        }
        if (getDomainMap() != null) {
          if (getDomainMap() == getRowMap()) {
            if (myImageID == 0) out << "\nDomain map is row map.";
          }
          else if (getDomainMap() == getColMap()) {
            if (myImageID == 0) out << "\nDomain map is row map.";
          }
          else {
            if (myImageID == 0) out << "\nDomain map: " << std::endl;
            getDomainMap()->describe(out,vl);
          }
        }
        if (getRangeMap() != null) {
          if (getRangeMap() == getDomainMap()) {
            if (myImageID == 0) out << "\nRange map is domain map." << std::endl;
          }
          else if (getRangeMap() == getRowMap()) {
            if (myImageID == 0) out << "\nRange map is row map." << std::endl;
          }
          else {
            if (myImageID == 0) out << "\nRange map: " << std::endl;
            getRangeMap()->describe(out,vl);
          }
        }
        if (myImageID == 0) out << std::endl;
      }
      // O(P) data
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            out << "Node ID = " << imageCtr << std::endl;
            // TODO: need a graph
            //               if (staticGraph_->indicesAreAllocated() == false) {
            //                 out << "Node not allocated" << std::endl;
            //               }
            //               else {
            //                 out << "Node number of allocated entries = " << staticGraph_->getNodeAllocationSize() << std::endl;
            //               }

            // TMP:
            //            const Epetra_CrsGraph & staticGraph_ = mtx_->Graph();
            // End of TMP

            out << "Node number of entries = " << getNodeNumEntries() << std::endl;
            if (isFillComplete()) {
              out << "Node number of diagonals = " << getNodeNumDiags() << std::endl;
            }
            out << "Node max number of entries = " << getNodeMaxNumRowEntries() << std::endl;
          }
          comm->barrier();
          comm->barrier();
          comm->barrier();
        }
      }
      // O(N) and O(NNZ) data
      if (vl == VERB_HIGH || vl == VERB_EXTREME) {
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            out << std::setw(width) << "Node ID"
                << std::setw(width) << "Global Row"
                << std::setw(width) << "Num Entries";
            if (vl == VERB_EXTREME) {
              out << std::setw(width) << "(Index,Value)";
            }
            out << std::endl;
            for (size_t r=0; r < getNodeNumRows(); ++r) {
              const size_t nE = getNumEntriesInLocalRow(r);
              GlobalOrdinal gid = getRowMap()->getGlobalElement(r);
              out << std::setw(width) << myImageID
                  << std::setw(width) << gid
                  << std::setw(width) << nE;
              if (vl == VERB_EXTREME) {
                if (isGloballyIndexed()) {
                  ArrayView<const GlobalOrdinal> rowinds;
                  ArrayView<const Scalar> rowvals;
                  getGlobalRowView(gid,rowinds,rowvals);
                  for (size_t j=0; j < nE; ++j) {
                    out << " (" << rowinds[j]
                        << ", " << rowvals[j]
                        << ") ";
                  }
                }
                else if (isLocallyIndexed()) {
                  ArrayView<const LocalOrdinal> rowinds;
                  ArrayView<const Scalar> rowvals;
                  getLocalRowView(r,rowinds,rowvals);
                  for (size_t j=0; j < nE; ++j) {
                    out << " (" << getColMap()->getGlobalElement(rowinds[j])
                        << ", " << rowvals[j]
                        << ") ";
                  }
                }
              }
              out << std::endl;
            }
          }
          comm->barrier();
          comm->barrier();
          comm->barrier();
        }
      }
    }

  }

  // TODO: use toEpetra()
  template<class EpetraGlobalOrdinal>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal>::doImport(const DistObject<char, int, GlobalOrdinal> &source,
                                 const Import<int, GlobalOrdinal> &importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsMatrixT::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsMatrixT<GlobalOrdinal>, source, tSource, "Xpetra::EpetraCrsMatrixT::doImport only accept Xpetra::EpetraCrsMatrixT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal>, importer, tImporter, "Xpetra::EpetraCrsMatrixT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    RCP<const Epetra_CrsMatrix> v = tSource.getEpetra_CrsMatrix();
    int err = mtx_->Import(*v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal>::doExport(const DistObject<char, int, GlobalOrdinal> &dest,
                                 const Import<int, GlobalOrdinal>& importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsMatrixT::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsMatrixT<GlobalOrdinal>, dest, tDest, "Xpetra::EpetraCrsMatrixT::doImport only accept Xpetra::EpetraCrsMatrixT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal>, importer, tImporter, "Xpetra::EpetraCrsMatrixT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    RCP<const Epetra_CrsMatrix> v = tDest.getEpetra_CrsMatrix();
    int err = mtx_->Export(*v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal>::doImport(const DistObject<char, int, GlobalOrdinal> &source,
                                 const Export<int, GlobalOrdinal>& exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsMatrixT::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsMatrixT<GlobalOrdinal>, source, tSource, "Xpetra::EpetraCrsMatrixT::doImport only accept Xpetra::EpetraCrsMatrixT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal>, exporter, tExporter, "Xpetra::EpetraCrsMatrixT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    RCP<const Epetra_CrsMatrix> v = tSource.getEpetra_CrsMatrix();
    int err = mtx_->Import(*v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");

  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal>::doExport(const DistObject<char, int, GlobalOrdinal> &dest,
                                 const Export<int, GlobalOrdinal>& exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsMatrixT::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsMatrixT<GlobalOrdinal>, dest, tDest, "Xpetra::EpetraCrsMatrixT::doImport only accept Xpetra::EpetraCrsMatrixT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal>, exporter, tExporter, "Xpetra::EpetraCrsMatrixT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    RCP<const Epetra_CrsMatrix> v = tDest.getEpetra_CrsMatrix();
    int err = mtx_->Export(*v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

    template<class EpetraGlobalOrdinal>
    void EpetraCrsMatrixT<EpetraGlobalOrdinal>::fillComplete(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &domainMap,
                                       const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rangeMap,
                                       const RCP< ParameterList > &params) {
      XPETRA_MONITOR("EpetraCrsMatrixT::fillComplete");

      // For Epetra matrices, resumeFill() is a fictive operation. There is no need for a fillComplete after some resumeFill() operations.
      if (isFillResumed_ == true) { isFillResumed_ = false; return; }

      bool doOptimizeStorage = true;
      if (params != null && params->get("Optimize Storage",true) == false) doOptimizeStorage = false;
      mtx_->FillComplete(toEpetra(domainMap), toEpetra(rangeMap), doOptimizeStorage);
    }

    template<class EpetraGlobalOrdinal>
    void EpetraCrsMatrixT<EpetraGlobalOrdinal>::fillComplete(const RCP< ParameterList > &params) {
      XPETRA_MONITOR("EpetraCrsMatrixT::fillComplete");

      // For Epetra matrices, resumeFill() is a fictive operation. There is no need for a fillComplete after some resumeFill() operations.
      if (isFillResumed_ == true) { isFillResumed_ = false; return; }

      bool doOptimizeStorage = true;
      if (params != null && params->get("Optimize Storage",true) == false) doOptimizeStorage = false;
      mtx_->FillComplete(doOptimizeStorage);
    }


  template<class EpetraGlobalOrdinal>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal>::replaceDomainMapAndImporter(const Teuchos::RCP< const  Map< LocalOrdinal, GlobalOrdinal, Node > >& newDomainMap, Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> >  & newImporter) {
      XPETRA_MONITOR("EpetraCrsMatrixT::replaceDomainMapAndImporter");
      XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal>, *newImporter, eImporter, "Xpetra::EpetraCrsMatrixT::replaceDomainMapAndImporter only accepts Xpetra::EpetraImportT.");

      const RCP<const Epetra_Import> & myImport = eImporter.getEpetra_Import();
      int rv=0;
      if(myImport==Teuchos::null)
        rv=mtx_->ReplaceDomainMapAndImporter( toEpetra(newDomainMap),0);
      else
        rv=mtx_->ReplaceDomainMapAndImporter( toEpetra(newDomainMap),&*myImport);
      TEUCHOS_TEST_FOR_EXCEPTION(rv != 0, std::runtime_error, "Xpetra::EpetraCrsMatrixT::replaceDomainMapAndImporter FAILED!");
  }


  template<class EpetraGlobalOrdinal>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal>::expertStaticFillComplete(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & domainMap,
                                                 const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & rangeMap,
                                                 const RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > &importer,
                                                 const RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > &exporter,
                                                 const RCP<ParameterList> & params) {
    XPETRA_MONITOR("EpetraCrsMatrixT::expertStaticFillComplete");
    int rv=0;
    const Epetra_Import * myimport =0;
    const Epetra_Export * myexport =0;

    if(!importer.is_null()) {
      XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal>, *importer, eImporter, "Xpetra::EpetraCrsMatrixT::expertStaticFillComplete only accepts Xpetra::EpetraImportT.");
      myimport = eImporter.getEpetra_Import().getRawPtr();
    }
    if(!exporter.is_null()) {
      XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal>, *exporter, eExporter, "Xpetra::EpetraCrsMatrixT::expertStaticFillComplete only accepts Xpetra::EpetraImportT.");
      myexport = eExporter.getEpetra_Export().getRawPtr();
    }

    rv=mtx_->ExpertStaticFillComplete(toEpetra(domainMap), toEpetra(rangeMap), myimport, myexport);

    TEUCHOS_TEST_FOR_EXCEPTION(rv != 0, std::runtime_error, "Xpetra::EpetraCrsMatrixT::expertStaticFillComplete FAILED!");
  }

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
template class EpetraCrsMatrixT<int>;
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
template class EpetraCrsMatrixT<long long>;
#endif

}
