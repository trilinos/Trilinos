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

  template<class EpetraGlobalOrdinal, class Node>
  EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::EpetraCrsMatrixT(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : mtx_(Teuchos::rcp(new Epetra_CrsMatrix(Copy, toEpetra<EpetraGlobalOrdinal,Node>(rowMap), maxNumEntriesPerRow, toEpetra(pftype)))), isFillResumed_(false)
#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
      , isInitializedLocalMatrix_(false)
#endif
  { }

  template<class EpetraGlobalOrdinal, class Node>
  EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::EpetraCrsMatrixT(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : isFillResumed_(false)
#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
      , isInitializedLocalMatrix_(false)
#endif
  {
    Teuchos::Array<int> numEntriesPerRowToAlloc(NumEntriesPerRowToAlloc.begin(), NumEntriesPerRowToAlloc.end()); // convert array of "size_t" to array of "int"
    mtx_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, toEpetra<EpetraGlobalOrdinal,Node>(rowMap), numEntriesPerRowToAlloc.getRawPtr(), toEpetra(pftype)));
  }

  template<class EpetraGlobalOrdinal, class Node>
  EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::EpetraCrsMatrixT(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : mtx_(Teuchos::rcp(new Epetra_CrsMatrix(Copy, toEpetra<EpetraGlobalOrdinal,Node>(rowMap), toEpetra<EpetraGlobalOrdinal,Node>(colMap), maxNumEntriesPerRow, toEpetra(pftype)))), isFillResumed_(false)
#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
      , isInitializedLocalMatrix_(false)
#endif
  { }

  template<class EpetraGlobalOrdinal, class Node>
  EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::EpetraCrsMatrixT(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : isFillResumed_(false)
#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
      , isInitializedLocalMatrix_(false)
#endif
  {
    Teuchos::Array<int> numEntriesPerRowToAlloc(NumEntriesPerRowToAlloc.begin(), NumEntriesPerRowToAlloc.end()); // convert array of "size_t" to array of "int"
    mtx_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, toEpetra<EpetraGlobalOrdinal,Node>(rowMap), toEpetra<EpetraGlobalOrdinal,Node>(colMap), numEntriesPerRowToAlloc.getRawPtr(), toEpetra(pftype)));
  }

  template<class EpetraGlobalOrdinal, class Node>
  EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::EpetraCrsMatrixT(const Teuchos::RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node > > &graph, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : mtx_(Teuchos::rcp(new Epetra_CrsMatrix(Copy, toEpetra<EpetraGlobalOrdinal,Node>(graph)))), isFillResumed_(false)
#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
      , isInitializedLocalMatrix_(false)
#endif
  { }

  template<class EpetraGlobalOrdinal, class Node>
  EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::EpetraCrsMatrixT(const EpetraCrsMatrixT& matrix)
    : mtx_(Teuchos::rcp(new Epetra_CrsMatrix(*(matrix.mtx_)))), isFillResumed_(false)
#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
      , isInitializedLocalMatrix_(false)
#endif
  { }

  template<class EpetraGlobalOrdinal, class Node>
  EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::EpetraCrsMatrixT(const Teuchos::RCP<const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& sourceMatrix,
                                   const Import<LocalOrdinal,GlobalOrdinal,Node> &importer,
                                   const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
                                   const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
                                   const Teuchos::RCP<Teuchos::ParameterList>& params):
    isFillResumed_(false)
#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
      , isInitializedLocalMatrix_(false)
#endif
  {
    XPETRA_DYNAMIC_CAST(const EpetraCrsMatrixT<GlobalOrdinal COMMA Node>, *sourceMatrix, tSourceMatrix, "Xpetra::EpetraCrsMatrixT constructor only accepts Xpetra::EpetraCrsMatrixT as an input argument.");
    XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal COMMA Node>, importer, tImporter, "Xpetra::EpetraCrsMatrixT constructor only accepts Xpetra::EpetraImportT as an input argument.");

    const Epetra_Map* myDomainMap = (domainMap!=Teuchos::null)? &toEpetra<EpetraGlobalOrdinal,Node>(domainMap): 0;
    const Epetra_Map* myRangeMap  = (rangeMap !=Teuchos::null)? &toEpetra<EpetraGlobalOrdinal,Node>(rangeMap) : 0;

    // Follows the Tpetra parameters
    bool restrictComm=false;
    if(!params.is_null()) restrictComm = params->get("Restrict Communicator",restrictComm);
    mtx_ = Teuchos::rcp(new Epetra_CrsMatrix(*tSourceMatrix.getEpetra_CrsMatrix(),*tImporter.getEpetra_Import(),myDomainMap,myRangeMap,restrictComm));
    if(restrictComm && mtx_->NumMyRows()==0)
      mtx_=Teuchos::null;
  }

  template<class EpetraGlobalOrdinal, class Node>
  EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::EpetraCrsMatrixT(const Teuchos::RCP<const Xpetra::CrsMatrix<Scalar,LocalOrdinal,EpetraGlobalOrdinal,Node> >& sourceMatrix,
                                   const Export<LocalOrdinal,EpetraGlobalOrdinal,Node> &exporter,
                                   const Teuchos::RCP<const Map<LocalOrdinal,EpetraGlobalOrdinal,Node> >& domainMap,
                                   const Teuchos::RCP<const Map<LocalOrdinal,EpetraGlobalOrdinal,Node> >& rangeMap,
                                   const Teuchos::RCP<Teuchos::ParameterList>& params):
    isFillResumed_(false)
#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
      , isInitializedLocalMatrix_(false)
#endif
  {
    XPETRA_DYNAMIC_CAST(const EpetraCrsMatrixT<GlobalOrdinal COMMA Node>, *sourceMatrix, tSourceMatrix, "Xpetra::EpetraCrsMatrixT constructor only accepts Xpetra::EpetraCrsMatrixT as an input argument.");
    XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal COMMA Node>, exporter, tExporter, "Xpetra::EpetraCrsMatrixT constructor only accepts Xpetra::EpetraExportT as an input argument.");

    const Epetra_Map* myDomainMap = (domainMap!=Teuchos::null)? &toEpetra<EpetraGlobalOrdinal,Node>(domainMap): 0;
    const Epetra_Map* myRangeMap  = (rangeMap !=Teuchos::null)? &toEpetra<EpetraGlobalOrdinal,Node>(rangeMap) : 0;

    // Follows the Tpetra parameters
    bool restrictComm=false;
    if(!params.is_null()) restrictComm = params->get("Restrict Communicator",restrictComm);

    mtx_ = Teuchos::rcp(new Epetra_CrsMatrix(*tSourceMatrix.getEpetra_CrsMatrix(),*tExporter.getEpetra_Export(),myDomainMap,myRangeMap,restrictComm));
  }

#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
  template<class EpetraGlobalOrdinal, class Node>
  EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::EpetraCrsMatrixT (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& colMap,
        const local_matrix_type& lclMatrix,
        const Teuchos::RCP<Teuchos::ParameterList>& params) {
      // local typedefs from local_matrix_type
      typedef typename local_matrix_type::size_type size_type;
      typedef typename local_matrix_type::value_type value_type;
      typedef typename local_matrix_type::ordinal_type ordinal_type;

      // The number of rows in the sparse matrix.
      ordinal_type lclNumRows = lclMatrix.numRows ();
      ordinal_type lclNumCols = lclMatrix.numCols ();  // do we need this?

      // plausibility checks
      TEUCHOS_TEST_FOR_EXCEPTION(lclNumRows != Teuchos::as<ordinal_type>(rowMap->getNodeNumElements()), Xpetra::Exceptions::RuntimeError, "Xpetra::EpetraCrsMatrixT: number of rows in local matrix and number of local entries in row map do not match!");
      TEUCHOS_TEST_FOR_EXCEPTION(lclNumCols != Teuchos::as<ordinal_type>(colMap->getNodeNumElements()), Xpetra::Exceptions::RuntimeError, "Xpetra::EpetraCrsMatrixT: number of columns in local matrix and number of local entries in column map do not match!");

      std::vector<GlobalOrdinal> domainMapGids;  // vector for collecting domain map GIDs

      Teuchos::ArrayRCP< size_t > NumEntriesPerRowToAlloc(lclNumRows);
      for (ordinal_type r = 0; r < lclNumRows; ++r) {
        // extract data from current row r
        Kokkos::SparseRowView<local_matrix_type,size_type> rowview = lclMatrix.template row<size_type>(r);
        NumEntriesPerRowToAlloc[r] = rowview.length;
      }

      // setup matrix
      isFillResumed_ = false;
      Teuchos::Array<int> numEntriesPerRowToAlloc(NumEntriesPerRowToAlloc.begin(), NumEntriesPerRowToAlloc.end()); // convert array of "size_t" to array of "int"
      mtx_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, toEpetra<EpetraGlobalOrdinal,Node>(rowMap), toEpetra<EpetraGlobalOrdinal,Node>(colMap), numEntriesPerRowToAlloc.getRawPtr(), toEpetra(DynamicProfile)));

      // loop over all rows and colums of local matrix and fill matrix
      for (ordinal_type r = 0; r < lclNumRows; ++r) {
        // extract data from current row r
        Kokkos::SparseRowView<local_matrix_type,size_type> rowview = lclMatrix.template row<size_type>(r);

        // arrays for current row data
        Teuchos::ArrayRCP<ordinal_type> indout(rowview.length,Teuchos::ScalarTraits<ordinal_type>::zero());
        Teuchos::ArrayRCP<value_type>   valout(rowview.length,Teuchos::ScalarTraits<value_type>::zero());

        for(ordinal_type c = 0; c < rowview.length; c++) {
          value_type   value  = rowview.value  (c);
          ordinal_type colidx = rowview.colidx (c);

          TEUCHOS_TEST_FOR_EXCEPTION(colMap->isNodeLocalElement(colidx) == false, Xpetra::Exceptions::RuntimeError, "Xpetra::EpetraCrsMatrixT: local matrix contains column elements which are not in the provided column map!");

          indout [c] = colidx;
          valout [c] = value;

          // collect GIDs for domain map
          GlobalOrdinal gcid = colMap->getGlobalElement(c);
          if(rowMap->isNodeGlobalElement(gcid)) domainMapGids.push_back(gcid);
        }
        insertLocalValues(r, indout.view(0,indout.size()), valout.view(0,valout.size()));
      }

      // sort entries in domainMapGids and remove duplicates
      const GlobalOrdinal INVALID = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
      std::sort(domainMapGids.begin(), domainMapGids.end());
      domainMapGids.erase(std::unique(domainMapGids.begin(), domainMapGids.end()), domainMapGids.end());
      Teuchos::ArrayView<GlobalOrdinal> domainMapGidsView(&domainMapGids[0], domainMapGids.size());
      Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > domainMap =
          Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(colMap->lib(), INVALID, domainMapGidsView, colMap->getIndexBase(), colMap->getComm());

      // call fill complete
      this->fillComplete(domainMap, rowMap, params);

      // AP (2015/10/22): Could probably be optimized using provided lclMatrix, but lets not worry about that
      isInitializedLocalMatrix_ = false;
    }
#endif

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::insertGlobalValues(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &cols, const ArrayView<const Scalar> &vals) {
    XPETRA_MONITOR("EpetraCrsMatrixT::insertGlobalValues");
    XPETRA_ERR_CHECK(mtx_->InsertGlobalValues(globalRow, vals.size(), vals.getRawPtr(), cols.getRawPtr()));
  }

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::insertLocalValues(LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &cols, const ArrayView<const Scalar> &vals) {
    XPETRA_MONITOR("EpetraCrsMatrixT::insertLocalValues");
    XPETRA_ERR_CHECK(mtx_->InsertMyValues(localRow, vals.size(), vals.getRawPtr(), cols.getRawPtr()));
  }

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::replaceGlobalValues(GlobalOrdinal globalRow, const ArrayView< const GlobalOrdinal > &indices, const ArrayView< const Scalar > &values) {
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

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::replaceLocalValues(LocalOrdinal localRow, const ArrayView< const LocalOrdinal > &indices, const ArrayView< const Scalar > &values) {
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

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::allocateAllValues(size_t numNonZeros, ArrayRCP<size_t>& rowptr, ArrayRCP<LocalOrdinal>& colind, ArrayRCP<Scalar>& values) {
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

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::setAllValues(const ArrayRCP<size_t>& rowptr, const ArrayRCP<LocalOrdinal>& colind, const ArrayRCP<Scalar>& values) {
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

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::getAllValues(ArrayRCP<const size_t>& rowptr, ArrayRCP<const LocalOrdinal>& colind, ArrayRCP<const Scalar>& values) const {
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

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::resumeFill(const RCP< ParameterList > &params) {
    XPETRA_MONITOR("EpetraCrsMatrixT::resumeFill");

    // According to Tpetra documentation, resumeFill() may be called repeatedly.
    isFillResumed_ = true;
  }

  template<class EpetraGlobalOrdinal, class Node>
  bool EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::isFillComplete() const { XPETRA_MONITOR("EpetraCrsMatrixT::isFillComplete"); if (isFillResumed_) return false; else return mtx_->Filled(); }

  template<class EpetraGlobalOrdinal, class Node>
  bool EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::isFillActive() const { XPETRA_MONITOR("EpetraCrsMatrixT::isFillActive"); return !isFillComplete(); }

  template<class EpetraGlobalOrdinal, class Node>
  bool EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::supportsRowViews() const { XPETRA_MONITOR("EpetraCrsMatrixT::supportsRowViews"); return true; }

  //TODO: throw same exception as Tpetra
  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::getLocalRowCopy(LocalOrdinal LocalRow, const ArrayView<LocalOrdinal> &Indices, const ArrayView<Scalar> &Values, size_t &NumEntries) const {
    XPETRA_MONITOR("EpetraCrsMatrixT::getLocalRowCopy");

    int numEntries = -1;
    XPETRA_ERR_CHECK(mtx_->ExtractMyRowCopy(LocalRow, Indices.size(), numEntries, Values.getRawPtr(), Indices.getRawPtr()));
    NumEntries = numEntries;
  }

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::getGlobalRowCopy(GlobalOrdinal GlobalRow, const ArrayView<GlobalOrdinal> &Indices, const ArrayView<Scalar> &Values, size_t &NumEntries) const {
    XPETRA_MONITOR("EpetraCrsMatrixT::getGlobalRowCopy");

    int numEntries = -1;
    XPETRA_ERR_CHECK(mtx_->ExtractGlobalRowCopy(GlobalRow, Indices.size(), numEntries, Values.getRawPtr(), Indices.getRawPtr()));
    NumEntries = numEntries;
  }

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &indices, ArrayView<const Scalar> &values) const {
    XPETRA_MONITOR("EpetraCrsMatrixT::getGlobalRowView");

    int      numEntries;
    double * eValues;
    GlobalOrdinal * eIndices;

    XPETRA_ERR_CHECK(mtx_->ExtractGlobalRowView(GlobalRow, numEntries, eValues, eIndices));
    if (numEntries == 0) { eValues = NULL; eIndices = NULL; } // Cf. TEUCHOS_TEST_FOR_EXCEPT( p == 0 && size_in != 0 ) in Teuchos ArrayView constructor.

    indices = ArrayView<const GlobalOrdinal>(eIndices, numEntries);
    values  = ArrayView<const double>(eValues, numEntries);
  }

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal> &indices, ArrayView<const Scalar> &values) const {
    XPETRA_MONITOR("EpetraCrsMatrixT::getLocalRowView");

    int      numEntries;
    double * eValues;
    int    * eIndices;

    XPETRA_ERR_CHECK(mtx_->ExtractMyRowView(LocalRow, numEntries, eValues, eIndices));
    if (numEntries == 0) { eValues = NULL; eIndices = NULL; } // Cf. TEUCHOS_TEST_FOR_EXCEPT( p == 0 && size_in != 0 ) in Teuchos ArrayView constructor.

    indices = ArrayView<const int>(eIndices, numEntries);
    values  = ArrayView<const double>(eValues, numEntries);
  }

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::apply(const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &X, MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta) const {
    XPETRA_MONITOR("EpetraCrsMatrixT::apply");

    //TEUCHOS_TEST_FOR_EXCEPTION((alpha != 1) || (beta != 0), Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraCrsMatrixT.multiply() only accept alpha==1 and beta==0");

    XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT<GlobalOrdinal COMMA Node>, X, eX, "Xpetra::EpetraCrsMatrixT->apply() only accept Xpetra::EpetraMultiVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(      EpetraMultiVectorT<GlobalOrdinal COMMA Node>, Y, eY, "Xpetra::EpetraCrsMatrixT->apply() only accept Xpetra::EpetraMultiVectorT as input arguments.");

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

  template<class EpetraGlobalOrdinal, class Node>
  std::string EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::description() const {
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

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
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
  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::doImport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &source,
                                 const Import<LocalOrdinal, GlobalOrdinal, Node> &importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsMatrixT::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsMatrixT<GlobalOrdinal COMMA Node>, source, tSource, "Xpetra::EpetraCrsMatrixT::doImport only accept Xpetra::EpetraCrsMatrixT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal COMMA Node>, importer, tImporter, "Xpetra::EpetraCrsMatrixT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    RCP<const Epetra_CrsMatrix> v = tSource.getEpetra_CrsMatrix();
    int err = mtx_->Import(*v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::doExport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &dest,
                                 const Import<LocalOrdinal, GlobalOrdinal, Node>& importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsMatrixT::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsMatrixT<GlobalOrdinal COMMA Node>, dest, tDest, "Xpetra::EpetraCrsMatrixT::doImport only accept Xpetra::EpetraCrsMatrixT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal COMMA Node>, importer, tImporter, "Xpetra::EpetraCrsMatrixT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    RCP<const Epetra_CrsMatrix> v = tDest.getEpetra_CrsMatrix();
    int err = mtx_->Export(*v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::doImport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &source,
                                 const Export<LocalOrdinal, GlobalOrdinal, Node>& exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsMatrixT::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsMatrixT<GlobalOrdinal COMMA Node>, source, tSource, "Xpetra::EpetraCrsMatrixT::doImport only accept Xpetra::EpetraCrsMatrixT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal COMMA Node>, exporter, tExporter, "Xpetra::EpetraCrsMatrixT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    RCP<const Epetra_CrsMatrix> v = tSource.getEpetra_CrsMatrix();
    int err = mtx_->Import(*v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");

  }

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::doExport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &dest,
                                 const Export<LocalOrdinal, GlobalOrdinal, Node>& exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsMatrixT::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsMatrixT<GlobalOrdinal COMMA Node>, dest, tDest, "Xpetra::EpetraCrsMatrixT::doImport only accept Xpetra::EpetraCrsMatrixT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal COMMA Node>, exporter, tExporter, "Xpetra::EpetraCrsMatrixT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    RCP<const Epetra_CrsMatrix> v = tDest.getEpetra_CrsMatrix();
    int err = mtx_->Export(*v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

    template<class EpetraGlobalOrdinal, class Node>
    void EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::fillComplete(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &domainMap,
                                       const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rangeMap,
                                       const RCP< ParameterList > &params) {
      XPETRA_MONITOR("EpetraCrsMatrixT::fillComplete");

      // For Epetra matrices, resumeFill() is a fictive operation. There is no need for a fillComplete after some resumeFill() operations.
      if (isFillResumed_ == true) { isFillResumed_ = false; return; }

      bool doOptimizeStorage = true;
      if (params != null && params->get("Optimize Storage",true) == false) doOptimizeStorage = false;
      mtx_->FillComplete(toEpetra<EpetraGlobalOrdinal,Node>(domainMap), toEpetra<EpetraGlobalOrdinal,Node>(rangeMap), doOptimizeStorage);
    }

    template<class EpetraGlobalOrdinal, class Node>
    void EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::fillComplete(const RCP< ParameterList > &params) {
      XPETRA_MONITOR("EpetraCrsMatrixT::fillComplete");

      // For Epetra matrices, resumeFill() is a fictive operation. There is no need for a fillComplete after some resumeFill() operations.
      if (isFillResumed_ == true) { isFillResumed_ = false; return; }

      bool doOptimizeStorage = true;
      if (params != null && params->get("Optimize Storage",true) == false) doOptimizeStorage = false;
      mtx_->FillComplete(doOptimizeStorage);
    }


  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::replaceDomainMapAndImporter(const Teuchos::RCP< const  Map< LocalOrdinal, GlobalOrdinal, Node > >& newDomainMap, Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> >  & newImporter) {
      XPETRA_MONITOR("EpetraCrsMatrixT::replaceDomainMapAndImporter");
      XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal COMMA Node>, *newImporter, eImporter, "Xpetra::EpetraCrsMatrixT::replaceDomainMapAndImporter only accepts Xpetra::EpetraImportT.");

      const RCP<const Epetra_Import> & myImport = eImporter.getEpetra_Import();
      int rv=0;
      if(myImport==Teuchos::null)
        rv=mtx_->ReplaceDomainMapAndImporter( toEpetra<EpetraGlobalOrdinal,Node>(newDomainMap),0);
      else
        rv=mtx_->ReplaceDomainMapAndImporter( toEpetra<EpetraGlobalOrdinal,Node>(newDomainMap),&*myImport);
      TEUCHOS_TEST_FOR_EXCEPTION(rv != 0, std::runtime_error, "Xpetra::EpetraCrsMatrixT::replaceDomainMapAndImporter FAILED!");
  }


  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsMatrixT<EpetraGlobalOrdinal, Node>::expertStaticFillComplete(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & domainMap,
                                                 const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & rangeMap,
                                                 const RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > &importer,
                                                 const RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > &exporter,
                                                 const RCP<ParameterList> & params) {
    XPETRA_MONITOR("EpetraCrsMatrixT::expertStaticFillComplete");
    int rv=0;
    const Epetra_Import * myimport =0;
    const Epetra_Export * myexport =0;

    if(!importer.is_null()) {
      XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal COMMA Node>, *importer, eImporter, "Xpetra::EpetraCrsMatrixT::expertStaticFillComplete only accepts Xpetra::EpetraImportT.");
      myimport = eImporter.getEpetra_Import().getRawPtr();
    }
    if(!exporter.is_null()) {
      XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal COMMA Node>, *exporter, eExporter, "Xpetra::EpetraCrsMatrixT::expertStaticFillComplete only accepts Xpetra::EpetraImportT.");
      myexport = eExporter.getEpetra_Export().getRawPtr();
    }

    rv=mtx_->ExpertStaticFillComplete(toEpetra<EpetraGlobalOrdinal,Node>(domainMap), toEpetra<EpetraGlobalOrdinal,Node>(rangeMap), myimport, myexport);

    TEUCHOS_TEST_FOR_EXCEPTION(rv != 0, std::runtime_error, "Xpetra::EpetraCrsMatrixT::expertStaticFillComplete FAILED!");
  }

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
#ifdef HAVE_XPETRA_SERIAL
template class EpetraCrsMatrixT<int, Kokkos::Compat::KokkosSerialWrapperNode >;
#endif
#ifdef HAVE_XPETRA_PTHREAD
template class EpetraCrsMatrixT<int, Kokkos::Compat::KokkosThreadsWrapperNode>;
#endif
#ifdef HAVE_XPETRA_OPENMP
template class EpetraCrsMatrixT<int, Kokkos::Compat::KokkosOpenMPWrapperNode >;
#endif
#ifdef HAVE_XPETRA_CUDA
typedef Kokkos::Compat::KokkosCudaWrapperNode default_node_type;
template class EpetraCrsMatrixT<int, default_node_type >;
#endif
#else
// Tpetra is disabled and Kokkos not available: use dummy node type
typedef Kokkos::Compat::KokkosSerialWrapperNode default_node_type;
template class EpetraCrsMatrixT<int, default_node_type >;
#endif // HAVE_XPETRA_TPETRA
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
#ifdef HAVE_XPETRA_SERIAL
template class EpetraCrsMatrixT<long long, Kokkos::Compat::KokkosSerialWrapperNode >;
#endif
#ifdef HAVE_XPETRA_PTHREAD
template class EpetraCrsMatrixT<long long, Kokkos::Compat::KokkosThreadsWrapperNode>;
#endif
#ifdef HAVE_XPETRA_OPENMP
template class EpetraCrsMatrixT<long long, Kokkos::Compat::KokkosOpenMPWrapperNode >;
#endif
#ifdef HAVE_XPETRA_CUDA
typedef Kokkos::Compat::KokkosCudaWrapperNode default_node_type;
template class EpetraCrsMatrixT<long long, default_node_type >;
#endif
#else
// Tpetra is disabled and Kokkos not available: use dummy node type
typedef Kokkos::Compat::KokkosSerialWrapperNode default_node_type;
template class EpetraCrsMatrixT<long long, default_node_type >;
#endif // HAVE_XPETRA_TPETRA
#endif

}
