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
#ifndef XPETRA_TPETRACRSMATRIX_DEF_HPP
#define XPETRA_TPETRACRSMATRIX_DEF_HPP

#include "Xpetra_TpetraCrsMatrix_decl.hpp"

namespace Xpetra {

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TpetraCrsMatrix(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &params)
      : mtx_ (Teuchos::rcp (new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> (toTpetra(rowMap), maxNumEntriesPerRow, toTpetra(pftype), params))) {  }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TpetraCrsMatrix(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &params)
      : mtx_(Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> (toTpetra(rowMap), NumEntriesPerRowToAlloc(), toTpetra(pftype), params))) {  }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TpetraCrsMatrix(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &params)
      : mtx_(Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(toTpetra(rowMap), toTpetra(colMap), maxNumEntriesPerRow, toTpetra(pftype), params))) {  }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TpetraCrsMatrix(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &params)
      : mtx_(Teuchos::rcp(new Tpetra::CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node >(toTpetra(rowMap), toTpetra(colMap), NumEntriesPerRowToAlloc(), toTpetra(pftype), params))) {  }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TpetraCrsMatrix(const Teuchos::RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node > > &graph, const Teuchos::RCP< Teuchos::ParameterList > &params)
      : mtx_(Teuchos::rcp(new Tpetra::CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node >(toTpetra(graph), params))) {  }



    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TpetraCrsMatrix(const Teuchos::RCP<const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& sourceMatrix,
                    const Import<LocalOrdinal,GlobalOrdinal,Node> & importer,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
                    const Teuchos::RCP<Teuchos::ParameterList>& params)
    {
      typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> MyTpetraCrsMatrix;
      XPETRA_DYNAMIC_CAST(const TpetraCrsMatrixClass, *sourceMatrix, tSourceMatrix, "Xpetra::TpetraCrsMatrix constructor only accepts Xpetra::TpetraCrsMatrix as the input argument.");//TODO: remove and use toTpetra()
      RCP< const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > v = tSourceMatrix.getTpetra_CrsMatrix();

      RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > myDomainMap = domainMap!=Teuchos::null ? toTpetra(domainMap) : Teuchos::null;
      RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > myRangeMap  = rangeMap!=Teuchos::null  ? toTpetra(rangeMap)  : Teuchos::null;
      mtx_=Tpetra::importAndFillCompleteCrsMatrix<MyTpetraCrsMatrix>(tSourceMatrix.getTpetra_CrsMatrix(),toTpetra(importer),myDomainMap,myRangeMap,params);
      bool restrictComm=false;
      if(!params.is_null()) restrictComm = params->get("Restrict Communicator",restrictComm);
      if(restrictComm && mtx_->getRowMap().is_null()) mtx_=Teuchos::null;

    }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TpetraCrsMatrix(const Teuchos::RCP<const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& sourceMatrix,
                    const Export<LocalOrdinal,GlobalOrdinal,Node> & exporter,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
                    const Teuchos::RCP<Teuchos::ParameterList>& params)
    {
      typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> MyTpetraCrsMatrix;
      XPETRA_DYNAMIC_CAST(const TpetraCrsMatrixClass, *sourceMatrix, tSourceMatrix, "Xpetra::TpetraCrsMatrix constructor only accepts Xpetra::TpetraCrsMatrix as the input argument.");//TODO: remove and use toTpetra()
      RCP< const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > v = tSourceMatrix.getTpetra_CrsMatrix();

      RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > myDomainMap = domainMap!=Teuchos::null ? toTpetra(domainMap) : Teuchos::null;
      RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > myRangeMap  = rangeMap!=Teuchos::null  ? toTpetra(rangeMap)  : Teuchos::null;
      mtx_=Tpetra::exportAndFillCompleteCrsMatrix<MyTpetraCrsMatrix>(tSourceMatrix.getTpetra_CrsMatrix(),toTpetra(exporter),myDomainMap,myRangeMap,params);

    }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TpetraCrsMatrix(const Teuchos::RCP<const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& sourceMatrix,
                    const Import<LocalOrdinal,GlobalOrdinal,Node> & RowImporter,
                    const Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > DomainImporter,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
                    const Teuchos::RCP<Teuchos::ParameterList>& params)
    {
      typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> MyTpetraCrsMatrix;
      XPETRA_DYNAMIC_CAST(const TpetraCrsMatrixClass, *sourceMatrix, tSourceMatrix, "Xpetra::TpetraCrsMatrix constructor only accepts Xpetra::TpetraCrsMatrix as the input argument.");//TODO: remove and use toTpetra()
      RCP< const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > v = tSourceMatrix.getTpetra_CrsMatrix();

      RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > myDomainMap = domainMap!=Teuchos::null ? toTpetra(domainMap) : Teuchos::null;
      RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > myRangeMap  = rangeMap!=Teuchos::null  ? toTpetra(rangeMap)  : Teuchos::null;

      mtx_=Tpetra::importAndFillCompleteCrsMatrix<MyTpetraCrsMatrix>(tSourceMatrix.getTpetra_CrsMatrix(),toTpetra(RowImporter),toTpetra(*DomainImporter),myDomainMap,myRangeMap,params);
      bool restrictComm=false;
      if(!params.is_null()) restrictComm = params->get("Restrict Communicator",restrictComm);
      if(restrictComm && mtx_->getRowMap().is_null()) mtx_=Teuchos::null;
    }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TpetraCrsMatrix(const Teuchos::RCP<const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& sourceMatrix,
                    const Export<LocalOrdinal,GlobalOrdinal,Node> & RowExporter,
                    const Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > DomainExporter,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
                    const Teuchos::RCP<Teuchos::ParameterList>& params)
    {
      typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> MyTpetraCrsMatrix;
      XPETRA_DYNAMIC_CAST(const TpetraCrsMatrixClass, *sourceMatrix, tSourceMatrix, "Xpetra::TpetraCrsMatrix constructor only accepts Xpetra::TpetraCrsMatrix as the input argument.");//TODO: remove and use toTpetra()
      RCP< const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > v = tSourceMatrix.getTpetra_CrsMatrix();

      RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > myDomainMap = domainMap!=Teuchos::null ? toTpetra(domainMap) : Teuchos::null;
      RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > myRangeMap  = rangeMap!=Teuchos::null  ? toTpetra(rangeMap)  : Teuchos::null;

      mtx_=Tpetra::exportAndFillCompleteCrsMatrix<MyTpetraCrsMatrix>(tSourceMatrix.getTpetra_CrsMatrix(),toTpetra(RowExporter),toTpetra(*DomainExporter),myDomainMap,myRangeMap,params);
    }
 
///////////////////////////////////////////////////////////////////////////////////////


#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
#ifdef HAVE_XPETRA_TPETRA
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TpetraCrsMatrix (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& colMap,
        const local_matrix_type& lclMatrix,
        const Teuchos::RCP<Teuchos::ParameterList>& params)
    : mtx_(Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(toTpetra(rowMap), toTpetra(colMap), lclMatrix, params))) {  }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TpetraCrsMatrix (
        const local_matrix_type& lclMatrix,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& colMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
        const Teuchos::RCP<Teuchos::ParameterList>& params)
    : mtx_(Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(lclMatrix, toTpetra(rowMap), toTpetra(colMap), toTpetra(domainMap), toTpetra(rangeMap), params))) {  }
#endif
#endif

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::~TpetraCrsMatrix() {  }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::insertGlobalValues(GlobalOrdinal globalRow, const ArrayView< const GlobalOrdinal > &cols, const ArrayView< const Scalar > &vals) { XPETRA_MONITOR("TpetraCrsMatrix::insertGlobalValues"); mtx_->insertGlobalValues(globalRow, cols, vals); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::insertLocalValues(LocalOrdinal localRow, const ArrayView< const LocalOrdinal > &cols, const ArrayView< const Scalar > &vals) { XPETRA_MONITOR("TpetraCrsMatrix::insertLocalValues"); mtx_->insertLocalValues(localRow, cols, vals); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceGlobalValues(GlobalOrdinal globalRow, const ArrayView< const GlobalOrdinal > &cols, const ArrayView< const Scalar > &vals) { XPETRA_MONITOR("TpetraCrsMatrix::replaceGlobalValues"); mtx_->replaceGlobalValues(globalRow, cols, vals); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void
    TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceLocalValues (LocalOrdinal localRow,
                        const ArrayView<const LocalOrdinal> &cols,
                        const ArrayView<const Scalar> &vals)
    {
      XPETRA_MONITOR("TpetraCrsMatrix::replaceLocalValues");
      typedef typename ArrayView<const LocalOrdinal>::size_type size_type;
      const LocalOrdinal numValid =
        mtx_->replaceLocalValues (localRow, cols, vals);
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_type> (numValid) != cols.size (), std::runtime_error,
        "replaceLocalValues returned " << numValid << " != cols.size() = " <<
        cols.size () << ".");
    }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setAllToScalar(const Scalar &alpha) { XPETRA_MONITOR("TpetraCrsMatrix::setAllToScalar"); mtx_->setAllToScalar(alpha); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::scale(const Scalar &alpha) { XPETRA_MONITOR("TpetraCrsMatrix::scale"); mtx_->scale(alpha); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::allocateAllValues(size_t numNonZeros,ArrayRCP<size_t> & rowptr, ArrayRCP<LocalOrdinal> & colind, ArrayRCP<Scalar> & values)
    { XPETRA_MONITOR("TpetraCrsMatrix::allocateAllValues"); rowptr.resize(getNodeNumRows()+1); colind.resize(numNonZeros); values.resize(numNonZeros);}

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setAllValues(const ArrayRCP<size_t> & rowptr, const ArrayRCP<LocalOrdinal> & colind, const ArrayRCP<Scalar> & values)
    { XPETRA_MONITOR("TpetraCrsMatrix::setAllValues"); mtx_->setAllValues(rowptr,colind,values); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getAllValues(ArrayRCP<const size_t>& rowptr, ArrayRCP<const LocalOrdinal>& colind, ArrayRCP<const Scalar>& values) const
    { XPETRA_MONITOR("TpetraCrsMatrix::getAllValues"); mtx_->getAllValues(rowptr,colind,values); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    bool TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::haveGlobalConstants() const
    { return mtx_->haveGlobalConstants();}

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::resumeFill(const RCP< ParameterList > &params) { XPETRA_MONITOR("TpetraCrsMatrix::resumeFill"); mtx_->resumeFill(params); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::fillComplete(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &domainMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rangeMap, const RCP< ParameterList > &params) { XPETRA_MONITOR("TpetraCrsMatrix::fillComplete"); mtx_->fillComplete(toTpetra(domainMap), toTpetra(rangeMap), params); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::fillComplete(const RCP< ParameterList > &params) { XPETRA_MONITOR("TpetraCrsMatrix::fillComplete"); mtx_->fillComplete(params); }


    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceDomainMapAndImporter(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >& newDomainMap, Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> >  & newImporter) {
      XPETRA_MONITOR("TpetraCrsMatrix::replaceDomainMapAndImporter");
      XPETRA_DYNAMIC_CAST( const TpetraImportClass , *newImporter, tImporter, "Xpetra::TpetraCrsMatrix::replaceDomainMapAndImporter only accepts Xpetra::TpetraImport.");
      RCP<const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > myImport = tImporter.getTpetra_Import();
            mtx_->replaceDomainMapAndImporter( toTpetra(newDomainMap),myImport);
    }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::expertStaticFillComplete(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & domainMap,
                                  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & rangeMap,
                                  const RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > &importer,
                                  const RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > &exporter,
                                  const RCP<ParameterList> &params) {
      XPETRA_MONITOR("TpetraCrsMatrix::expertStaticFillComplete");
      RCP<const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > myImport;
      RCP<const Tpetra::Export<LocalOrdinal,GlobalOrdinal,Node> > myExport;

      if(importer!=Teuchos::null) {
        XPETRA_DYNAMIC_CAST( const TpetraImportClass , *importer, tImporter, "Xpetra::TpetraCrsMatrix::expertStaticFillComplete only accepts Xpetra::TpetraImport.");
        myImport = tImporter.getTpetra_Import();
      }
      if(exporter!=Teuchos::null) {
        XPETRA_DYNAMIC_CAST( const TpetraExportClass , *exporter, tExporter, "Xpetra::TpetraCrsMatrix::expertStaticFillComplete only accepts Xpetra::TpetraExport.");
        myExport = tExporter.getTpetra_Export();
      }

      mtx_->expertStaticFillComplete(toTpetra(domainMap),toTpetra(rangeMap),myImport,myExport,params);
    }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getRowMap() const { XPETRA_MONITOR("TpetraCrsMatrix::getRowMap"); return toXpetra(mtx_->getRowMap()); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getColMap() const { XPETRA_MONITOR("TpetraCrsMatrix::getColMap"); return toXpetra(mtx_->getColMap()); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node > > TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getCrsGraph() const { XPETRA_MONITOR("TpetraCrsMatrix::getCrsGraph"); return toXpetra(mtx_->getCrsGraph()); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    global_size_t TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getGlobalNumRows() const { XPETRA_MONITOR("TpetraCrsMatrix::getGlobalNumRows"); return mtx_->getGlobalNumRows(); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    global_size_t TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getGlobalNumCols() const { XPETRA_MONITOR("TpetraCrsMatrix::getGlobalNumCols"); return mtx_->getGlobalNumCols(); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    size_t TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getNodeNumRows() const { XPETRA_MONITOR("TpetraCrsMatrix::getNodeNumRows"); return mtx_->getNodeNumRows(); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    size_t TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getNodeNumCols() const { XPETRA_MONITOR("TpetraCrsMatrix::getNodeNumCols"); return mtx_->getNodeNumCols(); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    global_size_t TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getGlobalNumEntries() const { XPETRA_MONITOR("TpetraCrsMatrix::getGlobalNumEntries"); return mtx_->getGlobalNumEntries(); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    size_t TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getNodeNumEntries() const { XPETRA_MONITOR("TpetraCrsMatrix::getNodeNumEntries"); return mtx_->getNodeNumEntries(); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    size_t TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getNumEntriesInLocalRow(LocalOrdinal localRow) const { XPETRA_MONITOR("TpetraCrsMatrix::getNumEntriesInLocalRow"); return mtx_->getNumEntriesInLocalRow(localRow); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    size_t TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const { XPETRA_MONITOR("TpetraCrsMatrix::getNumEntriesInGlobalRow"); return mtx_->getNumEntriesInGlobalRow(globalRow); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    size_t TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getGlobalMaxNumRowEntries() const { XPETRA_MONITOR("TpetraCrsMatrix::getGlobalMaxNumRowEntries"); return mtx_->getGlobalMaxNumRowEntries(); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    size_t TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getNodeMaxNumRowEntries() const { XPETRA_MONITOR("TpetraCrsMatrix::getNodeMaxNumRowEntries"); return mtx_->getNodeMaxNumRowEntries(); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    bool TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isLocallyIndexed() const { XPETRA_MONITOR("TpetraCrsMatrix::isLocallyIndexed"); return mtx_->isLocallyIndexed(); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    bool TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isGloballyIndexed() const { XPETRA_MONITOR("TpetraCrsMatrix::isGloballyIndexed"); return mtx_->isGloballyIndexed(); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    bool TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isFillComplete() const { XPETRA_MONITOR("TpetraCrsMatrix::isFillComplete"); return mtx_->isFillComplete(); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    bool TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isFillActive() const { XPETRA_MONITOR("TpetraCrsMatrix::isFillActive"); return mtx_->isFillActive(); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    typename ScalarTraits< Scalar >::magnitudeType TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getFrobeniusNorm() const { XPETRA_MONITOR("TpetraCrsMatrix::getFrobeniusNorm"); return mtx_->getFrobeniusNorm(); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    bool TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::supportsRowViews() const { XPETRA_MONITOR("TpetraCrsMatrix::supportsRowViews"); return mtx_->supportsRowViews(); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalRowCopy(LocalOrdinal LocalRow, const ArrayView< LocalOrdinal > &Indices, const ArrayView< Scalar > &Values, size_t &NumEntries) const { XPETRA_MONITOR("TpetraCrsMatrix::getLocalRowCopy"); mtx_->getLocalRowCopy(LocalRow, Indices, Values, NumEntries); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView< const GlobalOrdinal > &indices, ArrayView< const Scalar > &values) const { XPETRA_MONITOR("TpetraCrsMatrix::getGlobalRowView"); mtx_->getGlobalRowView(GlobalRow, indices, values); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getGlobalRowCopy(GlobalOrdinal GlobalRow, const ArrayView< GlobalOrdinal > &indices, const ArrayView< Scalar > &values, size_t &numEntries) const { XPETRA_MONITOR("TpetraCrsMatrix::getGlobalRowCopy"); mtx_->getGlobalRowCopy(GlobalRow, indices, values, numEntries); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalRowView(LocalOrdinal LocalRow, ArrayView< const LocalOrdinal > &indices, ArrayView< const Scalar > &values) const { XPETRA_MONITOR("TpetraCrsMatrix::getLocalRowView"); mtx_->getLocalRowView(LocalRow, indices, values); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::apply(const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &X, MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta) const { XPETRA_MONITOR("TpetraCrsMatrix::apply"); mtx_->apply(toTpetra(X), toTpetra(Y), mode, alpha, beta); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getDomainMap() const { XPETRA_MONITOR("TpetraCrsMatrix::getDomainMap"); return toXpetra(mtx_->getDomainMap()); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getRangeMap() const { XPETRA_MONITOR("TpetraCrsMatrix::getRangeMap"); return toXpetra(mtx_->getRangeMap()); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    std::string TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::description() const { XPETRA_MONITOR("TpetraCrsMatrix::description"); return mtx_->description(); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const { XPETRA_MONITOR("TpetraCrsMatrix::describe"); mtx_->describe(out, verbLevel); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setObjectLabel( const std::string &objectLabel ) {
      XPETRA_MONITOR("TpetraCrsMatrix::setObjectLabel");
      Teuchos::LabeledObject::setObjectLabel(objectLabel);
      mtx_->setObjectLabel(objectLabel);
    }



    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TpetraCrsMatrix(const TpetraCrsMatrix& matrix)
      :  mtx_(Teuchos::rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(*(matrix.mtx_),Teuchos::Copy))) {}
      
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalDiagCopy(Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &diag) const {
      XPETRA_MONITOR("TpetraCrsMatrix::getLocalDiagCopy");
      XPETRA_DYNAMIC_CAST(TpetraVectorClass, diag, tDiag, "Xpetra::TpetraCrsMatrix.getLocalDiagCopy() only accept Xpetra::TpetraVector as input arguments.");
      mtx_->getLocalDiagCopy(*tDiag.getTpetra_Vector());
    }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalDiagOffsets(Teuchos::ArrayRCP<size_t> &offsets) const {
      XPETRA_MONITOR("TpetraCrsMatrix::getLocalDiagOffsets");
      mtx_->getLocalDiagOffsets(offsets);
    }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalDiagCopy(Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &diag, const Teuchos::ArrayView<const size_t> &offsets) const {
      XPETRA_MONITOR("TpetraCrsMatrix::getLocalDiagCopy");
      mtx_->getLocalDiagCopy(*(toTpetra(diag)), offsets);
    }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceDiag(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &diag) {
      XPETRA_MONITOR("TpetraCrsMatrix::replaceDiag");
      Tpetra::replaceDiagonalCrsMatrix(*mtx_, *(toTpetra(diag)));
    }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::leftScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {
      XPETRA_MONITOR("TpetraCrsMatrix::leftScale");
      mtx_->leftScale(*(toTpetra(x)));
    }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::rightScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {
      XPETRA_MONITOR("TpetraCrsMatrix::rightScale");
      mtx_->rightScale(*(toTpetra(x)));
    }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getMap() const { XPETRA_MONITOR("TpetraCrsMatrix::getMap"); return rcp( new TpetraMap< LocalOrdinal, GlobalOrdinal, Node >(mtx_->getMap()) ); }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::doImport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &source,
                  const Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM) {
      XPETRA_MONITOR("TpetraCrsMatrix::doImport");

      XPETRA_DYNAMIC_CAST(const TpetraCrsMatrixClass, source, tSource, "Xpetra::TpetraCrsMatrix::doImport only accept Xpetra::TpetraCrsMatrix as input arguments.");//TODO: remove and use toTpetra()
      RCP< const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > v = tSource.getTpetra_CrsMatrix();
      //mtx_->doImport(toTpetraCrsMatrix(source), *tImporter.getTpetra_Import(), toTpetra(CM));
      mtx_->doImport(*v, toTpetra(importer), toTpetra(CM));
    }

    //! Export.
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::doExport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &dest,
                  const Import< LocalOrdinal, GlobalOrdinal, Node >& importer, CombineMode CM) {
      XPETRA_MONITOR("TpetraCrsMatrix::doExport");

      XPETRA_DYNAMIC_CAST(const TpetraCrsMatrixClass, dest, tDest, "Xpetra::TpetraCrsMatrix::doImport only accept Xpetra::TpetraCrsMatrix as input arguments.");//TODO: remove and use toTpetra()
      RCP< const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > v = tDest.getTpetra_CrsMatrix();
      mtx_->doExport(*v, toTpetra(importer), toTpetra(CM));

    }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::doImport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &source,
                  const Export< LocalOrdinal, GlobalOrdinal, Node >& exporter, CombineMode CM) {
      XPETRA_MONITOR("TpetraCrsMatrix::doImport");

      XPETRA_DYNAMIC_CAST(const TpetraCrsMatrixClass, source, tSource, "Xpetra::TpetraCrsMatrix::doImport only accept Xpetra::TpetraCrsMatrix as input arguments.");//TODO: remove and use toTpetra()
      RCP< const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > v = tSource.getTpetra_CrsMatrix();
      mtx_->doImport(*v, toTpetra(exporter), toTpetra(CM));

    }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::doExport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &dest,
                  const Export< LocalOrdinal, GlobalOrdinal, Node >& exporter, CombineMode CM) {
      XPETRA_MONITOR("TpetraCrsMatrix::doExport");

      XPETRA_DYNAMIC_CAST(const TpetraCrsMatrixClass, dest, tDest, "Xpetra::TpetraCrsMatrix::doImport only accept Xpetra::TpetraCrsMatrix as input arguments.");//TODO: remove and use toTpetra()
      RCP< const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > v = tDest.getTpetra_CrsMatrix();
      mtx_->doExport(*v, toTpetra(exporter), toTpetra(CM));

    }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::removeEmptyProcessesInPlace (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& newMap) {
      XPETRA_MONITOR("TpetraCrsMatrix::removeEmptyProcessesInPlace");
      mtx_->removeEmptyProcessesInPlace(toTpetra(newMap));
    }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    bool TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::hasMatrix() const {
      return ! mtx_.is_null ();
    }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TpetraCrsMatrix(const Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &mtx) : mtx_(mtx) {  }

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getTpetra_CrsMatrix() const { return mtx_; }

    //! Get the underlying Tpetra matrix
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getTpetra_CrsMatrixNonConst() const { return mtx_; } //TODO: remove

////////////////////////////////////////////
////////////////////////////////////////////
// End of TpetrCrsMatrix class definition //
////////////////////////////////////////////
////////////////////////////////////////////

} // Xpetra namespace

#endif // XPETRA_TPETRACRSMATRIX_DEF_HPP
