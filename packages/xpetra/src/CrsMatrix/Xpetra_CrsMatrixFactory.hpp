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
#ifndef XPETRA_CRSMATRIXFACTORY_HPP
#define XPETRA_CRSMATRIXFACTORY_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_CrsMatrix.hpp"

#include "Xpetra_TpetraCrsMatrix.hpp"
#include "Xpetra_TpetraBlockCrsMatrix.hpp"

#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraCrsMatrix.hpp"
#endif

#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class CrsMatrixFactory {
  private:
    //! Private constructor. This is a static class.
    CrsMatrixFactory() {}

  public:
    //! Constructor for empty matrix (intended use is an import/export target - can't insert entries directly)
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    Build (const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rowMap)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(rowMap->lib() == UseEpetra, std::logic_error,
          "Can't create Xpetra::EpetraCrsMatrix with these scalar/LO/GO types");
      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, 0) );

      XPETRA_FACTORY_END;
    }

    //! Constructor specifying fixed number of entries for each row.
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    Build (const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rowMap,
           size_t maxNumEntriesPerRow,
           const Teuchos::RCP<Teuchos::ParameterList>& plist = Teuchos::null)
    {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return Teuchos::rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, maxNumEntriesPerRow, plist) );

      XPETRA_FACTORY_ERROR_IF_EPETRA(rowMap->lib());
      XPETRA_FACTORY_END;
    }

    //! Constructor specifying (possibly different) number of entries in each row.
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    Build (const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >& rowMap,
           const ArrayRCP<const size_t>& NumEntriesPerRowToAlloc,
           const Teuchos::RCP<Teuchos::ParameterList>& plist = Teuchos::null)
    {
      if (rowMap->lib() == UseTpetra)
        return Teuchos::rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, NumEntriesPerRowToAlloc, plist) );

      XPETRA_FACTORY_ERROR_IF_EPETRA(rowMap->lib());
      XPETRA_FACTORY_END;
    }

    //! Constructor specifying column Map and fixed number of entries for each row.
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    Build (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& rowMap,
           const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& colMap,
           size_t maxNumEntriesPerRow,
           const Teuchos::RCP<Teuchos::ParameterList>& plist = Teuchos::null)
    {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, colMap, maxNumEntriesPerRow, plist) );

      XPETRA_FACTORY_ERROR_IF_EPETRA(rowMap->lib());
      XPETRA_FACTORY_END;
    }

    //! Constructor specifying column Map and number of entries in each row.
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, colMap, NumEntriesPerRowToAlloc, plist) );

      XPETRA_FACTORY_ERROR_IF_EPETRA(rowMap->lib());
      XPETRA_FACTORY_END;
    }

    //! Constructor specifying a previously constructed graph.
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(const Teuchos::RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node > > &graph, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (graph->getRowMap()->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(graph, plist) );

      XPETRA_FACTORY_ERROR_IF_EPETRA(graph->getRowMap()->lib());
      XPETRA_FACTORY_END;
    }


    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(
        const Teuchos::RCP< const CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > > &sourceMatrix,
        const Import<LocalOrdinal,GlobalOrdinal,Node> &importer,
        const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & domainMap = Teuchos::null,
        const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & rangeMap = Teuchos::null,
        const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (sourceMatrix->getRowMap()->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(sourceMatrix,importer,domainMap,rangeMap,params));

      XPETRA_FACTORY_ERROR_IF_EPETRA(sourceMatrix->getRowMap()->lib());
      XPETRA_FACTORY_END;
    }

    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(
        const Teuchos::RCP< const CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > > &sourceMatrix,
        const Export<LocalOrdinal,GlobalOrdinal,Node> &exporter,
        const RCP<Map<LocalOrdinal,GlobalOrdinal,Scalar> > & domainMap = Teuchos::null,
        const RCP<Map<LocalOrdinal,GlobalOrdinal,Scalar> > & rangeMap = Teuchos::null,
        const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (sourceMatrix->getRowMap()->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(sourceMatrix,exporter,domainMap,rangeMap,params));

      XPETRA_FACTORY_ERROR_IF_EPETRA(sourceMatrix->getRowMap()->lib());
      XPETRA_FACTORY_END;
    }

    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(
        const Teuchos::RCP< const CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > > &sourceMatrix,
        const Import<LocalOrdinal,GlobalOrdinal,Node> &RowImporter,
        const RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > DomainImporter,
        const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & domainMap,
        const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & rangeMap,
        const Teuchos::RCP<Teuchos::ParameterList>& params) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (sourceMatrix->getRowMap()->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(sourceMatrix,RowImporter,DomainImporter,domainMap,rangeMap,params));

      XPETRA_FACTORY_ERROR_IF_EPETRA(sourceMatrix->getRowMap()->lib());
      XPETRA_FACTORY_END;
    }

    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(
        const Teuchos::RCP< const CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > > &sourceMatrix,
        const Export<LocalOrdinal,GlobalOrdinal,Node> &RowExporter,
        const RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > DomainExporter,
        const RCP<Map<LocalOrdinal,GlobalOrdinal,Scalar> > & domainMap,
        const RCP<Map<LocalOrdinal,GlobalOrdinal,Scalar> > & rangeMap,
        const Teuchos::RCP<Teuchos::ParameterList>& params) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (sourceMatrix->getRowMap()->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(sourceMatrix,RowExporter,DomainExporter,domainMap,rangeMap,params));

      XPETRA_FACTORY_ERROR_IF_EPETRA(sourceMatrix->getRowMap()->lib());
      XPETRA_FACTORY_END;
    }

    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build (
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& colMap,
        const typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type& lclMatrix,
        const Teuchos::RCP<Teuchos::ParameterList>& params = null)  {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, colMap, lclMatrix, params));

      XPETRA_FACTORY_ERROR_IF_EPETRA(rowMap->lib());
      XPETRA_FACTORY_END;
    }

    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build (
        const typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type& lclMatrix,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& colMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap = Teuchos::null,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap = Teuchos::null,
        const Teuchos::RCP<Teuchos::ParameterList>& params = null)  {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(lclMatrix, rowMap, colMap, domainMap, rangeMap, params));

      XPETRA_FACTORY_ERROR_IF_EPETRA(rowMap->lib());
      XPETRA_FACTORY_END;
    }

    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build (
        const typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type& lclMatrix,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& colMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
        const Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node>>&  importer,
        const Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal,Node>>&  exporter,
        const Teuchos::RCP<Teuchos::ParameterList>& params = null)  {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(lclMatrix, rowMap, colMap, domainMap, rangeMap, importer, exporter, params));

      TEUCHOS_TEST_FOR_EXCEPTION(rowMap->lib() == UseEpetra, std::logic_error, "Epetra doesn't support this matrix constructor");

      XPETRA_FACTORY_END;
    }
    
    // Builds a BlockCrsMatrix
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > BuildBlock (
        const Teuchos::RCP<const Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >& blockGraph,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
        LocalOrdinal blockSize) {
  
      XPETRA_MONITOR("CrsMatrixFactory::BuildBlock");

      if (domainMap->lib() == UseTpetra) {
        return rcp(new Xpetra::TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(blockGraph,domainMap,rangeMap,blockSize) );
      }
      TEUCHOS_TEST_FOR_EXCEPTION(domainMap->lib() == UseEpetra, std::logic_error, "Epetra doesn't support this matrix constructor");

      XPETRA_FACTORY_END;
    }

  };

// we need the Epetra specialization only if Epetra is enabled
#if (defined(HAVE_XPETRA_EPETRA) && !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES))

  // Specializtion for SC=double, LO=int, GO=int and Node=EpetraNode
  // Used both for Epetra and Tpetra
  template <>
  class CrsMatrixFactory<double, int, int, EpetraNode> {
    typedef double Scalar;
    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef EpetraNode Node;

  private:
    //! Private constructor. This is a static class.
    CrsMatrixFactory() {}

  public:
    //! Constructor for empty matrix (intended use is an import/export target - can't insert entries directly)
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    Build (const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rowMap)
    {
      XPETRA_MONITOR("CrsMatrixFactory::Build");
      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, 0) );
      if(rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<int,Node>(rowMap));

      XPETRA_FACTORY_END;
    }

    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rowMap, size_t maxNumEntriesPerRow, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, maxNumEntriesPerRow, plist) );

      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<int,Node>(rowMap, maxNumEntriesPerRow, plist) );

      XPETRA_FACTORY_END;
    }

    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, NumEntriesPerRowToAlloc, plist) );

      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<int,Node>(rowMap, NumEntriesPerRowToAlloc, plist) );

      XPETRA_FACTORY_END;
    }

    //! Constructor specifying column Map and fixed number of entries for each row.
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, size_t maxNumEntriesPerRow, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, colMap, maxNumEntriesPerRow, plist) );

      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<int,Node>(rowMap, colMap, maxNumEntriesPerRow, plist) );

      XPETRA_FACTORY_END;
    }

    //! Constructor specifying column Map and number of entries in each row.
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, colMap, NumEntriesPerRowToAlloc, plist) );

      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<int,Node>(rowMap, colMap, NumEntriesPerRowToAlloc, plist) );

      XPETRA_FACTORY_END;
    }

    //! Constructor specifying a previously constructed graph.
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(const Teuchos::RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node > > &graph, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (graph->getRowMap()->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(graph, plist) );

      if (graph->getRowMap()->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<int,Node>(graph, plist) );

      XPETRA_FACTORY_END;
    }


    //! Constructor using FusedImport
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(
        const Teuchos::RCP< const CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > > &sourceMatrix,
        const Import<LocalOrdinal,GlobalOrdinal,Node> &importer,
        const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & domainMap = Teuchos::null,
        const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & rangeMap = Teuchos::null,
        const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (sourceMatrix->getRowMap()->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(sourceMatrix,importer,domainMap,rangeMap,params) );

      if (sourceMatrix->getRowMap()->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<int,Node>(sourceMatrix,importer,domainMap,rangeMap,params) );

      XPETRA_FACTORY_END;
    }

    //! Constructor using FusedExport
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(
        const Teuchos::RCP< const CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > > &sourceMatrix,
        const Export<LocalOrdinal,GlobalOrdinal,Node> &exporter,
        const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & domainMap = Teuchos::null,
        const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & rangeMap = Teuchos::null,
        const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (sourceMatrix->getRowMap()->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(sourceMatrix,exporter,domainMap,rangeMap,params) );

      if (sourceMatrix->getRowMap()->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<int,Node>(sourceMatrix,exporter,domainMap,rangeMap,params) );

      XPETRA_FACTORY_END;
    }

    //! Constructor using FusedImport
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(
        const Teuchos::RCP< const CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > > &sourceMatrix,
        const Import<LocalOrdinal,GlobalOrdinal,Node> & RowImporter,
        const RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > DomainImporter,
        const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & domainMap,
        const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & rangeMap,
        const Teuchos::RCP<Teuchos::ParameterList>& params) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (sourceMatrix->getRowMap()->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(sourceMatrix,RowImporter,DomainImporter,domainMap,rangeMap,params) );

      if (sourceMatrix->getRowMap()->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<int,Node>(sourceMatrix,RowImporter,DomainImporter,domainMap,rangeMap,params) );

      XPETRA_FACTORY_END;
    }

    //! Constructor using FusedExport
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(
        const Teuchos::RCP< const CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > > &sourceMatrix,
        const Export<LocalOrdinal,GlobalOrdinal,Node> &RowExporter,
        const RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > DomainExporter,
        const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & domainMap,
        const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & rangeMap,
        const Teuchos::RCP<Teuchos::ParameterList>& params) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (sourceMatrix->getRowMap()->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(sourceMatrix,RowExporter,DomainExporter,domainMap,rangeMap,params) );

      if (sourceMatrix->getRowMap()->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<int,Node>(sourceMatrix,RowExporter,DomainExporter,domainMap,rangeMap,params) );

      XPETRA_FACTORY_END;
    }


    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build (
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& colMap,
        const typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type& lclMatrix,
        const Teuchos::RCP<Teuchos::ParameterList>& params = null)  {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, colMap, lclMatrix, params));

      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<int,Node>(rowMap, colMap, lclMatrix, params) );

      XPETRA_FACTORY_END;
    }

    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build (
        const typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type& lclMatrix,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& colMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap = Teuchos::null,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap = Teuchos::null,
        const Teuchos::RCP<Teuchos::ParameterList>& params = null)  {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(lclMatrix, rowMap, colMap, domainMap, rangeMap, params));

      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<int,Node>(lclMatrix, rowMap, colMap, domainMap, rangeMap, params) );

      XPETRA_FACTORY_END;
    }

    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build (
        const typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type& lclMatrix,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& colMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
        const Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node>>&  importer,
        const Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal,Node>>&  exporter,
        const Teuchos::RCP<Teuchos::ParameterList>& params = null)  {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(lclMatrix, rowMap, colMap, domainMap, rangeMap, importer, exporter, params));

      TEUCHOS_TEST_FOR_EXCEPTION(rowMap->lib() == UseEpetra, std::logic_error, "Epetra doesn't support this matrix constructor");

      XPETRA_FACTORY_END;
    }

    //! Build a BlockCrsMatrix
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > BuildBlock (
        const Teuchos::RCP<const Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >& blockGraph,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
        LocalOrdinal blockSize) {
  
      XPETRA_MONITOR("CrsMatrixFactory::BuildBlock");
      if (domainMap->lib() == UseTpetra)
        return rcp(new Xpetra::TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(blockGraph,domainMap,rangeMap,blockSize) );
      TEUCHOS_TEST_FOR_EXCEPTION(domainMap->lib() == UseEpetra, std::logic_error, "Epetra doesn't support this matrix constructor");

      XPETRA_FACTORY_END;
    }

  };
#endif

// we need the Epetra specialization only if Epetra is enabled
#if (defined(HAVE_XPETRA_EPETRA) && !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES))

  template <>
  class CrsMatrixFactory<double, int, long long, EpetraNode> {
    typedef double Scalar;
    typedef int LocalOrdinal;
    typedef long long GlobalOrdinal;
    typedef EpetraNode Node;

  private:
    //! Private constructor. This is a static class.
    CrsMatrixFactory() {}

  public:
    //! Constructor for empty matrix (intended use is an import/export target - can't insert entries directly)
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    Build (const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rowMap)
    {
      XPETRA_MONITOR("CrsMatrixFactory::Build");
      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, 0) );
#ifdef HAVE_XPETRA_EPETRA
      if(rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<long long,Node>(rowMap,0));
#endif
      XPETRA_FACTORY_END;
    }

    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rowMap, size_t maxNumEntriesPerRow, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, maxNumEntriesPerRow, plist) );

      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<long long, Node>(rowMap, maxNumEntriesPerRow, plist) );

      XPETRA_FACTORY_END;
    }

    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, NumEntriesPerRowToAlloc, plist) );

      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<long long, Node>(rowMap, NumEntriesPerRowToAlloc, plist) );

      XPETRA_FACTORY_END;
    }

    //! Constructor specifying column Map and fixed number of entries for each row.
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, size_t maxNumEntriesPerRow, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, colMap, maxNumEntriesPerRow, plist) );

      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<long long, Node>(rowMap, colMap, maxNumEntriesPerRow, plist) );

      XPETRA_FACTORY_END;
    }

    //! Constructor specifying column Map and number of entries in each row.
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, colMap, NumEntriesPerRowToAlloc, plist) );

      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<long long, Node>(rowMap, colMap, NumEntriesPerRowToAlloc, plist) );

      XPETRA_FACTORY_END;
    }

    //! Constructor specifying a previously constructed graph.
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(const Teuchos::RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node > > &graph, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (graph->getRowMap()->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(graph, plist) );

      if (graph->getRowMap()->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<long long, Node>(graph, plist) );

      XPETRA_FACTORY_END;
    }


    //! Constructor using FusedImport
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(
        const Teuchos::RCP< const CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > > &sourceMatrix,
        const Import<LocalOrdinal,GlobalOrdinal,Node> &importer,
        const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & domainMap = Teuchos::null,
        const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & rangeMap = Teuchos::null,
        const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (sourceMatrix->getRowMap()->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(sourceMatrix,importer,domainMap,rangeMap,params) );

      if (sourceMatrix->getRowMap()->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<long long, Node>(sourceMatrix,importer,domainMap,rangeMap,params) );

      XPETRA_FACTORY_END;
    }

    //! Constructor using FusedExport
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(
        const Teuchos::RCP< const CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > > &sourceMatrix,
        const Export<LocalOrdinal,GlobalOrdinal,Node> &exporter,
        const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & domainMap = Teuchos::null,
        const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & rangeMap = Teuchos::null,
        const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (sourceMatrix->getRowMap()->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(sourceMatrix,exporter,domainMap,rangeMap,params) );

      if (sourceMatrix->getRowMap()->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<long long, Node>(sourceMatrix,exporter,domainMap,rangeMap,params) );

      XPETRA_FACTORY_END;
    }

    //! Constructor using FusedImport
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(
        const Teuchos::RCP< const CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > > &sourceMatrix,
        const Import<LocalOrdinal,GlobalOrdinal,Node> & RowImporter,
        const RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > DomainImporter,
        const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & domainMap,
        const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & rangeMap,
        const Teuchos::RCP<Teuchos::ParameterList>& params) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (sourceMatrix->getRowMap()->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(sourceMatrix,RowImporter,DomainImporter,domainMap,rangeMap,params) );

      if (sourceMatrix->getRowMap()->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<long long,Node>(sourceMatrix,RowImporter,DomainImporter,domainMap,rangeMap,params) );

      XPETRA_FACTORY_END;
    }

    //! Constructor using FusedExport
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(
        const Teuchos::RCP< const CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > > &sourceMatrix,
        const Export<LocalOrdinal,GlobalOrdinal,Node> &RowExporter,
        const RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > DomainExporter,
        const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & domainMap,
        const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & rangeMap,
        const Teuchos::RCP<Teuchos::ParameterList>& params) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (sourceMatrix->getRowMap()->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(sourceMatrix,RowExporter,DomainExporter,domainMap,rangeMap,params) );

      if (sourceMatrix->getRowMap()->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<long long,Node>(sourceMatrix,RowExporter,DomainExporter,domainMap,rangeMap,params) );

      XPETRA_FACTORY_END;
    }

    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build (
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& colMap,
        const typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type& lclMatrix,
        const Teuchos::RCP<Teuchos::ParameterList>& params = null)  {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, colMap, lclMatrix, params));

      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<long long, Node>(rowMap, colMap, lclMatrix, params) );

      XPETRA_FACTORY_END;
    }
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build (
        const typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type& lclMatrix,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& colMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap = Teuchos::null,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap = Teuchos::null,
        const Teuchos::RCP<Teuchos::ParameterList>& params = null)  {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(lclMatrix, rowMap, colMap, domainMap, rangeMap, params));

      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrixT<long long, Node>(lclMatrix, rowMap, colMap, domainMap, rangeMap, params) );

      XPETRA_FACTORY_END;
    }


    //! Build a BlockCrsMatrix
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > BuildBlock (
        const Teuchos::RCP<const Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >& blockGraph,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
        LocalOrdinal blockSize) {
  
      XPETRA_MONITOR("CrsMatrixFactory::BuildBlock");

      if (domainMap->lib() == UseTpetra) {
        return rcp(new Xpetra::TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(blockGraph,domainMap,rangemap,blockSize) );
      }
      TEUCHOS_TEST_FOR_EXCEPTION(domainMap->lib() == UseEpetra, std::logic_error, "Epetra doesn't support this matrix constructor");

      XPETRA_FACTORY_END;
    }



  };
#endif

}

#define XPETRA_CRSMATRIXFACTORY_SHORT
#endif
