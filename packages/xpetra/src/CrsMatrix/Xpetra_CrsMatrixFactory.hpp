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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef XPETRA_CRSMATRIXFACTORY_HPP
#define XPETRA_CRSMATRIXFACTORY_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_CrsMatrix.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraCrsMatrix.hpp"
#endif

#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraCrsMatrix.hpp"
#endif

#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

  template <class Scalar, class LocalOrdinal  = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps>
  class CrsMatrixFactory {

  private:
    //! Private constructor. This is a static class.
    CrsMatrixFactory() {}

  public:

    //! Constructor specifying fixed number of entries for each row.
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rowMap, size_t maxNumEntriesPerRow, Xpetra::ProfileType pftype = Xpetra::DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(rowMap, maxNumEntriesPerRow, pftype, plist) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(rowMap->lib());
      XPETRA_FACTORY_END;
    }

    //! Constructor specifying (possibly different) number of entries in each row.
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype = Xpetra::DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {

#ifdef HAVE_XPETRA_TPETRA
      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(rowMap, NumEntriesPerRowToAlloc, pftype, plist) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(rowMap->lib());
      XPETRA_FACTORY_END;
    }

    //! Constructor specifying column Map and fixed number of entries for each row.
    RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(rowMap, colMap, maxNumEntriesPerRow, pftype, plist) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(rowMap->lib());
      XPETRA_FACTORY_END;
    }

    //! Constructor specifying column Map and number of entries in each row.
    RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(rowMap, colMap, NumEntriesPerRowToAlloc, pftype, plist) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(rowMap->lib());
      XPETRA_FACTORY_END;
    }

    //! Constructor specifying a previously constructed graph.
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Build(const Teuchos::RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node, LocalMatOps > > &graph, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (graph->getRowMap()->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(graph, plist) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(graph->getRowMap()->lib());
      XPETRA_FACTORY_END;
    }

  };

  template <>
  class CrsMatrixFactory<double, int, int> {

    typedef double Scalar;
    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;
    typedef Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps LocalMatOps;

  private:
    //! Private constructor. This is a static class.
    CrsMatrixFactory() {}

  public:

    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rowMap, size_t maxNumEntriesPerRow, Xpetra::ProfileType pftype = Xpetra::DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(rowMap, maxNumEntriesPerRow, pftype, plist) );
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrix(rowMap, maxNumEntriesPerRow, pftype, plist) );
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype = Xpetra::DynamicProfile,  const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(rowMap, NumEntriesPerRowToAlloc, pftype, plist) );
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrix(rowMap, NumEntriesPerRowToAlloc, pftype, plist) );
#endif

      XPETRA_FACTORY_END;
    }

    //! Constructor specifying column Map and fixed number of entries for each row.
    RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(rowMap, colMap, maxNumEntriesPerRow, pftype, plist) );
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrix(rowMap, colMap, maxNumEntriesPerRow, pftype, plist) );
#endif

      XPETRA_FACTORY_END;
    }

    //! Constructor specifying column Map and number of entries in each row.
    RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(rowMap, colMap, NumEntriesPerRowToAlloc, pftype, plist) );
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrix(rowMap, colMap, NumEntriesPerRowToAlloc, pftype, plist) );
#endif

      XPETRA_FACTORY_END;
    }

    //! Constructor specifying a previously constructed graph.
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Build(const Teuchos::RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node, LocalMatOps > > &graph, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (graph->getRowMap()->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(graph, plist) );
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (graph->getRowMap()->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrix(graph, plist) );
#endif

      XPETRA_FACTORY_END;
    }
  };

}

#define XPETRA_CRSMATRIXFACTORY_SHORT
#endif
