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
#ifndef XPETRA_CRSGRAPHFACTORY_HPP
#define XPETRA_CRSGRAPHFACTORY_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_CrsGraph.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraCrsGraph.hpp"
#endif

#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraCrsGraph.hpp"
#endif

#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

  template <class LocalOrdinal/* = CrsGraph<>::local_ordinal_type*/,
            class GlobalOrdinal/* =
              typename CrsGraph<LocalOrdinal>::global_ordinal_type*/,
            class Node/* =
              typename CrsGraph<LocalOrdinal, GlobalOrdinal>::node_type*/>
  class CrsGraphFactory {
  private:
    //! Private constructor. This is a static class.
    CrsGraphFactory() {}

  public:
    //! Constructor specifying the number of non-zeros for all rows.
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, size_t NumVectors, ProfileType pftype=DynamicProfile) {
      XPETRA_MONITOR("CrsGraphFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node> (map, NumVectors, pftype) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(map->lib());
      XPETRA_FACTORY_END;
      TEUCHOS_UNREACHABLE_RETURN(null);
    }

    //! Constructor specifying column Map and number of entries in each row.
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsGraphFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(rowMap, colMap, NumEntriesPerRowToAlloc, pftype, plist) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(rowMap->lib());
      XPETRA_FACTORY_END;
      TEUCHOS_UNREACHABLE_RETURN(null);
    }
  };

// we need the Epetra specialization only if Epetra is enabled
#if (defined(HAVE_XPETRA_EPETRA) && !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES))

  template <>
  class CrsGraphFactory<int, int, EpetraNode> {

    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef EpetraNode Node;

  private:
    //! Private constructor. This is a static class.
    CrsGraphFactory() {}

  public:

    static RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, size_t NumVectors, ProfileType pftype=DynamicProfile) {
      XPETRA_MONITOR("CrsGraphFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node> (map, NumVectors, pftype) );
#endif

      if (map->lib() == UseEpetra)
        return rcp( new EpetraCrsGraphT<int, Node>(map, NumVectors, pftype) );

      XPETRA_FACTORY_END;
      TEUCHOS_UNREACHABLE_RETURN(null);
    }

    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsGraphFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(rowMap, colMap, NumEntriesPerRowToAlloc, pftype, plist) );
#endif

      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsGraphT<int, Node>(rowMap, colMap, NumEntriesPerRowToAlloc, pftype, plist) );

      XPETRA_FACTORY_END;
      TEUCHOS_UNREACHABLE_RETURN(null);
    }


  };
#endif

// we need the Epetra specialization only if Epetra is enabled
#if (defined(HAVE_XPETRA_EPETRA) && !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES))

  template <>
  class CrsGraphFactory<int, long long, EpetraNode> {

    typedef int LocalOrdinal;
    typedef long long GlobalOrdinal;
    typedef EpetraNode Node;

  private:
    //! Private constructor. This is a static class.
    CrsGraphFactory() {}

  public:

    static RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, size_t NumVectors, ProfileType pftype=DynamicProfile) {
      XPETRA_MONITOR("CrsGraphFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node> (map, NumVectors, pftype) );
#endif

      if (map->lib() == UseEpetra)
        return rcp( new EpetraCrsGraphT<long long, Node>(map, NumVectors, pftype) );

      XPETRA_FACTORY_END;
      TEUCHOS_UNREACHABLE_RETURN(null);
    }

    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsGraphFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(rowMap, colMap, NumEntriesPerRowToAlloc, pftype, plist) );
#endif

      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsGraphT<long long, Node>(rowMap, colMap, NumEntriesPerRowToAlloc, pftype, plist) );

      XPETRA_FACTORY_END;
      TEUCHOS_UNREACHABLE_RETURN(null);
    }

  };
#endif
}

#define XPETRA_CRSGRAPHFACTORY_SHORT
#endif
