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
/*
 * Xpetra_MapExtractorFactory.hpp
 *
 *  Created on: Aug 22, 2011
 *      Author: wiesner
 */

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_MAPEXTRACTORFACTORY_HPP_
#define XPETRA_MAPEXTRACTORFACTORY_HPP_

#include <Kokkos_DefaultNode.hpp>

#include <Xpetra_MapExtractor.hpp>

#ifdef HAVE_XPETRA_TPETRA
#include <Xpetra_TpetraMapExtractor.hpp>
#endif
#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraMapExtractor.hpp>
#endif

namespace Xpetra
{
  // factory class
  template <class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class MapExtractorFactory {
    typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractorClass;
    typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> MapClass;
#ifdef HAVE_XPETRA_TPETRA
    typedef TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMapClass;
    typedef TpetraMapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraMapExtractorClass;
#endif
  private:
    //! Private construtor, since this is a static class
    MapExtractorFactory() {}

  public:
    //! Constructor specifying the maps (and indirectly the used LINALG library (Tpetra vs Epetra))
    static Teuchos::RCP<MapExtractorClass> Build(const Teuchos::RCP<const MapClass>& fullmap, const std::vector<Teuchos::RCP<const MapClass> >& maps)
    {
#ifdef HAVE_XPETRA_TPETRA
      const Teuchos::RCP<const TpetraMapClass> &tMap = Teuchos::rcp_dynamic_cast<const TpetraMapClass>(fullmap);
      if (tMap != null)
      {
        Teuchos::RCP<TpetraMapExtractorClass> tMapExtractor = Teuchos::rcp(new TpetraMapExtractorClass(tMap,maps));
        return Teuchos::rcp_dynamic_cast<MapExtractorClass>(tMapExtractor);
      }
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(fullmap->lib());

      XPETRA_FACTORY_END;
    }

    //! Constructor for Tpetra
#ifdef HAVE_XPETRA_TPETRA
    static Teuchos::RCP<TpetraMapExtractorClass> Build(const Teuchos::RCP<const TpetraMapClass>& fullmap, const std::vector<Teuchos::RCP<const TpetraMapClass> >& maps)
    {
      const Teuchos::RCP<const TpetraMapClass> &tMap = Teuchos::rcp_dynamic_cast<const TpetraMapClass>(fullmap);
      if (tMap != null)
        //return Teuchos::rcp( new TpetraMultiVectorClass(map, NumVectors, zeroOut) );
        return Teuchos::rcp(new TpetraMapExtractorClass(tMap,maps));

      XPETRA_FACTORY_END;
    }
#endif

  };


  // factory class
  template <>
  class MapExtractorFactory<double,int,int> {
      typedef int LocalOrdinal;
      typedef int GlobalOrdinal;
      typedef double Scalar;
      typedef Kokkos::DefaultNode::DefaultNodeType Node;
    typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractorClass;
    typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> MapClass;
#ifdef HAVE_XPETRA_TPETRA
    typedef TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMapClass;
    typedef TpetraMapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraMapExtractorClass;
#endif
  private:
    //! Private construtor, since this is a static class
    MapExtractorFactory() {}

  public:
    //! Constructor specifying the maps (and indirectly the used LINALG library (Tpetra vs Epetra))
    static Teuchos::RCP<MapExtractorClass> Build(const Teuchos::RCP<const MapClass>& fullmap, const std::vector<Teuchos::RCP<const MapClass> >& maps)
    {
#ifdef HAVE_XPETRA_EPETRA
      const RCP<const EpetraMap> &eMap = Teuchos::rcp_dynamic_cast<const EpetraMap>(fullmap);
      if (eMap != null)
        return Teuchos::rcp(new Xpetra::EpetraMapExtractor(eMap,maps));
#endif
#ifdef HAVE_XPETRA_TPETRA
      const Teuchos::RCP<const TpetraMapClass> &tMap = Teuchos::rcp_dynamic_cast<const TpetraMapClass>(fullmap);
      if (tMap != null)
      {
        Teuchos::RCP<TpetraMapExtractorClass> tMapExtractor = Teuchos::rcp(new TpetraMapExtractorClass(tMap,maps));
        return Teuchos::rcp_dynamic_cast<MapExtractorClass>(tMapExtractor);
      }
#endif

      XPETRA_FACTORY_END;
    }

    //! Constructor for Tpetra
#ifdef HAVE_XPETRA_TPETRA
    static Teuchos::RCP<TpetraMapExtractorClass> Build(const Teuchos::RCP<const TpetraMapClass>& fullmap, const std::vector<Teuchos::RCP<const TpetraMapClass> >& maps)
    {
      const Teuchos::RCP<const TpetraMapClass> &tMap = Teuchos::rcp_dynamic_cast<const TpetraMapClass>(fullmap);
      if (tMap != null)
        //return Teuchos::rcp( new TpetraMultiVectorClass(map, NumVectors, zeroOut) );
        return Teuchos::rcp(new TpetraMapExtractorClass(tMap,maps));

      XPETRA_FACTORY_END;
    }
#endif

    //! Constructor for Epetra
#ifdef HAVE_XPETRA_EPETRA
    static Teuchos::RCP<Xpetra::EpetraMapExtractor> Build(const Teuchos::RCP<const Xpetra::EpetraMap>& fullmap, const std::vector<Teuchos::RCP<const Xpetra::EpetraMap> >& maps)
    {
      const RCP<const EpetraMap> &eMap = Teuchos::rcp_dynamic_cast<const EpetraMap>(fullmap);
      if (eMap != null)
        return Teuchos::rcp(new Xpetra::EpetraMapExtractor(eMap,maps));

      XPETRA_FACTORY_END;
    }
#endif
  };
}


#endif /* XPETRA_MAPEXTRACTORFACTORY_HPP_ */
