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
#ifndef XPETRA_MULTIVECTORFACTORY_HPP
#define XPETRA_MULTIVECTORFACTORY_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_MultiVector.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraMultiVector.hpp"
#endif

#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraMultiVector.hpp"
#include "Xpetra_EpetraIntMultiVector.hpp"
#endif

#include "Xpetra_BlockedMap.hpp"
#include "Xpetra_BlockedMultiVector.hpp"

#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of BlockedMultiVector, needed to prevent circular inclusions
  //template<class SC, class LO, class GO, class N> class BlockedMultiVector;
#endif

  template <class Scalar = MultiVector<>::scalar_type,
            class LocalOrdinal =
              typename MultiVector<Scalar>::local_ordinal_type,
            class GlobalOrdinal =
              typename MultiVector<Scalar, LocalOrdinal>::global_ordinal_type,
            class Node =
              typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
  class MultiVectorFactory {
  private:
    //! Private constructor. This is a static class.
    MultiVectorFactory() {}

  public:

    //! Constructor specifying the number of non-zeros for all rows.
    static Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    Build (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map,
           size_t NumVectors,
           bool zeroOut=true)
    {
      XPETRA_MONITOR("MultiVectorFactory::Build");

      RCP<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node> > bmap = Teuchos::rcp_dynamic_cast<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node> >(map);
      if(!bmap.is_null()) {
        return rcp(new Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(bmap, NumVectors, zeroOut));
      }

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> (map, NumVectors, zeroOut) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(map->lib());
      XPETRA_FACTORY_END;
    }

    //! Set multi-vector values from array of pointers using Teuchos memory management classes. (copy).
    static Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &map,
          const Teuchos::ArrayView< const Teuchos::ArrayView< const Scalar > > &ArrayOfPtrs,
          size_t NumVectors) {
      XPETRA_MONITOR("MultiVectorFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> (map, ArrayOfPtrs, NumVectors) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(map->lib());
      XPETRA_FACTORY_END;
    }

  };


// we need the Epetra specialization only if Epetra is enabled
#if (defined(HAVE_XPETRA_EPETRA) && !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES))


  // Specialization for Scalar=double, LO=GO=int and Serial node
  // Used both for Epetra and Tpetra
  // For any other node definition the general default implementation is used which allows Tpetra only
  template <>
  class MultiVectorFactory<double, int, int, EpetraNode> {

    typedef double Scalar;
    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef EpetraNode Node;

  private:
    //! Private constructor. This is a static class.
    MultiVectorFactory() {}

  public:

    static RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map,
          size_t NumVectors,
          bool zeroOut=true) {
      XPETRA_MONITOR("MultiVectorFactory::Build");

      RCP<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node> > bmap = Teuchos::rcp_dynamic_cast<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node> >(map);
      if(!bmap.is_null()) {
        return rcp(new BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(bmap, NumVectors, zeroOut));
      }

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> (map, NumVectors, zeroOut) );
#endif

      if (map->lib() == UseEpetra)
        return rcp( new EpetraMultiVectorT<int,Node>(map, NumVectors, zeroOut) );

      XPETRA_FACTORY_END;
    }

    //! Set multi-vector values from array of pointers using Teuchos memory management classes. (copy).
    static Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &map,
          const Teuchos::ArrayView< const Teuchos::ArrayView< const Scalar > > &ArrayOfPtrs,
          size_t NumVectors) {
      XPETRA_MONITOR("MultiVectorFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> (map, ArrayOfPtrs, NumVectors) );
#endif

      if (map->lib() == UseEpetra)
        return rcp( new EpetraMultiVectorT<int,Node>(map, ArrayOfPtrs, NumVectors) );

      XPETRA_FACTORY_END;
    }

  };

  template <>
  class MultiVectorFactory<int, int, int, EpetraNode> {

    typedef int Scalar;
    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef EpetraNode Node;

  private:
    //! Private constructor. This is a static class.
    MultiVectorFactory() {}

  public:

    static RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map,
          size_t NumVectors,
          bool zeroOut=true) {
      XPETRA_MONITOR("MultiVectorFactory::Build");

      RCP<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node> > bmap = Teuchos::rcp_dynamic_cast<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node> >(map);
      if(!bmap.is_null()) {
        return rcp(new BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(bmap, NumVectors, zeroOut));
      }

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> (map, NumVectors, zeroOut) );
#endif

      if (map->lib() == UseEpetra)
        return rcp( new EpetraIntMultiVectorT<int,Node>(map, NumVectors, zeroOut) );

      XPETRA_FACTORY_END;
    }

    //! Set multi-vector values from array of pointers using Teuchos memory management classes. (copy).
    static Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &map,
          const Teuchos::ArrayView< const Teuchos::ArrayView< const Scalar > > &ArrayOfPtrs,
          size_t NumVectors) {
      XPETRA_MONITOR("MultiVectorFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> (map, ArrayOfPtrs, NumVectors) );
#endif

      if (map->lib() == UseEpetra)
        return rcp( new EpetraIntMultiVectorT<int,Node>(map, ArrayOfPtrs, NumVectors) );

      XPETRA_FACTORY_END;
    }

  };

// we need the Epetra specialization only if Epetra is enabled
#if (defined(HAVE_XPETRA_EPETRA) && !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES))

  template <>
  class MultiVectorFactory<double, int, long long, EpetraNode> {

    typedef double Scalar;
    typedef int LocalOrdinal;
    typedef long long GlobalOrdinal;
    typedef EpetraNode Node;

  private:
    //! Private constructor. This is a static class.
    MultiVectorFactory() {}

  public:

    static RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map,
          size_t NumVectors,
          bool zeroOut = true) {
      XPETRA_MONITOR("MultiVectorFactory::Build");

      RCP<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node> > bmap = Teuchos::rcp_dynamic_cast<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node> >(map);
      if(!bmap.is_null()) {
        return rcp(new BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(bmap, NumVectors, zeroOut));
      }

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> (map, NumVectors, zeroOut) );
#endif

      if (map->lib() == UseEpetra)
        return rcp( new EpetraMultiVectorT<long long,Node>(map, NumVectors, zeroOut) );

      XPETRA_FACTORY_END;
    }

    //! Set multi-vector values from array of pointers using Teuchos memory management classes. (copy).
    static Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &map,
          const Teuchos::ArrayView< const Teuchos::ArrayView< const Scalar > > &ArrayOfPtrs,
          size_t NumVectors) {
      XPETRA_MONITOR("MultiVectorFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> (map, ArrayOfPtrs, NumVectors) );
#endif

      if (map->lib() == UseEpetra)
        return rcp( new EpetraMultiVectorT<long long,Node>(map, ArrayOfPtrs, NumVectors) );

      XPETRA_FACTORY_END;
    }

  };

  template <>
  class MultiVectorFactory<int, int, long long, EpetraNode> {

    typedef int Scalar;
    typedef int LocalOrdinal;
    typedef long long GlobalOrdinal;
    typedef EpetraNode Node;

  private:
    //! Private constructor. This is a static class.
    MultiVectorFactory() {}

  public:

    static RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map,
          size_t NumVectors,
          bool zeroOut = true) {
      XPETRA_MONITOR("MultiVectorFactory::Build");

      RCP<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node> > bmap = Teuchos::rcp_dynamic_cast<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node> >(map);
      if(!bmap.is_null()) {
        return rcp(new BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(bmap, NumVectors, zeroOut));
      }

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> (map, NumVectors, zeroOut) );
#endif

      if (map->lib() == UseEpetra)
        return rcp( new EpetraIntMultiVectorT<long long,Node>(map, NumVectors, zeroOut) );

      XPETRA_FACTORY_END;
    }

    //! Set multi-vector values from array of pointers using Teuchos memory management classes. (copy).
    static Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &map,
          const Teuchos::ArrayView< const Teuchos::ArrayView< const Scalar > > &ArrayOfPtrs,
          size_t NumVectors) {
      XPETRA_MONITOR("MultiVectorFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> (map, ArrayOfPtrs, NumVectors) );
#endif

      if (map->lib() == UseEpetra)
        return rcp( new EpetraIntMultiVectorT<long long,Node>(map, ArrayOfPtrs, NumVectors) );

      XPETRA_FACTORY_END;
    }

  };
#endif // Epetra64
#endif // HAVE_XPETRA_EPETRA

}

#define XPETRA_MULTIVECTORFACTORY_SHORT
#endif
