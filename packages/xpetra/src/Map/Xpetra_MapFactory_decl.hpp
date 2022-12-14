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
#ifndef XPETRA_MAPFACTORY_DECL_HPP
#define XPETRA_MAPFACTORY_DECL_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_Map_decl.hpp"
#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

/// \class MapFactory
/// \brief Create an Xpetra::Map instance.
///
/// Users must specify the exact class of the object that they want
/// to create (either an Xpetra::TpetraMap or an Xpetra::EpetraMap).
template<class LocalOrdinal,
         class GlobalOrdinal,
         class Node = typename Map<LocalOrdinal, GlobalOrdinal>::node_type>
class MapFactory
{


  private:


    //! Private constructor. This is a static class.
    MapFactory() {}


  public:


    //! Map constructor with Xpetra-defined contiguous uniform distribution.


    static Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node>>
    Build(UnderlyingLib                                 lib,
          global_size_t                                 numGlobalElements,
          GlobalOrdinal                                 indexBase,
          const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
          LocalGlobal lg = Xpetra::GloballyDistributed);


    //! Map constructor with a user-defined contiguous distribution.


    static Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node>>
    Build(UnderlyingLib                                 lib,
          global_size_t                                 numGlobalElements,
          size_t                                        numLocalElements,
          GlobalOrdinal                                 indexBase,
          const Teuchos::RCP<const Teuchos::Comm<int>>& comm);


    //! Map constructor with user-defined non-contiguous (arbitrary) distribution.


    static Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node>>
    Build(UnderlyingLib                                  lib,
          global_size_t                                  numGlobalElements,
          const Teuchos::ArrayView<const GlobalOrdinal>& elementList,
          GlobalOrdinal                                  indexBase,
          const Teuchos::RCP<const Teuchos::Comm<int>>&  comm);


    /*!
      @brief Map constructor generating degrees of freedom with numDofPerNode for given nodeMap

      @param[in] nodeMap Existing (node) map
      @param[in] numDofPerNode Number of DOFs per node for output map
      @param[in] gidOffset GID offset for output map
      @return Map

      \note This acts like a deep copy.
    */
    static Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node>>
    Build(const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& nodeMap,
          const LocalOrdinal numDofPerNode,
          const GlobalOrdinal gidOffset = Teuchos::ScalarTraits<GlobalOrdinal>::zero());


    static Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node>>
    Build(UnderlyingLib                                                         lib,
          global_size_t                                                         numGlobalElements,
          const Kokkos::View<const GlobalOrdinal*, typename Node::device_type>& indexList,
          GlobalOrdinal                                                         indexBase,
          const Teuchos::RCP<const Teuchos::Comm<int>>&                         comm);


    //! Create a locally replicated Map with the default node.
    static Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>
    createLocalMap(UnderlyingLib                                 lib,
                   size_t                                        numElements,
                   const Teuchos::RCP<const Teuchos::Comm<int>>& comm);


    //! Create a locally replicated Map with a specified node.


    static Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>
    createLocalMapWithNode(UnderlyingLib                                 lib,
                           size_t                                        numElements,
                           const Teuchos::RCP<const Teuchos::Comm<int>>& comm);


    //! Create a uniform, contiguous Map with a user-specified node.


    static Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>
    createUniformContigMapWithNode(UnderlyingLib                                 lib,
                                   global_size_t                                 numElements,
                                   const Teuchos::RCP<const Teuchos::Comm<int>>& comm);


    //! Create a uniform, contiguous Map with the default node.
    static Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>
    createUniformContigMap(UnderlyingLib                                 lib,
                           global_size_t                                 numElements,
                           const Teuchos::RCP<const Teuchos::Comm<int>>& comm);


    //! Create a (potentially) non-uniform, contiguous Map with the default node.
    static Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>
    createContigMap(UnderlyingLib                                 lib,
                    global_size_t                                 numElements,
                    size_t                                        localNumElements,
                    const Teuchos::RCP<const Teuchos::Comm<int>>& comm);


    //! Create a (potentially) non-uniform, contiguous Map with a user-specified node.


    static Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>
    createContigMapWithNode(UnderlyingLib                                 lib,
                            global_size_t                                 numElements,
                            size_t                                        localNumElements,
                            const Teuchos::RCP<const Teuchos::Comm<int>>& comm);

    //! Create a copy of the map, only using the new Comm object *if* the Comm would be valid
    // for this map.
    static Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
    copyMapWithNewComm(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> & oldmap,
                       const Teuchos::RCP<const Teuchos::Comm<int>>& newComm);
  
};      // class MapFactory



//////////////////////////////////////////////////////////////
///  X P E T R A   E P E T R A   S P E C I A L I Z A T I O N
//////////////////////////////////////////////////////////////


#if defined(HAVE_XPETRA_EPETRA)


#if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)


  template <>
  class MapFactory<int, int, EpetraNode>
  {

    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef EpetraNode Node;

  private:

    //! Private constructor. This is a static class.
    MapFactory();

  public:




    static RCP<Map<LocalOrdinal,GlobalOrdinal, Node> >
    Build (UnderlyingLib lib,
           global_size_t numGlobalElements,
           int indexBase,
           const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
           LocalGlobal lg=GloballyDistributed);




    static RCP<Map<LocalOrdinal,GlobalOrdinal, Node> >
    Build (UnderlyingLib lib,
           global_size_t numGlobalElements,
           size_t numLocalElements,
           int indexBase,
           const Teuchos::RCP<const Teuchos::Comm<int> > &comm);




    static RCP<Map<LocalOrdinal,GlobalOrdinal, Node> >
    Build(UnderlyingLib lib,
          global_size_t numGlobalElements,
          const Teuchos::ArrayView<const GlobalOrdinal> &elementList,
          int indexBase,
          const Teuchos::RCP<const Teuchos::Comm<int> > &comm);

    /*!
      @brief Map constructor transforming degrees of freedom for numDofPerNode

      @param[in] map Existing (node) map
      @param[in] numDofPerNode Number of DOFs per node for output map
      @param[in] gidOffset GID offset for output map
      @return Map

      \note This acts like a deep copy.
    */
    static Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node>>
    Build(const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& map,
          const LocalOrdinal numDofPerNode,
          const GlobalOrdinal gidOffset = Teuchos::ScalarTraits<GlobalOrdinal>::zero());


    static Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal, Node> >
    createLocalMap(UnderlyingLib lib,
                   size_t numElements,
                   const Teuchos::RCP< const Teuchos::Comm< int > > &comm);


    // TODO remove this


    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createLocalMapWithNode(UnderlyingLib lib,
                           size_t numElements,
                           const Teuchos::RCP< const Teuchos::Comm< int > > &comm);


    // TODO remove this



    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createUniformContigMapWithNode (UnderlyingLib lib, global_size_t numElements,
                                    const Teuchos::RCP< const Teuchos::Comm< int > > &comm);


    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node> >
    createUniformContigMap(UnderlyingLib lib,
                           global_size_t numElements,
                           const Teuchos::RCP< const Teuchos::Comm< int > > &comm);


    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createContigMap(UnderlyingLib lib,
                    global_size_t numElements,
                    size_t localNumElements,
                    const Teuchos::RCP< const Teuchos::Comm< int > > &comm);




    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createContigMapWithNode(UnderlyingLib lib,
                            global_size_t numElements,
                            size_t localNumElements,
                            const Teuchos::RCP< const Teuchos::Comm< int > > &comm);

    static Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
    copyMapWithNewComm(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> & oldmap,
                       const Teuchos::RCP<const Teuchos::Comm<int>>& newComm);


  };    // class MapFactory<int, int ... > specialization


#endif  // #if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)






// we need the Epetra specialization only if Epetra is enabled
#if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)


  template <>
  class MapFactory<int, long long, EpetraNode>
  {

    typedef int LocalOrdinal;
    typedef long long GlobalOrdinal;
    typedef EpetraNode Node;

  private:

    //! Private constructor. This is a static class.
    MapFactory();

  public:



    static RCP<Map<LocalOrdinal,GlobalOrdinal, Node> >
    Build (UnderlyingLib lib,
           global_size_t numGlobalElements,
           int indexBase,
           const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
           LocalGlobal lg=GloballyDistributed);




    static RCP<Map<LocalOrdinal,GlobalOrdinal, Node> >
    Build (UnderlyingLib lib,
           global_size_t numGlobalElements,
           size_t numLocalElements,
           int indexBase,
           const Teuchos::RCP<const Teuchos::Comm<int> > &comm);




    static RCP<Map<LocalOrdinal,GlobalOrdinal, Node> >
    Build(UnderlyingLib lib,
          global_size_t numGlobalElements,
          const Teuchos::ArrayView<const GlobalOrdinal> &elementList,
          int indexBase,
          const Teuchos::RCP<const Teuchos::Comm<int> > &comm);


    //! Map constructor transforming degrees of freedom
    //! for numDofPerNode this acts like a deep copy
    static Teuchos::RCP<Map<LocalOrdinal,GlobalOrdinal, Node> >
    Build(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
          LocalOrdinal numDofPerNode);


    static Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal, Node> >
    createLocalMap(UnderlyingLib lib,
                   size_t numElements,
                   const Teuchos::RCP< const Teuchos::Comm< int > > &comm);




    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createLocalMapWithNode(UnderlyingLib lib,
                           size_t numElements,
                           const Teuchos::RCP< const Teuchos::Comm< int > > &comm);




    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createUniformContigMapWithNode (UnderlyingLib lib, global_size_t numElements,
                                    const Teuchos::RCP< const Teuchos::Comm< int > > &comm);


    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node> >
    createUniformContigMap(UnderlyingLib lib,
                           global_size_t numElements,
                           const Teuchos::RCP< const Teuchos::Comm< int > > &comm);


    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createContigMap(UnderlyingLib lib,
                    global_size_t numElements,
                    size_t localNumElements,
                    const Teuchos::RCP< const Teuchos::Comm< int > > &comm);




    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createContigMapWithNode(UnderlyingLib lib,
                            global_size_t numElements,
                            size_t localNumElements,
                            const Teuchos::RCP< const Teuchos::Comm< int > > &comm);

    static Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
    copyMapWithNewComm(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> & oldmap,
                       const Teuchos::RCP<const Teuchos::Comm<int>>& newComm);


  };    // class MapFactory<int, long long, EpetraNode> specialization


#endif  // #if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)


#endif // #if defined(HAVE_XPETRA_EPETRA)


}      // namespace Xpetra


#define XPETRA_MAPFACTORY_SHORT

#endif      // XPETRA_MAPFACTORY_DECL_HPP

// TODO: removed unused methods
