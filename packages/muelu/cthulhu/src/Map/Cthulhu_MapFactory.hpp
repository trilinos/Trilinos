#ifndef CTHULHU_MAPFACTORY_HPP
#define CTHULHU_MAPFACTORY_HPP

#include <Kokkos_DefaultNode.hpp>

#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_Exceptions.hpp"

#include "Cthulhu_Map.hpp"

#ifdef HAVE_CTHULHU_TPETRA
#include "Cthulhu_TpetraMap.hpp"
#endif
#ifdef HAVE_CTHULHU_EPETRA
#include "Cthulhu_EpetraMap.hpp"
#endif

#include "Cthulhu_Debug.hpp"

// This factory creates Cthulhu::Map. User have to specify the exact class of object that he want to create (ie: a Cthulhu::TpetraMap or a Cthulhu::EpetraMap).

namespace Cthulhu {
  
  template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class MapFactory {
    
    typedef Map<LocalOrdinal, GlobalOrdinal, Node> Map;

#ifdef HAVE_CTHULHU_TPETRA
    typedef TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMap;
#endif

  private:
    //! Private constructor. This is a static class. 
    MapFactory() {}
    
  public:
    
    /** \brief Map constructor with Cthulhu-defined contiguous uniform distribution.
     *   The elements are distributed among nodes so that the subsets of global elements
     *   are non-overlapping and contiguous and as evenly distributed across the nodes as 
     *   possible.
     */
    //TODO: replace Tpetra::LocalGlobal by Cthulhu::LocalGlobal
    static RCP<Map> Build(UnderlyingLib lib, global_size_t numGlobalElements, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, 
                          Cthulhu::LocalGlobal lg=Cthulhu::GloballyDistributed, const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode()) {
#ifdef HAVE_CTHULHU_TPETRA
      if (lib == UseTpetra)
        return rcp( new TpetraMap(numGlobalElements, indexBase, comm, lg, node) );
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::RuntimeError,"Cannot create a map"); return Teuchos::null;
    }

    /** \brief Map constructor with a user-defined contiguous distribution.
     *  The elements are distributed among the nodes so that the subsets of global elements
     *  are non-overlapping and contiguous 
     *  
     *  If numGlobalElements == Teuchos::OrdinalTraits<global_size_t>::invalid(), it will be computed via a global communication.
     *  Otherwise, it must be equal to the sum of the local elements across all 
     *  nodes. This will only be verified if Trilinos was compiled with --enable-teuchos-debug.
     *  If this verification fails, a std::invalid_argument exception will be thrown.
     */
    static RCP<Map> Build(UnderlyingLib lib, global_size_t numGlobalElements, size_t numLocalElements, GlobalOrdinal indexBase, 
                          const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode()) {
#ifdef HAVE_CTHULHU_TPETRA
      if (lib == UseTpetra)       
        return rcp( new TpetraMap(numGlobalElements, numLocalElements, indexBase, comm, node) );
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::RuntimeError,"Cannot create a map"); return Teuchos::null;
    }
        
    /** \brief Map constructor with user-defined non-contiguous (arbitrary) distribution.
     *  
     *  If numGlobalElements == Teuchos::OrdinalTraits<global_size_t>::invalid(), it will be computed via a global communication.
     *  Otherwise, it must be equal to the sum of the local elements across all 
     *  nodes. This will only be verified if Trilinos was compiled with --enable-teuchos-debug.
     *  If this verification fails, a std::invalid_argument exception will be thrown.
     */
    static RCP<Map> Build(UnderlyingLib lib, global_size_t numGlobalElements, const Teuchos::ArrayView<const GlobalOrdinal> &elementList, GlobalOrdinal indexBase, 
                          const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode()) {
#ifdef HAVE_CTHULHU_TPETRA
      if (lib == UseTpetra) 
        return rcp( new TpetraMap(numGlobalElements, elementList, indexBase, comm, node) );
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::RuntimeError,"Cannot create a map"); return Teuchos::null;
    }

    /** \brief Non-member function to create a locally replicated Map with the default node.

    This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

    The Map is configured to use zero-based indexing.

    \relates TpetraMap
    */
    Teuchos::RCP<const Map>
    createLocalMap(UnderlyingLib lib, size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) { CTHULHU_DEBUG_ME;
#ifdef HAVE_CTHULHU_TPETRA
      if (lib == UseTpetra) 
        return rcp(new Cthulhu::TpetraMap<LocalOrdinal,GlobalOrdinal>(Tpetra::createLocalMap<LocalOrdinal,GlobalOrdinal>(numElements, comm)));
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::RuntimeError,"Cannot create a map"); return Teuchos::null;
    }

    /** \brief Non-member function to create a locally replicated Map with a specified node.

    The Map is configured to use zero-based indexing.

    \relates TpetraMap
    */
    Teuchos::RCP< const Map >
    createLocalMapWithNode(UnderlyingLib lib, size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) { CTHULHU_DEBUG_ME;
#ifdef HAVE_CTHULHU_TPETRA
      if (lib == UseTpetra) 
        return rcp(new TpetraMap(Tpetra::createLocalMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(numElements, comm, node)));
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::RuntimeError,"Cannot create a map"); return Teuchos::null;
    }

    /** \brief Non-member function to create a uniform, contiguous Map with a user-specified node.

    The Map is configured to use zero-based indexing.

    \relates TpetraMap
    */
    Teuchos::RCP< const Map >
    createUniformContigMapWithNode(UnderlyingLib lib, global_size_t numElements,
                                   const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) { CTHULHU_DEBUG_ME;
#ifdef HAVE_CTHULHU_TPETRA
      if (lib == UseTpetra) 
        return rcp(new TpetraMap(Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(numElements, comm, node)));
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::RuntimeError,"Cannot create a map"); return Teuchos::null;
    }

    /** \brief Non-member function to create a uniform, contiguous Map with the default node.

    This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

    The Map is configured to use zero-based indexing.

    \relates TpetraMap
    */
    Teuchos::RCP< const Map>
    createUniformContigMap(UnderlyingLib lib, global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) { CTHULHU_DEBUG_ME;
#ifdef HAVE_CTHULHU_TPETRA
      if (lib == UseTpetra) 
        return rcp(new Cthulhu::TpetraMap<LocalOrdinal,GlobalOrdinal>(Tpetra::createUniformContigMap<LocalOrdinal,GlobalOrdinal>(numElements, comm)));
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::RuntimeError,"Cannot create a map"); return Teuchos::null;
    }

    /** \brief Non-member function to create a (potentially) non-uniform, contiguous Map with the default node.

    This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

    The Map is configured to use zero-based indexing.

    \relates TpetraMap
    */
    Teuchos::RCP< const Map >
    createContigMap(UnderlyingLib lib, global_size_t numElements, size_t localNumElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) { CTHULHU_DEBUG_ME;
#ifdef HAVE_CTHULHU_TPETRA
      if (lib == UseTpetra) 
        return rcp(new Cthulhu::TpetraMap<LocalOrdinal,GlobalOrdinal>(Tpetra::createContigMap<LocalOrdinal,GlobalOrdinal>(numElements, localNumElements, comm)));
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::RuntimeError,"Cannot create a map"); return Teuchos::null;
    }

    /** \brief Non-member function to create a (potentially) non-uniform, contiguous Map with a user-specified node.

    The Map is configured to use zero-based indexing.

    \relates TpetraMap
    */
    Teuchos::RCP< const Map >
    createContigMapWithNode(UnderlyingLib lib, global_size_t numElements, size_t localNumElements, 
                            const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) { CTHULHU_DEBUG_ME;
#ifdef HAVE_CTHULHU_TPETRA
      if (lib == UseTpetra) 
        return rcp(new TpetraMap(Tpetra::createContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(numElements, localNumElements, comm, node)));
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::RuntimeError,"Cannot create a map"); return Teuchos::null;
    }

#ifdef CTHULHU_NOT_IMPLEMENTED
    /** \brief Non-member function to create a contiguous Map with user-defined weights and a user-specified node.

    The Map is configured to use zero-based indexing.

    \relates TpetraMap
    */
    Teuchos::RCP< const Map >
    createWeightedContigMapWithNode(UnderlyingLib lib, int thisNodeWeight, global_size_t numElements, 
                                    const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) { CTHULHU_DEBUG_ME;

    }
#endif
  };

  template <>
  class MapFactory<int, int, Kokkos::DefaultNode::DefaultNodeType> {
      
    typedef Map<int, int> Map;

#ifdef HAVE_CTHULHU_TPETRA
    typedef TpetraMap<int, int> TpetraMap;
#endif

  private:
    //! Private constructor. This is a static class. 
    MapFactory() {}
    
  public:
    
    /** \brief Map constructor with Cthulhu-defined contiguous uniform distribution.
     *   The elements are distributed among nodes so that the subsets of global elements
     *   are non-overlapping and contiguous and as evenly distributed across the nodes as 
     *   possible.
     */
    //TODO: replace Tpetra::LocalGlobal by Cthulhu::LocalGlobal
    static RCP<Map> Build(UnderlyingLib lib, global_size_t numGlobalElements, int indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, 
                          LocalGlobal lg=GloballyDistributed, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode()) {
#ifdef HAVE_CTHULHU_TPETRA
      if (lib == UseTpetra)
        return rcp( new TpetraMap(numGlobalElements, indexBase, comm, lg, node) );
#endif
#ifdef HAVE_CTHULHU_EPETRA
      if (lib == UseEpetra)
        return rcp( new EpetraMap(numGlobalElements, indexBase, comm, lg, node) );
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::RuntimeError,"Cannot create a map"); return Teuchos::null;
    }

    /** \brief Map constructor with a user-defined contiguous distribution.
     *  The elements are distributed among the nodes so that the subsets of global elements
     *  are non-overlapping and contiguous 
     *  
     *  If numGlobalElements == Teuchos::OrdinalTraits<global_size_t>::invalid(), it will be computed via a global communication.
     *  Otherwise, it must be equal to the sum of the local elements across all 
     *  nodes. This will only be verified if Trilinos was compiled with --enable-teuchos-debug.
     *  If this verification fails, a std::invalid_argument exception will be thrown.
     */
    static RCP<Map> Build(UnderlyingLib lib, global_size_t numGlobalElements, size_t numLocalElements, int indexBase, 
                          const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode()) {
#ifdef HAVE_CTHULHU_TPETRA
      if (lib == UseTpetra)       
        return rcp( new TpetraMap(numGlobalElements, numLocalElements, indexBase, comm, node) );

#endif
#ifdef HAVE_CTHULHU_EPETRA
      if (lib == UseEpetra)
        return rcp( new EpetraMap(numGlobalElements, numLocalElements, indexBase, comm, node) );
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::RuntimeError,"Cannot create a map"); return Teuchos::null;
    }
        
    /** \brief Map constructor with user-defined non-contiguous (arbitrary) distribution.
     *  
     *  If numGlobalElements == Teuchos::OrdinalTraits<global_size_t>::invalid(), it will be computed via a global communication.
     *  Otherwise, it must be equal to the sum of the local elements across all 
     *  nodes. This will only be verified if Trilinos was compiled with --enable-teuchos-debug.
     *  If this verification fails, a std::invalid_argument exception will be thrown.
     */
    static RCP<Map> Build(UnderlyingLib lib, global_size_t numGlobalElements, const Teuchos::ArrayView<const int> &elementList, int indexBase, 
                          const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode()) {
#ifdef HAVE_CTHULHU_TPETRA
      if (lib == UseTpetra) 
        return rcp( new TpetraMap(numGlobalElements, elementList, indexBase, comm, node) );

#endif
#ifdef HAVE_CTHULHU_EPETRA
      if (lib == UseEpetra)
        return rcp( new EpetraMap(numGlobalElements, elementList, indexBase, comm, node) );
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::RuntimeError,"Cannot create a map"); return Teuchos::null;
    }

    /** \brief Non-member function to create a locally replicated Map with the default node.

    This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

    The Map is configured to use zero-based indexing.

    \relates TpetraMap
    */
    static Teuchos::RCP<const Map>
    createLocalMap(UnderlyingLib lib, size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) { CTHULHU_DEBUG_ME;
#ifdef HAVE_CTHULHU_TPETRA
      if (lib == UseTpetra) 
        return rcp( new TpetraMap(Tpetra::createLocalMap<int,int>(numElements, comm)));
#endif
#ifdef HAVE_CTHULHU_EPETRA
      if (lib == UseEpetra)
        return MapFactory<int, int, Kokkos::DefaultNode::DefaultNodeType>::createLocalMapWithNode(lib, numElements, comm, Kokkos::DefaultNode::getDefaultNode());
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::RuntimeError,"Cannot create a map"); return Teuchos::null;
    }

    /** \brief Non-member function to create a locally replicated Map with a specified node.

    The Map is configured to use zero-based indexing.

    \relates TpetraMap
    */
    static Teuchos::RCP< const Map >
    createLocalMapWithNode(UnderlyingLib lib, size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node) { CTHULHU_DEBUG_ME;
#ifdef HAVE_CTHULHU_TPETRA
      if (lib == UseTpetra) 
        return rcp(new TpetraMap(Tpetra::createLocalMapWithNode<int,int,Kokkos::DefaultNode::DefaultNodeType>(numElements, comm, node)));
#endif
#ifdef HAVE_CTHULHU_EPETRA
      if (lib == UseEpetra)
        {
          Teuchos::RCP< EpetraMap > map;
          map = Teuchos::rcp( new EpetraMap((Cthulhu::global_size_t)numElements, // num elements, global and local
                                            0,                                   // index base is zero
                                            comm, LocallyReplicated, node));
          return map.getConst();
        }
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::RuntimeError,"Cannot create a map"); return Teuchos::null;
    }

    /** \brief Non-member function to create a uniform, contiguous Map with a user-specified node.

    The Map is configured to use zero-based indexing.

    \relates TpetraMap
    */
    static Teuchos::RCP< const Map >
    createUniformContigMapWithNode(UnderlyingLib lib, global_size_t numElements,
                                   const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node) { CTHULHU_DEBUG_ME;
#ifdef HAVE_CTHULHU_TPETRA
      if (lib == UseTpetra) 
        return rcp(new TpetraMap(Tpetra::createUniformContigMapWithNode<int,int,Kokkos::DefaultNode::DefaultNodeType>(numElements, comm, node)));
#endif
#ifdef HAVE_CTHULHU_EPETRA
      if (lib == UseEpetra)
        {
          Teuchos::RCP< EpetraMap > map;
          map = Teuchos::rcp( new EpetraMap(numElements,        // num elements, global and local
                                            0,                  //index base is zero
                                            comm, GloballyDistributed, node));
          return map.getConst();
        }
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::RuntimeError,"Cannot create a map"); return Teuchos::null;
    }

    /** \brief Non-member function to create a uniform, contiguous Map with the default node.

    This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

    The Map is configured to use zero-based indexing.

    \relates TpetraMap
    */
    static Teuchos::RCP< const Map>
    createUniformContigMap(UnderlyingLib lib, global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) { CTHULHU_DEBUG_ME;
#ifdef HAVE_CTHULHU_TPETRA
      if (lib == UseTpetra) 
        return rcp( new TpetraMap(Tpetra::createUniformContigMap<int,int>(numElements, comm)));
#endif
#ifdef HAVE_CTHULHU_EPETRA
      if (lib == UseEpetra)
        return MapFactory<int, int, Kokkos::DefaultNode::DefaultNodeType>::createUniformContigMapWithNode(lib, numElements, comm, Kokkos::DefaultNode::getDefaultNode());
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::RuntimeError,"Cannot create a map"); return Teuchos::null;
    }

    /** \brief Non-member function to create a (potentially) non-uniform, contiguous Map with the default node.

    This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

    The Map is configured to use zero-based indexing.

    \relates TpetraMap
    */
    static Teuchos::RCP< const Map >
    createContigMap(UnderlyingLib lib, global_size_t numElements, size_t localNumElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) { CTHULHU_DEBUG_ME;
#ifdef HAVE_CTHULHU_TPETRA
      if (lib == UseTpetra) 
        return rcp( new TpetraMap(Tpetra::createContigMap<int,int>(numElements, localNumElements, comm)));
#endif
#ifdef HAVE_CTHULHU_EPETRA
      if (lib == UseEpetra)
        return MapFactory<int, int, Kokkos::DefaultNode::DefaultNodeType>::createContigMapWithNode(lib, numElements, localNumElements, comm, Kokkos::DefaultNode::getDefaultNode() );
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::RuntimeError,"Cannot create a map"); return Teuchos::null;
    }

    /** \brief Non-member function to create a (potentially) non-uniform, contiguous Map with a user-specified node.

    The Map is configured to use zero-based indexing.

    \relates TpetraMap
    */
    static Teuchos::RCP< const Map >
    createContigMapWithNode(UnderlyingLib lib, global_size_t numElements, size_t localNumElements, 
                            const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node) { CTHULHU_DEBUG_ME;
#ifdef HAVE_CTHULHU_TPETRA
      if (lib == UseTpetra) 
        return rcp( new TpetraMap(Tpetra::createContigMapWithNode<int,int,Kokkos::DefaultNode::DefaultNodeType>(numElements, localNumElements, comm, node)));
#endif
#ifdef HAVE_CTHULHU_EPETRA
      if (lib == UseEpetra)
        {
          Teuchos::RCP< EpetraMap > map;
          map = Teuchos::rcp( new EpetraMap(numElements,localNumElements,
                                            0,  // index base is zero
                                            comm, node) );
          return map.getConst();
        }
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::RuntimeError,"Cannot create a map"); return Teuchos::null;
    }

#ifdef CTHULHU_NOT_IMPLEMENTED
    /** \brief Non-member function to create a contiguous Map with user-defined weights and a user-specified node.

    The Map is configured to use zero-based indexing.

    \relates TpetraMap
    */
    static Teuchos::RCP< const Map >
    createWeightedContigMapWithNode(UnderlyingLib lib, int thisNodeWeight, global_size_t numElements, 
                                    const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node) { CTHULHU_DEBUG_ME;
      
    }
#endif

  };

}

#define CTHULHU_MAPFACTORY_SHORT
#endif
