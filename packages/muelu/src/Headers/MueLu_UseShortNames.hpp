// Helper to get ride of template parameters

// This file can be use for two purpose:
// 1) As an header of a user program.
//    In this case, this file must be include *after* other headers
//    and type ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps must be defined.
//    Note also that there is no #ifndef/#endif to protect again the multiple inclusion of this file.
//    User should create is own header file including this one:
//
//    Example:
//     #ifndef MY_HEADER
//     #define MY_HEADER
//     #include <MueLu_UseDefaultTypes.hpp>
//     #include <MueLu_UseShortNames.hpp>
//     #endif
//
// 2) Inside of MueLu to enhance the readability.
//
// template <class Scalar, class LocalOrdinal=int, class GlobalOrdinal=LocalOrdinal, class Node=Kokkos::DefaultNode::DefaultNodeType>
//  class TpetraMultiVector : public virtual Cthulhu::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
// 
//  #include <MueLu_UseShortNames.hpp>
//
//  myMethod(RCP<const Map> & map) { [...] } // instead of myMethod(RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map)
//
//  [...]
//
// } 
//

#include "Cthulhu_Classes.hpp"
#include <Cthulhu_UseShortNames.hpp>
#include "MueLu_Classes"

// New definition of types using the types ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps of the current context.
#ifdef MUELU_SAPFACTORY_HPP
typedef MueLu::SaPFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>      SaPFactory;
#endif

#ifdef MUELU_TRANSPFACTORY_HPP
typedef MueLu::TransPFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>   TransPFactory;
#endif

#ifdef MUELU_RAPFACTORY_HPP
typedef MueLu::RAPFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>      RAPFactory;
#endif

#ifdef MUELU_SMOOTHERFACTORY_HPP
typedef MueLu::SmootherFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> SmootherFactory;
#endif
