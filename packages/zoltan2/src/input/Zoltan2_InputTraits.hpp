// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_InputTraits.hpp

    \brief Traits for application input objects.
*/

#ifndef ZOLTAN2_INPUTTRAITS_HPP
#define ZOLTAN2_INPUTTRAITS_HPP

#include <Zoltan2_Standards.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Epetra_CrsMatrix.h>
#include <Epetra_CrsGraph.h>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsGraph.hpp>

namespace Zoltan2{

/*! \brief A class that can be the User template argument for an InputAdapter.
 *
 *  BasicUserTypes is a convenience class that provides a simple way
 *  to supply the template argument when constructing an InputAdapter.
 *
 *  Typically, a C++ user will have a templated or non-templated class that
 *  represents his or her input data.  He or she will create an InputTraits
 *  specialization for this class, and then supply the class as the template
 *  argument of the InputAdapter class.  (Commonly used Trilinos classes
 *  that represent vectors, graphs and matrices already have InputTraits
 *  specializations.)
 *
 *  This makes sense if you are a C++ programmer who uses templates, but
 *  is a great deal of trouble if you are not.  In this case you can
 *  construct your InputAdapter as follows.
 *
 *  Suppose you want to construct at Zoltan2::BasicVectorInput object and
 *  you use \c float for vector values in your application, \c long for 
 *  global identifiers, and \c int for local indices. 
 *
 *  You need to determine an integral data type that Zoltan2 can use internally
 *  for global identifiers. Often this is the same data type that you use for
 *  this purpose, but if you have non-integral global identifiers (such as
 *  std::pair<int, int>) then Zoltan2 needs you to supply an integral data
 *  type that is large enough to enumerate your global number of objects.
 *  In this example let's say that you want Zoltan2 to use \c unsigned \c long
 *  for its global Identifiers.  Then you would do the following:
 *
\code
   typedef BasicUserTypes<float, long, int, unsigned long> myTypes;
   Zoltan2::BasicVectorInput<myTypes> myInput({constructor arguments})
\endcode
 *
 * In particular, the BasicUserTypes template parameters are:

    \li \c scalar is the data type for element values, weights and coordinates.
    \li \c gid is the data type used by the application for global Ids.  If the application's global Id data type is a Teuchos Ordinal, then \c gid and \c gno can the same.  Otherwise, the application global Ids will be mapped to Teuchos Ordinals for use by Zoltan2 internally.  (Teuchos Ordinals are those data types for which traits are defined in Trilinos/packages/teuchos/src/Teuchos_OrdinalTraits.hpp.)
    \li \c lno is the integral data type used by the application and by Zoltan2 for local indices and local counts.
    \li \c gno is the integral data type used by Zoltan2 to represent global indices and global counts.
 */  

template <typename scalar, typename gid, typename lno, typename gno>
class BasicUserTypes{
};

template <typename User>
struct InputTraits {
  // Input Adapter implementations must provide the following typedefs
  // for use in Zoltan2:
  //   scalar_t :  weights and coordinates
  //   lno_t    :  ordinal (e.g., int, long, int64_t) that can represent
  //               the number of local data items.
  //   gno_t    :  ordinal (e.g., int, long, int64_t) that can represent
  //               the number of global data items.
  //   gid_t    :  user type that represents a globally unique identifier
  //               for data items.
  //   node_t   :  Kokkos node.
  //   name     :  Name of the User data being used; currently used
  //               only for debugging.
  // 
  // Default typedefs are included here. If a specialization of User is
  // not provided, these types will be used.
  typedef float scalar_t;
  typedef int   lno_t;
  typedef long  gno_t;
  typedef long  gid_t;
  typedef Kokkos::DefaultNode::DefaultNodeType node_t;
  static inline std::string name() {return "InputAdapter";}
};

// Specialization for generic user input.
template <typename Scalar,
          typename GID,
          typename LNO,
          typename GNO>
struct InputTraits<BasicUserTypes<Scalar, GID, LNO, GNO> >
{
  typedef Scalar        scalar_t;
  typedef LNO lno_t;
  typedef GNO gno_t;
  typedef GID gid_t;
  typedef Zoltan2::default_node_t node_t;
  static inline std::string name() {return "BasicUserTypes";}
};

// Specialization for Xpetra::CrsMatrix.
// KDDKDD Do we need specializations for Xpetra::EpetraCrsMatrix and
// KDDKDD Xpetra::TpetraCrsMatrix?
template <typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct InputTraits<Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
{
  typedef Scalar        scalar_t;
  typedef LocalOrdinal  lno_t;
  typedef GlobalOrdinal gno_t;
  typedef GlobalOrdinal gid_t;
  typedef Node          node_t;
  static inline std::string name() {return "Xpetra::CrsMatrix";}
};

// Specialization for Tpetra::CrsMatrix
template <typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct InputTraits<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
{
  typedef Scalar        scalar_t;
  typedef LocalOrdinal  lno_t;
  typedef GlobalOrdinal gno_t;
  typedef GlobalOrdinal gid_t;
  typedef Node          node_t;
  static inline std::string name() {return "Tpetra::CrsMatrix";}
};

// Epetra_CrsMatrix
template < >
struct InputTraits<Epetra_CrsMatrix>
{
  typedef double scalar_t;
  typedef int lno_t;
  typedef int gno_t;
  typedef int gid_t;
  typedef Kokkos::DefaultNode::DefaultNodeType node_t;
  static inline std::string name() {return "Epetra_CrsMatrix";}
};

// Xpetra::CrsGraph
// KDDKDD:  Do we need specializations for Xpetra::EpetraCrsGraph and
// KDDKDD:  Xpetra::TpetraCrsGraph
template <typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct InputTraits<Xpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >
{
  typedef float         scalar_t;
  typedef LocalOrdinal  lno_t;
  typedef GlobalOrdinal gno_t;
  typedef GlobalOrdinal gid_t;
  typedef Node          node_t;
  static inline std::string name() {return "Xpetra::CrsGraph";}
};

// Tpetra::CrsGraph
template <typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct InputTraits<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >
{
  typedef float         scalar_t;
  typedef LocalOrdinal  lno_t;
  typedef GlobalOrdinal gno_t;
  typedef GlobalOrdinal gid_t;
  typedef Node          node_t;
  static inline std::string name() {return "Tpetra::CrsGraph";}
};

// Epetra_CrsGraph
template < >
struct InputTraits<Epetra_CrsGraph>
{
  typedef float scalar_t;
  typedef int   lno_t;
  typedef int   gno_t;
  typedef int   gid_t;
  typedef Kokkos::DefaultNode::DefaultNodeType node_t;
  static inline std::string name() {return "Epetra_CrsGraph";}
};

template <typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct InputTraits<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
{
  typedef Scalar        scalar_t;
  typedef LocalOrdinal  lno_t;
  typedef GlobalOrdinal gno_t;
  typedef GlobalOrdinal gid_t;
  typedef Node          node_t;
  static inline std::string name() {return "Xpetra::Vector";}
};

template <typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct InputTraits<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
{
//TODO A Tpetra::Vector is a Tpetra::MultiVector - can we just
//  define MultiVector traits only?  Ditto with Xpetra.  Test this
  typedef Scalar        scalar_t;
  typedef LocalOrdinal  lno_t;
  typedef GlobalOrdinal gno_t;
  typedef GlobalOrdinal gid_t;
  typedef Node          node_t;
  static inline std::string name() {return "Tpetra::Vector";}
};

// Epetra_Vector
template < >
struct InputTraits<Epetra_Vector>
{
  typedef double scalar_t;
  typedef int   lno_t;
  typedef int   gno_t;
  typedef int   gid_t;
  typedef Kokkos::DefaultNode::DefaultNodeType node_t;
  static inline std::string name() {return "Epetra_Vector";}
};


template <typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct InputTraits<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
{
  typedef Scalar        scalar_t;
  typedef LocalOrdinal  lno_t;
  typedef GlobalOrdinal gno_t;
  typedef GlobalOrdinal gid_t;
  typedef Node          node_t;
  static inline std::string name() {return "Xpetra::MultiVector";}
};

template <typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct InputTraits<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
{
  typedef Scalar        scalar_t;
  typedef LocalOrdinal  lno_t;
  typedef GlobalOrdinal gno_t;
  typedef GlobalOrdinal gid_t;
  typedef Node          node_t;
  static inline std::string name() {return "Tpetra::MultiVector";}
};


template < >
struct InputTraits<Epetra_MultiVector>
{
  typedef double scalar_t;
  typedef int   lno_t;
  typedef int   gno_t;
  typedef int   gid_t;
  typedef Kokkos::DefaultNode::DefaultNodeType node_t;
  static inline std::string name() {return "Epetra_MultiVector";}
};

}  // namespace
#endif // ZOLTAN2_INPUTTRAITS_HPP
