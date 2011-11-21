#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Epetra_CrsMatrix.h>
#include <Epetra_CrsGraph.h>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsGraph.hpp>
#include <Zoltan2_Standards.hpp>

#ifndef ZOLTAN2_INPUTTRAITS_HPP
#define ZOLTAN2_INPUTTRAITS_HPP

namespace Zoltan2{

// Users who do not have templated input may find it convenient to create a:
// 
//   Zoltan2UserTypes<float, int, std::pair<int,int>, int, long> UserTypes
//
// and use it to instantiate the InputAdapter.
//
template <typename scalar, typename lid, typename gid, typename lno, typename gno>
class Zoltan2UserTypes{
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
  //   lid_t    :  user type that represents a locally unique identifier
  //               for data items.
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
  typedef int   lid_t;
  typedef long  gid_t;
  typedef Kokkos::DefaultNode::DefaultNodeType node_t;
  static inline std::string name() {return "InputAdapter";}
};

// Specialization for generic user input.
template <typename Scalar,
          typename LID,
          typename GID,
          typename LNO,
          typename GNO>
struct InputTraits<Zoltan2UserTypes<Scalar, LID, GID, LNO, GNO> >
{
  typedef Scalar        scalar_t;
  typedef LNO lno_t;
  typedef GNO gno_t;
  typedef LID lid_t;
  typedef GID gid_t;
  typedef Zoltan2::default_node_t node_t;
  static inline std::string name() {return "Zoltan2UserTypes";}
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
  typedef LocalOrdinal  lid_t;
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
  typedef LocalOrdinal  lid_t;
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
  typedef int lid_t;
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
  typedef LocalOrdinal  lid_t;
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
  typedef LocalOrdinal  lid_t;
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
  typedef int   lid_t;
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
  typedef LocalOrdinal  lid_t;
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
  typedef LocalOrdinal  lid_t;
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
  typedef int   lid_t;
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
  typedef LocalOrdinal  lid_t;
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
  typedef LocalOrdinal  lid_t;
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
  typedef int   lid_t;
  typedef int   gid_t;
  typedef Kokkos::DefaultNode::DefaultNodeType node_t;
  static inline std::string name() {return "Epetra_MultiVector";}
};

}  // namespace
#endif // ZOLTAN2_INPUTTRAITS_HPP
