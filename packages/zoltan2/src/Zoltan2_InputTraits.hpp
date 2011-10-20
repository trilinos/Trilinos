#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Epetra_CrsMatrix.h>
#include <Epetra_CrsGraph.h>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsGraph.hpp>

#ifndef ZOLTAN2_INPUTTRAITS_HPP
#define ZOLTAN2_INPUTTRAITS_HPP

namespace Zoltan2{

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

#ifdef KDDKDD_CAN_WE_LEAVE_SPECIALIZATIONS_IN_SPECIFIC_INPUT_ADAPTER_FILES
// KDDKDD For now, I left the specializations in the same files as their
// KDDKDD specific input adapters.  That way, I could add methods (e.g.,
// KDDKDD convertToXpetra) to the InputTraits
// KDDKDD for specific specializations (Epetra, Tpetra, Xpetra) without
// KDDKDD this interface file becoming confusing for users (who might 
// KDDKDD think they need to implement convertToXpetra for every input adapter.
// KDDKDD Also, I envision users would implement their own Traits and Adapters
// KDDKDD in one file, rather than edit this file.
// KDDKDD If this idea is faulty, we can restore this code, adding the 
// KDDKDD extra methods for Epetra, Tpetra and Xpetra specializations.

template < >
struct InputTraits<Epetra_CrsGraph>
{
  typedef float scalar_t;
  typedef int   lno_t;
  typedef int   gno_t;
  typedef int   lid_t;
  typedef int   gid_t;
  typedef Kokkos::DefaultNode::DefaultNodeType node_t;
};

template < >
struct InputTraits<Epetra_CrsMatrix>
{
  typedef double scalar_t;
  typedef int lno_t; 
  typedef int gno_t; 
  typedef int lid_t;
  typedef int gid_t;
  typedef Kokkos::DefaultNode::DefaultNodeType node_t; 
};  

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
};  

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
};

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
};

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
};

#endif  //KDDKDD_CAN_WE_LEAVE_SPECIALIZATIONS_IN_SPECIFIC_INPUT_ADAPTER_FILES
}

#endif
