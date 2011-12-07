// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// TODO: a test where user creates his/her traits and input adapter.

#include <Zoltan2_Standards.hpp>    // for Zoltan2::default_node_t
#include <Zoltan2_InputTraits.hpp>

// user data structure 
template<typename AppGID>
struct TestData{
  Teuchos::ArrayRCP<AppGID> gids;
};

// the InputTraits of our structure for Zoltan2
namespace Zoltan2{
template<>
template<typename AppGID>
struct InputTraits<struct TestData<AppGID> >
{
  typedef float scalar_t;
  typedef int lno_t;
  typedef long gno_t;
  typedef AppGID gid_t;
  typedef Zoltan2::default_node_t node_t;
};
}

#include <Zoltan2_IdentifierMap.hpp>

