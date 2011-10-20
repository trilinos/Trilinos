// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_EpetraCrsGraphInput.hpp

    \brief An input adapter for a Epetra_CrsGraph.
*/

#ifndef _ZOLTAN2_EPETRACRSGRAPHINPUT_HPP_
#define _ZOLTAN2_EPETRACRSGRAPHINPUT_HPP_

#include <Zoltan2_XpetraCrsGraphInput.hpp>
#include <Epetra_CrsGraph.h>

namespace Zoltan2 {

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
  static inline RCP<const Xpetra::CrsGraph<lno_t,gno_t,node_t> >
    convertToXpetra(const RCP<const Epetra_CrsGraph> &a)
    {
      return rcp(new Xpetra::EpetraCrsGraph(
                             rcp_const_cast<Epetra_CrsGraph>(a)));
    }
};

/*! Zoltan2::EpetraCrsGraphInput
    \brief Provides access for Zoltan2 to Epetra_CrsGraph data and weights.

    EpetraCrsGraphInput provides access for Zoltan2 to an Epetra_CrsGraph and
    vertex and edges weights.  Most of the methods used by Zoltan2 are in the 
    XpetraGraphInput parent class.  The methods used by the Zoltan2 caller 
    to initialize the objects are in this class.

    The template parameter is the Epetra_CrsGraph (the User data type).  
    Epetra local and global IDs are ints, and weights are double.
*/

template <typename User>
class EpetraCrsGraphInput : public XpetraCrsGraphInput<User>{
private:

  RCP<const User> ingraph_;

#if 0
  int vtxWeightDim_;
  int edgeWeightDim_;
  int coordinateDim_;
  std::vector<double> edgeWgt_;
  std::vector<double> vertexWgt_;
  std::vector<double> xyz_;
#endif

public:

  typedef typename InputAdapter<User>::scalar_t scalar_t;
  typedef typename InputAdapter<User>::lno_t    lno_t;
  typedef typename InputAdapter<User>::gno_t    gno_t;
  typedef typename InputAdapter<User>::lid_t    lid_t;
  typedef typename InputAdapter<User>::gid_t    gid_t;
  typedef typename InputAdapter<User>::node_t   node_t;

  // TODO - make the name part of the Traits definition
  std::string inputAdapterName()const {return std::string("EpetraCrsGraph");}

  //~EpetraCrsGraphInput() { }

  /*! Constructor with graph only
      \param graph  the graph represented by this input adapter
   */
  EpetraCrsGraphInput(const RCP<const User> &graph):
    XpetraCrsGraphInput<User>(graph)
  {
    ingraph_ = graph;
  }

  /*! Access to graph that instantiated adapter
   */
  RCP<const User> getGraph() const 
  { 
    return ingraph_;
  }
};

} //namespace Zoltan2

#endif
