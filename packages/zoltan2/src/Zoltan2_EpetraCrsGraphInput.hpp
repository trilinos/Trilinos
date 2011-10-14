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
};

/*! Zoltan2::EpetraCrsGraphInput
    \brief Provides access for Zoltan2 to Epetra_CrsGraph data and weights.

    EpetraCrsGraphInput provides access for Zoltan2 to an Epetra_CrsGraph and
    vertex and edges weights.  Most of the methods used by Zoltan2 are in the 
    XpetraGraphInput parent class.  The methods used by the Zoltan2 caller 
    to initialize the objects are in this class.

    The template parameter is the weight type.  Epetra local and global IDs 
    are ints.
*/

template <typename User>
class EpetraCrsGraphInput : public XpetraCrsGraphInput<User>{
private:

  RCP<const Epetra_CrsGraph> ingraph_;

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

  std::string inputAdapterName()const {return std::string("EpetraCrsGraph");}

  //~EpetraCrsGraphInput() { }

  /*! Constructor with graph only
      \param graph  the graph represented by this input adapter
   */
  EpetraCrsGraphInput(const RCP<const Epetra_CrsGraph> &graph):
    XpetraCrsGraphInput<User> (
      Teuchos::rcp(new Xpetra::EpetraCrsGraph(
        Teuchos::rcp_const_cast<Epetra_CrsGraph>(graph))))
  {
    ingraph_ = graph;
  }

#if 0
  /*! Constructor with weights
      \param graph  the graph represented by this input adapter
      \param vertexWeights are given in vertex local number order
      \param edgeWeights when queried, the Epetra_CrsGraph gives edges in a 
               certain order.  The edgeWeights follow this order, omitting
               any self edges.
   */
  EpetraCrsGraphInput(const RCP<const Epetra_CrsGraph> &graph,
                      const ArrayRCP<double> &vertexWeights,
                      const ArrayRCP<double> &edgeWeights):
    XpetraCrsGraphInput<int, int> (graph) {}
#endif

  /*! Access to graph that instantiated adapter
   */
  RCP<const Epetra_CrsGraph> getGraph() const 
  { 
    return ingraph_;
  }
};

} //namespace Zoltan2

#endif
