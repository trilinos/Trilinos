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

/*! Zoltan2::EpetraCrsGraphInput
    \brief Provides access for Zoltan2 to Epetra_CrsGraph data and weights.

    EpetraCrsGraphInput provides access for Zoltan2 to an Epetra_CrsGraph and
    vertex and edges weights.  Most of the methods used by Zoltan2 are in the 
    XpetraGraphInput parent class.  The methods used by the Zoltan2 caller 
    to initialize the objects are in this class.

    The template parameter is the weight type.  Epetra local and global IDs 
    are ints.
*/

class EpetraCrsGraphInput : public XpetraCrsGraphInput<int, int>{
private:

  RCP<const Epetra_CrsGraph> _ingraph;

#if 0
  int _vtxWeightDim;
  int _edgeWeightDim;
  int _coordinateDim;
  std::vector<double> _edgeWgt;
  std::vector<double> _vertexWgt;
  std::vector<double> _xyz;
#endif

public:

  std::string inputAdapterName()const {return std::string("EpetraCrsGraph");}

  //~EpetraCrsGraphInput() { }

  /*! Constructor with graph only
      \param graph  the graph represented by this input adapter
   */
  EpetraCrsGraphInput(const RCP<const Epetra_CrsGraph> graph):
    XpetraCrsGraphInput<int, int> (
      Teuchos::rcp(new Xpetra::EpetraCrsGraph(
        Teuchos::rcp_const_cast<Epetra_CrsGraph>(graph))))
  {
    _ingraph = graph;
  }

#if 0
  /*! Constructor with weights
      \param graph  the graph represented by this input adapter
      \param vertexWeights are given in vertex local number order
      \param edgeWeights when queried, the Epetra_CrsGraph gives edges in a 
               certain order.  The edgeWeights follow this order, omitting
               any self edges.
   */
  EpetraCrsGraphInput(RCP<const Epetra_CrsGraph> graph,
                      ArrayRCP<double> vertexWeights,
                      ArrayRCP<double> edgeWeights):
    XpetraCrsGraphInput<int, int> (graph) {}
#endif

  /*! Access to graph that instantiated adapter
   */
  RCP<const Epetra_CrsGraph> getGraph() const 
  { 
    return _ingraph;
  }
};

} //namespace Zoltan2

#endif
