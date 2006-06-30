//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Alan Williams (william@sandia.gov)
                or Erik Boman    (egboman@sandia.gov)

************************************************************************
*/
//@HEADER

#ifndef _Isorropia_CostDescriber_hpp_
#define _Isorropia_CostDescriber_hpp_

#include <Isorropia_configdefs.hpp>
#include <Teuchos_ParameterList.hpp>

/** Isorropia is the namespace that contains isorropia's declarations
  for classes and functions.
*/
namespace Isorropia {

/** Interface (abstract base class) for describing the weights/costs
  associated with the vertices and/or edges of the object to be
  repartitioned and redistributed.

  This interface is intended to be queried by the partitioner.

  In many cases, a user may provide some weights but not others. E.g.,
  they may provide vertex weights but not edge weights, etc. Any weights
  which are not provided are assumed to have a default value of 1.0.
  The caller which needs the weights (i.e., the partitioner) can
  save time/work by querying whether "non-default" weights are
  specified, and if not, then simply use weights of 1.0 rather than
  making the call to obtain the weights arrays. The queries for
  whether weights are provided are the methods haveVertexWeights(),
  haveGraphEdgeWeights() and haveHypergraphEdgeWeights().
*/
class CostDescriber {
public:

  /** Destructor */
  virtual ~CostDescriber() {}

  /** Set parameters for the CostDescriber instance. The contents of the
      input paramlist object are copied into an internal ParameterList
      attribute. Instances of this interface should not retain a reference
      to the input ParameterList after this method returns.
  */
  virtual void setParameters(const Teuchos::ParameterList& paramlist) = 0;

  /** Query whether non-default vertex weights are present. If this
    function returns false, the caller can assume that vertex weights
    are all 1.0.
  */
  virtual bool haveVertexWeights() const = 0;

  /** Get the number of vertices. Vertices typically correspond to
    matrix rows.
  */
  virtual int getNumVertices() const = 0;

  /** Get the lists of vertex ids and weights.
  */
  virtual void getVertexWeights(int numVertices,
                                int* global_ids,
                                float* weights) const = 0;

  /** Query whether non-default graph edge weights are present.
     If this function returns false, the caller can assume that
     graph edge weights are all 1.0.
  */
  virtual bool haveGraphEdgeWeights() const = 0;

  /** Get the number of graph edges for a specified vertex.
     Graph edges typically correspond to matrix nonzeros.
  */
  virtual int getNumGraphEdges(int vertex_global_id) const = 0;

  /** Get the graph edge weights for a specified vertex.
  */
  virtual void getGraphEdgeWeights(int vertex_global_id,
                                   int num_neighbors,
                                   int* neighbor_global_ids,
                                   float* weights) const = 0;

  /** Query whether non-default hypergraph edge weights are present.
     If this function returns false, the caller can assume that
     hypergraph edge weights are all 1.0.
  */
  virtual bool haveHypergraphEdgeWeights() const = 0;

  /** Get the number of Hypergraph edges. Hypergraph edges typically
     correspond to matrix columns.
  */
  virtual int getNumHypergraphEdges() const = 0;

  /** Get the hypergraph edge weights.
  */
  virtual void getHypergraphEdgeWeights(int numEdges,
                                        int* global_ids,
                                        float* weights) const = 0;
};//class CostDescriber

}//namespace Isorropia

#endif

