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

************************************************************************
*/
//@HEADER

#ifndef _Isorropia_CostDescriber_hpp_
#define _Isorropia_CostDescriber_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_ParameterList.hpp>

/** Isorropia is the namespace that contains general definitions that
    apply to all partitioners and that contains abstract classes that 
    declare the methods and data to be supplied by specific partitioners.
*/

namespace Isorropia {

/** Interface (abstract base class) for describing the weights or costs
  associated with the vertices and/or edges of the object to be
  partitioned, ordered or colored.

  A CostDescriber object is created by the application.  If no CostDescriber is supplied by the
  application, sensible default weights should be used.

  The queries haveVertexWeights(), haveGraphEdgeWeights() and haveHypergraphEdgeWeights() 
  indicate whether the application has supplied weights to Isorropia.
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
    function returns false, and the operation is being performed by Zoltan, 
    the caller can assume that vertex weights are all 1.0.  (When Zoltan is
    turned off and Isorropia's simple internal linear partitioning is
    employed, then the default weight of a vertex is the number of non-zeros in the
    matrix row representing the vertex, or equivalently the number of neighbors the
    vertex has in the graph.)
  */
  virtual bool haveVertexWeights() const = 0;

  /** Get the number of vertices. Vertices typically correspond to
    matrix rows.
  */
  virtual int getNumVertices() const = 0;

  /** Get the lists of vertex ids and weights.  global_ids and weights are
      allocated by the caller.
  */
  virtual void getVertexWeights(int numVertices,
                                int* global_ids,
                                float* weights) const = 0;

  /** Query whether non-default graph edge weights are present.
     If this function returns false, the caller can assume that
     graph edge weights are all 1.0.  (If the operation is being performed
     by Isorropia's internal partitioner because Zoltan has been turned
     off with a parameter, then edge weights are ignored.)
  */
  virtual bool haveGraphEdgeWeights() const = 0;

  /** Get the number of graph edges for a specified vertex.
     Graph edges typically correspond to matrix nonzeros.
  */
  virtual int getNumGraphEdges(int vertex_global_id) const = 0;

  /** Get the graph edge weights for a specified vertex.  neighbor_global_ids and
      weights are allocated by the caller.
  */
  virtual void getGraphEdgeWeights(int vertex_global_id,
                                   int num_neighbors,
                                   int* neighbor_global_ids,
                                   float* weights) const = 0;

  /** Query whether non-default hypergraph edge weights are present.
     If this function returns false, the caller can assume that
     hypergraph edge weights are all 1.0.  (If Zoltan is not being used,
     but rather the Isorropia internal partitioner is being used, there
     is no concept of hyperedges.)
  */
  virtual bool haveHypergraphEdgeWeights() const = 0;

  /** Get the number of Hypergraph edges. Hypergraph edges typically
     correspond to matrix columns.
  */
  virtual int getNumHypergraphEdgeWeights() const = 0;

  /** Get the hypergraph edge weights.  The global_ids and weights are
      allocated by the caller.
  */
  virtual void getHypergraphEdgeWeights(int numEdges,
                                        int* global_ids,
                                        float* weights) const = 0;
};//class CostDescriber

}//namespace Isorropia

#endif

