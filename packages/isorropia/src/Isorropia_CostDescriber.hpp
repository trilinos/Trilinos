//@HEADER
//************************************************************************
//
//              Isorropia: Partitioning and Load Balancing Package
//                Copyright (2006) Sandia Corporation
//
//Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
//license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//************************************************************************
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
  associated with the vertices and/or edges or hyperedges of the object to be
  partitioned, ordered or colored.

  A CostDescriber object is created by the application.  If no CostDescriber is supplied by the
  application, sensible default weights should be used.

*/
class CostDescriber {
public:

  /** Destructor */
  virtual ~CostDescriber() {}

private:
  /** Set parameters for the CostDescriber instance. The contents of the
      input paramlist object are copied into an internal ParameterList
      attribute. Instances of this interface should not retain a reference
      to the input ParameterList after this method returns.
  */
  virtual void setParameters(const Teuchos::ParameterList& paramlist) = 0;

  /** Query whether vertex weights have been supplied by the application.

      \return returns true if the application has supplied vertex weights
              with the CostDescriber, false otherwise
  */
  virtual bool haveVertexWeights() const = 0;

  /** Get the number of vertices for which this process supplied
      vertex weights. Vertices typically correspond to matrix rows.

      \return returns the number of vertices on this process for which
              the CostDescriber has vertex weights
  */
  virtual int getNumVertices() const = 0;

  /** Get lists of the vertex ids and weights supplied by this process.  

     \param numVertices size of global_ids and weights arrays

     \param global_ids pointer to an array of vertex global IDs,
                         allocated by the caller.
     \param weights pointer to an array of vertex weights,
                         allocated by the caller.
  */
  virtual void getVertexWeights(int numVertices,
                                int* global_ids,
                                float* weights) const = 0;

  /** Query whether graph edge weights have been supplied by the application.

      \return returns true if the application has supplied graph edge weights
              with the CostDescriber, false otherwise
  */
  virtual bool haveGraphEdgeWeights() const = 0;

  /** Get the number of graph edges for a specified vertex.
     Graph edges typically correspond to matrix nonzeros.

     \param vertex_global_id  The global ID for the vertex (on this process)
                               for which the number of edges is desired

     \return the number of graph edges supplied by this process for this vertex
  */
  virtual int getNumGraphEdges(int vertex_global_id) const = 0;

  /** Get the graph edge weights for a specified vertex.  

      \param vertex_global_id  vertex global ID (on this process) for which
                               edge information is requested

      \param num_neighbors     size for which neighbor_global_ids and weights
                               had been preallocated

      \param neighbor_global_ids  buffer allocated by caller, on return will
                                   contain a list of neighbor vertex global IDs

      \param weights  buffer allocated by caller, on return will contain a
                      weight for each edge indicated in neighbor_global_ids

  */
  virtual void getGraphEdgeWeights(int vertex_global_id,
                                   int num_neighbors,
                                   int* neighbor_global_ids,
                                   float* weights) const = 0;

  /** Query whether hypergraph edge weights have been supplied by the application.

      \return returns true if the application has supplied hypergraph edge weights
              with the CostDescriber, false otherwise
  */
  virtual bool haveHypergraphEdgeWeights() const = 0;

  /** Get the number of Hypergraph edges. Hypergraph edges typically
     correspond to matrix columns.

    \return returns the number of hypergraph edge weights supplied by this process
  */
  virtual int getNumHypergraphEdgeWeights() const = 0;

  /** Get the hypergraph edge weights that were supplied by this process

      \param numEdges size for which global_ids and weights
                               had been preallocated

      \param  global_ids  buffer allocated by caller, on return will
                                   contain a list of hyperedge global IDs

      \param weights  buffer allocated by caller, on return will contain a
                      weight for each hyperedge indicated in global_ids
  */

  virtual void getHypergraphEdgeWeights(int numEdges,
                                        int* global_ids,
                                        float* weights) const = 0;
};//class CostDescriber

}//namespace Isorropia

#endif


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

