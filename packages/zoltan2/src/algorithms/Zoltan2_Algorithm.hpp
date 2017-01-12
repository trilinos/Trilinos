// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_BaseAdapter.hpp
    \brief Defines the Adapter interface for accessing user data.
*/

#ifndef _ZOLTAN2_ALGORITHM_HPP_
#define _ZOLTAN2_ALGORITHM_HPP_

namespace Zoltan2 {
template <typename Adapter>
class Algorithm;
}

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_ColoringSolution.hpp>
#include <Zoltan2_OrderingSolution.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_MappingSolution.hpp>
#include <Zoltan2_CoordinatePartitioningGraph.hpp>
#include <Zoltan2_PartitionTree.hpp>


namespace Zoltan2 {

//! \brief Algorithm defines the base class for all algorithms.
//
//  Algorithms do not have to implement all methods in the Algorithm base
//  class.  They should implement only those methods that are relevant.
//  For example AlgScotch might implement partition() and order(), while
//  AlgMJ might implement partition() and boxAssign().
//  Default implementations throw a "not implemented" error

template <typename Adapter>
class Algorithm {

public:

  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::part_t part_t;

  // Virtual destructor needed to avoid undefined behavior and compiler warnings
  virtual ~Algorithm() {}

  //! \brief Ordering method
  virtual int localOrder(const RCP<LocalOrderingSolution<lno_t> > &solution)
  {
    Z2_THROW_NOT_IMPLEMENTED
  }

  //! \brief Ordering method
  virtual int globalOrder(const RCP<GlobalOrderingSolution<gno_t> > &solution)
  {
    Z2_THROW_NOT_IMPLEMENTED 
  }
  
  //! \brief Coloring method
  virtual void color(const RCP<ColoringSolution<Adapter> > &solution) 
  {
    Z2_THROW_NOT_IMPLEMENTED
  }
  
  //! \brief Matching method
  virtual void match() { 
    Z2_THROW_NOT_IMPLEMENTED 
  }

  //! \brief Partitioning method
  virtual void partition(const RCP<PartitioningSolution<Adapter> > &solution) 
  {
    Z2_THROW_NOT_IMPLEMENTED
  }

  //! \brief Mapping method
  virtual void map(const RCP<MappingSolution<Adapter> > &solution) 
  {
    Z2_THROW_NOT_IMPLEMENTED
  }

  //! \brief  for partitioning methods, return partition tree of the
  //          computed parts
  //
  virtual const std::vector<partitionTreeNode<part_t> > &
  getPartitionTreeNodes() const
  {
    return partitionTree;
  }

  //! \brief  for partitioning methods, set partition tree of the
  //          computed parts.
  //
  void setPartitionTreeNodes(
    const std::vector<partitionTreeNode<part_t> > & setPartitionTree)
  {
    partitionTree = setPartitionTree;

    // now by convention the last node is the root
    // this would not be the natural result of rcb
    // however other algorithms may also need this swap step so
    // I made it general here
    int findRoot = -1;
    for(int n = 0; n < static_cast<int>(partitionTree.size()); ++n) {
      if(partitionTree[n].parent == 0) { // root has parent=0 by convention
        if(findRoot != -1) {
          throw std::logic_error("setPartitionTreeNodes found more than one root. "
            " This is not expected.");
        }
        findRoot = n;
      }
    }

    if(findRoot == -1) {
      throw std::logic_error( "setPartitionTreeNodes did not find a root. "
        " This is not expected." );
    }

    // potentially we would not do this and simply apply the 'swap' on the
    // data write out when requested - but that could create confusion
    // To do - verify how we want to apply this convention
    swapPartitionTreeNodes(findRoot, partitionTree.size()-1);

    // Swap the order of the two nodes under the root
    // This would intentionally put permPartNums out of sequence for rcb for ex.
    // Just for validating code right now...
    /*
    int newRoot = partitionTree.size()-1;
    int save1 = partitionTree[newRoot].children[1];
    partitionTree[newRoot].children[1] = partitionTree[newRoot].children[0];
    partitionTree[newRoot].children[0] = save1;
    */
  }

  //! \brief  swap to node indices and update all the corresponding parent and
  //          children indices to preserve the tree structure.
  //
  void swapPartitionTreeNodes(int a, int b) {
    if(a != b) {
      partitionTreeNode<part_t> saveOld = partitionTree[a];
      partitionTree[a] = partitionTree[b];
      partitionTree[b] = saveOld;

      // now we have to remap all parent/child indices
      for(int n = 0; n < static_cast<int>(partitionTree.size()); ++n) {
        partitionTreeNode<part_t> & node = partitionTree[n];
        if(node.parent == a+1) { // index+1 convention
          node.parent = b+1; // index+1 convention
        }
        else if(node.parent == b+1) { // index+1 convention
          node.parent = a+1; // index+1 convention
        }

        for(int c = 0; c < static_cast<int>(node.children.size()); ++c) {
          if(node.children[c] > 0 && node.children[c] == a+1) { // index+1 convention
            node.children[c] = b+1; // index+1 convention
          }
          else if(node.children[c] == b+1) { // index+1 convention
            node.children[c] = a+1; // index+1 convention
          }
        }
      }
    }
  }

  //! \brief  for partitioning methods, return bounding boxes of the 
  //          computed parts
  //          Not all partitioning algorithms will support
  //          this method.
  //
  virtual std::vector<coordinateModelPartBox<scalar_t, part_t> > &
  getPartBoxesView() const
  {
    Z2_THROW_NOT_IMPLEMENTED
  }

  //! \brief pointAssign method: Available only for some partitioning algorithms
  //          when a point lies on a part boundary, the lowest part
  //          number on that boundary is returned.
  //          Not all partitioning algorithms will support
  //          this method.
  //
  //   \param dim : the number of dimensions specified for the point in space
  //   \param point : the coordinates of the point in space; array of size dim
  //   \return the part number of a part overlapping the given point
  virtual part_t pointAssign(int dim, scalar_t *point) const
  {
    Z2_THROW_NOT_IMPLEMENTED
  }

  //! \brief boxAssign method:  Available only for some partitioning algorithms
  //   Return an array of all parts overlapping a given box in space.
  //   This method allocates memory for the return argument, but does not
  //   control that memory.  The user is responsible for freeing the 
  //   memory.
  //
  //   \param dim : (in) the number of dimensions specified for the box
  //   \param ptLower : (in) the coordinates of the lower corner of the box; 
  //                     array of size dim
  //   \param ptUpper : (in) the coordinates of the upper corner of the box; 
  //                     array of size dim
  //   \param nParts : (out) the number of parts overlapping the box
  //   \param parts :  (out) array of parts overlapping the box
  virtual void boxAssign(int dim, scalar_t *lower, scalar_t *upper,
                         size_t &nParts, part_t **partsFound) const
  {
    Z2_THROW_NOT_IMPLEMENTED
  }

  //! \brief returns serial communication graph of a computed partition
  //  Returned graph is identical on all processors, and represents the
  //  global communication pattern in the partition.
  //  
  //  \param comXAdj:  (out) the offset array:  offsets into comAdj
  //                         Format is standard CSR format:
  //                         # nbor parts of part i = comXAdj[i+1]-comXAdj[i]
  //                         That is, comXAdj[i] = Sum of # nbor parts of parts
  //                                               0 through i-1
  //  \param comAdj    (out) the neighboring parts
  virtual void getCommunicationGraph(
    const PartitioningSolution<Adapter> *solution,
    ArrayRCP<part_t> &comXAdj,
    ArrayRCP<part_t> &comAdj)
    // TODO:  Should the return args be ArrayViews?
  {
    Z2_THROW_NOT_IMPLEMENTED
  }

  //! \brief In mapping, returns the rank to which a part is assigned
  //  \param p: (in) the part for which the rank is sought
  //  This method need not be implemented by every algorithm or, indeed,
  //  for every mapping algorithm.  Mapping algorithms may provide this 
  //  function to prevent additional memory use in MappingSolution.
  //  For example, AlgContiguousMapping can compute this function implicitly, 
  //  with no additional storage.  However, Mapping algorithms can skip this
  //  function and, instead, register their results in MappingSolution.
  virtual int getRankForPart(part_t p)
  {
    Z2_THROW_NOT_IMPLEMENTED
  }

  //! \brief In mapping, returns a view of parts assigned to the current rank
  //  \param numParts: (out) the number of parts assigned to the current rank
  //  \param parts: (out) a view of the assigned parts
  //
  //  This method need not be implemented by every algorithm or, indeed,
  //  for every mapping algorithm.  Mapping algorithms may provide this 
  //  function to prevent additional memory use in MappingSolution.
  //  For example, AlgContiguousMapping can compute this function implicitly, 
  //  with no additional storage.  However, Mapping algorithms can skip this
  //  function and, instead, register their results in MappingSolution.
  virtual void getMyPartsView(part_t &numParts, part_t *&parts)
  {
    Z2_THROW_NOT_IMPLEMENTED
  }


private:
  std::vector<partitionTreeNode<part_t> > partitionTree;
};
  
}  //namespace Zoltan2
  
#endif
