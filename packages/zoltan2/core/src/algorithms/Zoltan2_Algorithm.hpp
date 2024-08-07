// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include <Zoltan2_MatrixPartitioningSolution.hpp>
#include <Zoltan2_MappingSolution.hpp>
#include <Zoltan2_CoordinatePartitioningGraph.hpp>


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
  virtual int localOrder(const RCP<LocalOrderingSolution<lno_t> > &/* solution */)
  {
    Z2_THROW_NOT_IMPLEMENTED
  }

  //! \brief Ordering method
  virtual int globalOrder(const RCP<GlobalOrderingSolution<gno_t> > &/* solution */)
  {
    Z2_THROW_NOT_IMPLEMENTED 
  }
  
  //! \brief Coloring method
  virtual void color(const RCP<ColoringSolution<Adapter> > &/* solution */) 
  {
    Z2_THROW_NOT_IMPLEMENTED
  }
  
  //! \brief Matching method
  virtual void match() { 
    Z2_THROW_NOT_IMPLEMENTED 
  }

  //! \brief Partitioning method
  virtual void partition(const RCP<PartitioningSolution<Adapter> > &/* solution */) 
  {
    Z2_THROW_NOT_IMPLEMENTED
  }

  //! \brief Matrix Partitioning method
  virtual void partitionMatrix(const RCP<MatrixPartitioningSolution<Adapter> > &/* solution */) 
  {
    Z2_THROW_NOT_IMPLEMENTED
  }

  //! \brief Mapping method
  virtual void map(const RCP<MappingSolution<Adapter> > &/* solution */) 
  {
    Z2_THROW_NOT_IMPLEMENTED
  }

  //! \brief  return if algorithm determins tree to be binary
  virtual bool isPartitioningTreeBinary() const
  {
    Z2_THROW_NOT_IMPLEMENTED
  }

  //! \brief  for partitioning methods, fill arrays with partition tree info
  virtual void getPartitionTree(part_t /* numParts */,
                        part_t & /* numTreeVerts */,
                        std::vector<part_t> & /* permPartNums */,
                        std::vector<part_t> & /* splitRangeBeg */,
                        std::vector<part_t> & /* splitRangeEnd */,
                        std::vector<part_t> & /* treeVertParents */) const
  {
    Z2_THROW_NOT_IMPLEMENTED
  }

  //! \brief  for partitioning methods, return bounding boxes of the 
  //          computed parts
  //          Not all partitioning algorithms will support
  //          this method.
  //
  virtual std::vector<coordinateModelPartBox> &
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
  virtual part_t pointAssign(int /* dim */, scalar_t * /* point */) const
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
  virtual void boxAssign(int /* dim */, scalar_t * /* lower */, scalar_t * /* upper */,
                         size_t &/* nParts */, part_t ** /* partsFound */) const
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
    const PartitioningSolution<Adapter> * /* solution */,
    ArrayRCP<part_t> &/* comXAdj */,
    ArrayRCP<part_t> &/* comAdj */)
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
  virtual int getRankForPart(part_t /* p */)
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
  virtual void getMyPartsView(part_t &/* numParts */, part_t *&/* parts */)
  {
    Z2_THROW_NOT_IMPLEMENTED
  }


private:
};
  
}  //namespace Zoltan2
  
#endif
