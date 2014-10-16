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

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_ColoringSolution.hpp>
#include <Zoltan2_OrderingSolution.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

#define Z2_THROW_NOT_IMPLEMENTED_IN_ALGORITHM \
  { \
    std::ostringstream emsg; \
    emsg << __FILE__ << "," << __LINE__ \
         << " error:  " <<  __func__zoltan2__ \
         << " is not implement in selected algorithm " \
         << std::endl; \
    throw std::runtime_error(emsg.str()); \
  }
  

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
  typedef typename Adapter::zgid_t zgid_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::part_t part_t;

  // Virtual destructor needed to avoid undefined behavior and compiler warnings
  virtual ~Algorithm() {}

  //! \brief Ordering method
  virtual int order(const RCP<OrderingSolution<zgid_t, lno_t> > &solution) 
  {
    Z2_THROW_NOT_IMPLEMENTED_IN_ALGORITHM 
  }
  
  //! \brief Coloring method
  virtual void color(const RCP<ColoringSolution<Adapter> > &solution) 
  {
    Z2_THROW_NOT_IMPLEMENTED_IN_ALGORITHM
  }
  
  //! \brief Matching method
  virtual void match() { Z2_THROW_NOT_IMPLEMENTED_IN_ALGORITHM }

  //! \brief Partitioning method
  virtual void partition(const RCP<PartitioningSolution<Adapter> > &solution) 
  {
    Z2_THROW_NOT_IMPLEMENTED_IN_ALGORITHM
  }

  //! \brief pointAssign method: Available only for some partitioning algorithms
  //          when a point lies on a part boundary, the lowest part
  //          number on that boundary is returned.
  //          Note that not all partitioning algorithms will support
  //          this method.
  //
  //   \param dim : the number of dimensions specified for the point in space
  //   \param point : the coordinates of the point in space; array of size dim
  //   \return the part number of a part overlapping the given point
  virtual part_t pointAssign(int dim, scalar_t *point) const
  {
    Z2_THROW_NOT_IMPLEMENTED_IN_ALGORITHM
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
    Z2_THROW_NOT_IMPLEMENTED_IN_ALGORITHM
  }

private:
};
  
}  //namespace Zoltan2
  
#endif
