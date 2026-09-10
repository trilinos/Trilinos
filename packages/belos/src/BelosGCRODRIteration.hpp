// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_GCRODR_ITERATION_HPP
#define BELOS_GCRODR_ITERATION_HPP

/*! \file BelosGCRODRIteration.hpp
    \brief Common interface for GCRODR-like iteration classes.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosIteration.hpp"

namespace Belos {

template <class ScalarType, class MV, class DM>
struct GCRODRIterState;

/*! \class Belos::GCRODRIteration
    \brief Common base interface for GCRODRIter and FGCRODRIter.

    This mirrors the role of GmresIteration for BlockGmresSolMgr.  It
    lets GCRODRSolMgr store either the standard or flexible iterator in
    the same RCP variable.
*/
template<class ScalarType, class MV, class OP, class DM = DefaultDenseMatrix<int, ScalarType> >
class GCRODRIteration : virtual public Iteration<ScalarType,MV,OP,DM> {
public:
  virtual ~GCRODRIteration() {}

  using Iteration<ScalarType,MV,OP,DM>::initialize;
  virtual void initialize(GCRODRIterState<ScalarType,MV,DM>& newstate) = 0;

  virtual GCRODRIterState<ScalarType,MV,DM> getState() const = 0;

  virtual void updateLSQR(int dim = -1) = 0;

  virtual int getCurSubspaceDim() const = 0;

  virtual int getMaxSubspaceDim() const = 0;

  virtual void setSize(int recycledBlocks, int numBlocks) = 0;
};

} // namespace Belos

#endif // BELOS_GCRODR_ITERATION_HPP
