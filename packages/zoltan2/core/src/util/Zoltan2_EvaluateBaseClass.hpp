// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_EvaluateBaseClass.hpp
 *  \brief Base class for the EvaluatePartition and EvaluateOrdering classes.
 */

#ifndef ZOLTAN2_EVALUATEBASECLASS_HPP
#define ZOLTAN2_EVALUATEBASECLASS_HPP

namespace Zoltan2{

/*! \brief A base class for EvaluatePartition, EvaluateOrdering, ...
 */
class EvaluateBaseClassRoot{
  public:
    virtual ~EvaluateBaseClassRoot() {} // required virtual declaration

    /*! \brief Print all metrics */
    virtual void printMetrics(std::ostream &/* os */) const {};
};

template <typename Adapter>
class EvaluateBaseClass : public EvaluateBaseClassRoot {
public:
};

}   // namespace Zoltan2

#endif
