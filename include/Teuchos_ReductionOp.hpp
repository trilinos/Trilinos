// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_REDUCTION_OP_HPP
#define TEUCHOS_REDUCTION_OP_HPP

#include "Teuchos_Describable.hpp"


namespace Teuchos {


/** \brief Base interface class for user-defined reduction operations for
 * objects that use value semantics.
 *
 * ToDo: Finish documentation!
 *
 * \note Ordinal refers to the template parameter of \c Comm, not to
 *   the type of array indices.
 */
template<typename Ordinal, typename T>
class ValueTypeReductionOp : public Describable {
public:
  /** \brief . */
  virtual void reduce(
    const Ordinal count,
    const T inBuffer[],
    T inoutBuffer[]
    ) const = 0;
};


/** \brief Base interface class for user-defined reduction operations for
 * objects that use reference semantics.
 *
 * ToDo: Finish documentation!
 *
 * \note Ordinal refers to the template parameter of \c Comm, not to
 *   the type of array indices.
 */
template<typename Ordinal, typename T>
class ReferenceTypeReductionOp : public Describable {
public:
  /** \brief . */
  virtual void reduce(
    const Ordinal count,
    const T*const inBuffer[],
    T*const inoutBuffer[]
    ) const = 0;
};


} // namespace Teuchos


#endif // TEUCHOS_REDUCTION_OP_HPP
