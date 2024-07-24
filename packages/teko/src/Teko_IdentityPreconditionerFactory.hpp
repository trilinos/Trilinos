// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_IdentityPreconditionerFactory_hpp__
#define __Teko_IdentityPreconditionerFactory_hpp__

#include "Teuchos_RCP.hpp"

#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_BlockInvDiagonalStrategy.hpp"

namespace Teko {

/** Uses a scaled identity operator as the approximate inverse.
  * This is a last resort!
  *
  \code
       <Parameter name="Type" type="string" value="Identity"/>
       <Parameter name="Scaling" type="double" value="<Some value>"/>
  \endcode
  */
class IdentityPreconditionerFactory : public PreconditionerFactory {
 public:
  //! @name Constructors.
  //@{

  /** Build an empty Identity preconditioner factory.
   */
  IdentityPreconditionerFactory();

  //@}

  /** \brief Create the Identity preconditioner operator.
   *
   * This method breaks apart the BlockLinearOp and builds a block
   * diagonal preconditioner. The inverse of the diagonals are specified
   * by the BlockInvDiagonalStrategy object.
   */
  LinearOp buildPreconditionerOperator(LinearOp& lo, PreconditionerState& state) const;

 protected:
  //! some members
  double scaling_;

  //! Initialize from a parameter list
  virtual void initializeFromParameterList(const Teuchos::ParameterList& pl);
};

}  // end namespace Teko

#endif
