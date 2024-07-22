// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_BLOCKED_LINEAR_OP_WITH_SOLVE_BASE_HPP
#define THYRA_BLOCKED_LINEAR_OP_WITH_SOLVE_BASE_HPP

#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"


namespace Thyra {


/** \brief Base interface for linear operators with a solve that can be
 * accessed as sub-blocks.
 *
 * \ingroup Thyra_Op_Solve_extended_interfaces_code_grp
 *
 * ToDo: Finish Documentation.
 */
template<class Scalar>
class BlockedLinearOpWithSolveBase
  : virtual public LinearOpWithSolveBase<Scalar>,
    virtual public BlockedLinearOpBase<Scalar>
  
{
public:

  /** \brief . */
  virtual Teuchos::RCP<LinearOpWithSolveBase<Scalar> >
  getNonconstLOWSBlock(const int i, const int j) = 0; 

  /** \brief . */
  virtual Teuchos::RCP<const LinearOpWithSolveBase<Scalar> >
  getLOWSBlock(const int i, const int j) const = 0; 

private:
  
  // Not defined and not to be called
  BlockedLinearOpWithSolveBase<Scalar>&
  operator=(const BlockedLinearOpWithSolveBase<Scalar>&);

};


} // namespace Thyra


#endif // THYRA_BLOCKED_LINEAR_OP_WITH_SOLVE_BASE_HPP
