// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_DiagonalPreconditionerFactory_hpp__
#define __Teko_DiagonalPreconditionerFactory_hpp__

// Teko includes
#include "Teko_PreconditionerState.hpp"
#include "Teko_PreconditionerFactory.hpp"

class EpetraExt_PointToBlockDiagPermute;

namespace Teko {

using Thyra::DefaultPreconditioner;
using Thyra::LinearOpBase;

/** Preconditioning factories must supply a 'State' class which
 * is where data specific to the preconditioner construction
 * is stored. The constructor will be invoked elsewhere.
 */
class DiagonalPrecondState : public Teko::PreconditionerState {
 public:
  DiagonalPrecondState();

  Teuchos::RCP<EpetraExt_PointToBlockDiagPermute> BDP_;
};

/** \brief Preconditioner factory for building explcit inverse of diagonal operators.
  *        This includes block operators.
  *
  * Preconditioner factory for building explicit inverse diagonal operators, including
  * block diagonals. These operators need to be Epetra_CrsMatrices under the hood or
  * this will bomb. Uses EpetraExt_PointToBlockDiagPermute.
  *
  * To invoke this preconditioner using the XML file a diagonal inverse
  * needs to be specified. For example the following XML code creates
  * a AbsRowSum preconditioner
  *
    \verbatim
    <ParameterList name="AbsRowSum">
       <Parameter name="Type" type="string" value="Explicit Diagonal Preconditioner"/>
       <Parameter name="Diagonal Type" type="string" value="AbsRowSum"/> <!-- Could also be "Lumped"
  or "Diagonal" -->
    </ParameterList>
    \endverbatim
  *
  * The block diagonal inverse operator is constructed from the following XML code
  *
    \verbatim
    <ParameterList name="4-Block Diagonal">
       <Parameter name="Type" type="string" value="Explicit Diagonal Preconditioner"/>
       <Parameter name="Diagonal Type" type="string" value="BlkDiag"/>
       <Parameter name="contiguous block size" type="int" value="4"/> <!-- block size of 4 -->
    </ParameterList>
    \endverbatim
  */
class DiagonalPreconditionerFactory : public virtual Teko::PreconditionerFactory {
 public:
  //! @name Constructors.
  //@{

  DiagonalPreconditionerFactory();

  //@}

  //! Builds a preconditioner state object
  Teuchos::RCP<PreconditionerState> buildPreconditionerState() const;

  /** Create the diagonal preconditioner operator.
   */
  LinearOp buildPreconditionerOperator(LinearOp& lo, PreconditionerState& state) const;

  //! Initialize from a parameter list
  virtual void initializeFromParameterList(const Teuchos::ParameterList& pl);

 protected:
  //! some members
  mutable Teuchos::ParameterList List_;

  DiagonalType diagonalType_;
};

}  // end namespace Teko

#endif
