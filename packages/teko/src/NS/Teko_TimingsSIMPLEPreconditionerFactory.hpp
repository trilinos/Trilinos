// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_TimingsSIMPLEPreconditionerFactory_hpp__
#define __Teko_TimingsSIMPLEPreconditionerFactory_hpp__

#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_DiagnosticLinearOp.hpp"

#include "Teuchos_Time.hpp"

namespace Teko {
namespace NS {

// Declaration of the preconditioner factory
/** The basic XML parameter list for SIMPLE looks like.

   \code
    <ParameterList name="SIMPLEC">
      <Parameter name="Type" type="string" value="NS SIMPLE-Timed"/>

      <!-- Inverse operations to use -->
      <Parameter name="Inverse Velocity Type" type="string" value="ML_GS-Velocity"/>
      <Parameter name="Inverse Pressure Type" type="string" value="ML_GS-Pressure"/>

      <Parameter name="Explicit Velocity Inverse Type" type="string" value="AbsRowSum"/>

      <Parameter name="Alpha" type="double" value="0.9"/>
    </ParameterList>
  \endcode
  */
class TimingsSIMPLEPreconditionerFactory : public BlockPreconditionerFactory {
 public:
  // Constructor
  TimingsSIMPLEPreconditionerFactory(const Teuchos::RCP<InverseFactory>& inverse, double alpha);

  // Constructor
  TimingsSIMPLEPreconditionerFactory(const Teuchos::RCP<InverseFactory>& invVelFactory,
                                     const Teuchos::RCP<InverseFactory>& invPrsFactory,
                                     double alpha);

  //! Default constructor
  TimingsSIMPLEPreconditionerFactory();

  //! Destructor that outputs construction timings
  virtual ~TimingsSIMPLEPreconditionerFactory();

  // Function inherited from BlockPreconditionerFactory
  LinearOp buildPreconditionerOperator(BlockedLinearOp& blo, BlockPreconditionerState& state) const;

  //! Set the mass matrix for this factory
  virtual void setMassMatrix(Teko::LinearOp& mass) { massMatrix_ = mass; }

  //! For assisting in construction of the preconditioner
  virtual Teuchos::RCP<Teuchos::ParameterList> getRequestedParameters() const;

  //! For assisting in construction of the preconditioner
  virtual bool updateRequestedParameters(const Teuchos::ParameterList& pl);

  //! Initialize from a parameter list
  virtual void initializeFromParameterList(const Teuchos::ParameterList& pl);

 protected:
  using Teko::BlockPreconditionerFactory::buildPreconditionerOperator;

  // class members
  Teuchos::RCP<InverseFactory> customHFactory_;
  Teuchos::RCP<InverseFactory> invVelFactory_;
  Teuchos::RCP<InverseFactory> invPrsFactory_;
  double alpha_;
  DiagonalType fInverseType_;
  // enum FInverseType {Diagonal,Lumped,AbsRowSum,Custom} fInverseType_;

  bool useMass_;
  Teko::LinearOp massMatrix_;

  // Info for the block-diagonal approximation to H if used.
  mutable Teuchos::ParameterList BlkDiagList_;

  mutable Teuchos::Time constrTotal_;
  mutable Teuchos::Time subTotal_;
  mutable int constrCount_;

  mutable Teuchos::RCP<DiagnosticLinearOp> timed_HBt_;
  mutable Teuchos::RCP<DiagnosticLinearOp> timed_B_;
  mutable Teuchos::RCP<DiagnosticLinearOp> timed_invF_;
  mutable Teuchos::RCP<DiagnosticLinearOp> timed_invS_;
  mutable Teuchos::RCP<DiagnosticLinearOp> timed_iU_t_iL_;
};

}  // end namespace NS
}  // end namespace Teko

#endif
