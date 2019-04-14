// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_AuxModel_hpp
#define Tempus_AuxModel_hpp


#include "Tempus_ODEParameters.hpp"


namespace Tempus {


/** \brief A model to compute auxiliary variables.
 *
 *  This could grow into the Tempus version of ModelEvaluator.
 */
template <typename Scalar>
class AuxModel
{
public:

  virtual void evaluateAuxModel(
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & y,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & yDot,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & x,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & xDot,
  const Scalar time,
  const Teuchos::RCP<ODEParameters<Scalar> > & p ) = 0;

};

} // namespace Tempus

#endif // Tempus_AuxModel_hpp
