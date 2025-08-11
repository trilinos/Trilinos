// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __Panzer_ReadOnlyVector_GlobalEvaluationData_hpp__
#define   __Panzer_ReadOnlyVector_GlobalEvaluationData_hpp__

// Panzer
#include "PanzerDiscFE_config.hpp"

#include "Panzer_GlobalEvaluationData.hpp"

// Teuchos
#include "Teuchos_RCP.hpp"

// Thyra
#include "Thyra_VectorBase.hpp"

namespace panzer {

/** This class encapsulates the needs of a gather operation to do a halo
  * exchange.
  */
class ReadOnlyVector_GlobalEvaluationData : public GlobalEvaluationData {
public:

  //! Virtual d
  virtual ~ReadOnlyVector_GlobalEvaluationData() {}

  //! Is this object initialized
  virtual bool isInitialized() const = 0;

  /** For this class, this method does the halo exchange for the
    * vector.
    */
  virtual void globalToGhost(int mem) = 0;

  /** For this class, this method does nothing.
    */
  virtual void ghostToGlobal(int /* mem */) {}

  //! Set the owned vector
  virtual void setOwnedVector(const Teuchos::RCP<const Thyra::VectorBase<double> > & ownedVector) = 0;

  //! Get the owned vector
  virtual Teuchos::RCP<const Thyra::VectorBase<double> > getOwnedVector() const = 0;

  //! Get the ghosted vector
  virtual Teuchos::RCP<Thyra::VectorBase<double> > getGhostedVector() const = 0;

}; // end of class ReadOnlyVector_GlobalEvaluationData

} // end of namespace panzer

#endif // __Panzer_ReadOnlyVector_GlobalEvaluationData_hpp__
