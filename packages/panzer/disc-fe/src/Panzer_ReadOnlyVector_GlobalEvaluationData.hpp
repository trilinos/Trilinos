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
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Panzer_KokkosUtils_VectorToView.hpp"
#endif

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

#ifdef PANZER_HAVE_EPETRA_STACK
  /**
   *  \brief Element access.
   *
   *  Get the `lid`-th element in this `GlobalEvaluationData`.
   *
   *  \note This will pull the appropriate element out of either the owned or
   *        ghosted vector, depending on the value of `lid`.
   *
   *  \param[in] lid The local ID of the element you'd like to get.
   *
   *  \returns The `lid`-th element in this `GlobalEvaluationData`.
   */
  const double&
  operator[](
    const int& lid) const
  {
    if (lid < static_cast<int>(ownedView_.extent(0)))
      return ownedView_(lid);
    else // if (lid >= static_cast<int>(ownedView_.extent(0)))
      return ghostedView_(lid - ownedView_.extent(0));
  } // end of operator[]()

protected:

  /**
   *  \brief The `PHX::View` of the owned vector.
   */
  typename panzer::kokkos_utils::VectorToViewTraits<const Epetra_Vector>::View
  ownedView_;

  /**
   *  \brief The `PHX::View` of the ghosted vector.
   */
  typename panzer::kokkos_utils::VectorToViewTraits<Epetra_Vector>::View
  ghostedView_;
#endif

}; // end of class ReadOnlyVector_GlobalEvaluationData

} // end of namespace panzer

#endif // __Panzer_ReadOnlyVector_GlobalEvaluationData_hpp__
