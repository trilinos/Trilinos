// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __Panzer_EpetraVector_Write_GlobalEvaluationData_hpp__
#define   __Panzer_EpetraVector_Write_GlobalEvaluationData_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Epetra
#include "Epetra_Export.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

// Panzer
#include "Panzer_WriteVector_GlobalEvaluationData.hpp"

// Teuchos
#include "Teuchos_RCP.hpp"

// Thyra
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorBase.hpp"

namespace panzer
{

  /**
   *  \brief This class provides a boundary exchange communication mechanism
   *         for vectors.
   *
   *  \note This provides a "write" interface for parameters.
   */
  class EpetraVector_Write_GlobalEvaluationData
    :
    public WriteVector_GlobalEvaluationData
  {
    public:

      /**
       *  \brief Default Constructor.
       */
      EpetraVector_Write_GlobalEvaluationData()
        :
        isInitialized_(false)
      {
      } // end of Default Constructor

      /**
       *  \brief Copy Constructor.
       *
       *  \param[in] src The object to be copied.
       */
      EpetraVector_Write_GlobalEvaluationData(
        const EpetraVector_Write_GlobalEvaluationData& src)
        :
        isInitialized_(false)
      {
        initialize(src.exporter_, src.ghostedMap_, src.ownedMap_);
      } // end of Copy Constructor

      /**
       *  \brief Initializing Constructor.
       *
       *  Initialize this object with some Epetra communication objects. This
       *  method must be called before an object of this type can be used.
       *
       *  \param[in] exporter   Exporter for doing communication from the owned
       *                        to the ghosted vector.
       *  \param[in] ghostedMap Map describing the ghosted vector.
       *  \param[in] ownedMap   Map describing the owned vector.
       */
      EpetraVector_Write_GlobalEvaluationData(
        const Teuchos::RCP<const Epetra_Export>& exporter,
        const Teuchos::RCP<const Epetra_Map>&    ghostedMap,
        const Teuchos::RCP<const Epetra_Map>&    ownedMap)
        :
        isInitialized_(false)
      {
        initialize(exporter, ghostedMap, ownedMap);
      } // end of Initializing Constructor

      /**
       *  \brief Initialize this object with some `Epetra` communication
       *         objects.
       *
       *  This method must be called before an object of this type can be used.
       *
       *  \param[in] exporter   Exporter for doing communication from the owned
       *                        to the ghosted vector.
       *  \param[in] ghostedMap Map describing the ghosted vector.
       *  \param[in] ownedMap   Map describing the owned vector.
       */
      void
      initialize(
        const Teuchos::RCP<const Epetra_Export>& exporter,
        const Teuchos::RCP<const Epetra_Map>&    ghostedMap,
        const Teuchos::RCP<const Epetra_Map>&    ownedMap);

      /**
       *  \brief Communicate the ghosted data to the owned vector.
       *
       *  For this class, this method does the halo exchange for the vector.
       *
       *  \param[in] mem Not needed for this class, but part of the
       *                 `GlobalEvaluationData` interface.
       */
      virtual void
      ghostToGlobal(
        int mem);

      /**
       *  \brief Clear out the ghosted vector.                                   // JMG:  Is this right?
       */
      virtual void
      initializeData();

      /**
       *  \brief Determine if a Dirichlet adjustment is necessary.
       *
       *  \returns False.                                                        // JMG:  But why?
       */
      virtual bool
      requiresDirichletAdjustment() const
      {
        return false;
      } // end of requiresDirichletAdjustment()

      /**
       *  \brief Set the owned vector (`Epetra` version).
       *
       *  \param[in] ownedVector An `Epetra_Vector` that you would like to set
       *                         as the owned vector.
       */
      void
      setOwnedVector_Epetra(
        const Teuchos::RCP<Epetra_Vector>& ownedVector);

      /**
       *  \brief Get the ghosted vector (`Epetra` version).
       *
       *  \returns The ghosted vector as an `Epetra_Vector`.
       */
      Teuchos::RCP<Epetra_Vector>
      getGhostedVector_Epetra() const;

      /**
       *  \brief Set the owned vector (`Thyra` version).
       *
       *  \param[in] ownedVector A `Thyra::VectorBase<double>` that you would
       *                         like to set as the owned vector.
       */
      void
      setOwnedVector(
        const Teuchos::RCP<Thyra::VectorBase<double>>& ownedVector);

      /**
       *  \brief Get the owned vector (`Thyra` version).
       *
       *  \returns The owned vector as a `Thyra::VectorBase<double>`.
       */
      Teuchos::RCP<Thyra::VectorBase<double>>
      getOwnedVector() const;

      /**
       *  \brief Get the ghosted vector (`Thyra` version).
       *
       *  \returns The ghosted vector as a `Thyra::VectorBase<double>`.
       */
      Teuchos::RCP<Thyra::VectorBase<double>>
      getGhostedVector() const;

      /**
       *  \brief Is this object initialized?
       *
       *  \returns Whether or not this object is initialized.
       */
      virtual bool
      isInitialized() const
      {
        return isInitialized_;
      } // end of isInitialized()

      /**
       *  \brief Print the object.
       *
       *  This is a diagnostic function for debugging purposes.
       *
       *  \param[in,out] os The output stream to which the data should be
       *                    printed.
       */
      void
      print(
        std::ostream& os) const;

    private:

      /**
       *  \brief A flag indicating whether or not the object has been
       *         initialized.
       */
      bool isInitialized_;

      /**
       *  \brief The map corresponding to the ghosted vector.
       */
      Teuchos::RCP<const Epetra_Map> ghostedMap_;

      /**
       *  \brief The map corresponding to the owned vector.
       */
      Teuchos::RCP<const Epetra_Map> ownedMap_;

      /**
       *  \brief The vector space corresponding to the ghosted vector.
       */
      Teuchos::RCP<const Thyra::VectorSpaceBase<double>> ghostedSpace_;

      /**
       *  \brief The vector space corresponding to the owned vector.
       */
      Teuchos::RCP<const Thyra::VectorSpaceBase<double>> ownedSpace_;

      /**
       *  \brief The exporter used to communicate between the owned and ghosted
       *         vectors.
       */
      Teuchos::RCP<const Epetra_Export> exporter_;

      /**
       *  \brief The ghosted vector.
       */
      Teuchos::RCP<Epetra_Vector> ghostedVector_;

      /**
       *  \brief The owned vector.
       */
      Teuchos::RCP<Thyra::VectorBase<double>> ownedVector_;

  }; // end of class EpetraVector_Write_GlobalEvaluationData

} // end of namespace panzer

#endif // __Panzer_EpetraVector_Write_GlobalEvaluationData_hpp__