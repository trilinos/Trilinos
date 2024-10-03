// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __Panzer_EpetraVector_ReadOnly_GlobalEvaluationData_hpp__
#define   __Panzer_EpetraVector_ReadOnly_GlobalEvaluationData_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Epetra
#include "Epetra_Import.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

// Panzer
#include "Panzer_ReadOnlyVector_GlobalEvaluationData.hpp"

// Teuchos
#include "Teuchos_RCP.hpp"

///////////////////////////////////////////////////////////////////////////////
//
//  Forward Declarations
//
///////////////////////////////////////////////////////////////////////////////

namespace panzer
{
  /**
   *  \brief This class provides a boundary exchange communication mechanism
   *         for vectors.
   *
   *  \note This provides a read-only (RO) interface for parameters (so vectors
   *        are write protected).
   */
  class EpetraVector_ReadOnly_GlobalEvaluationData
    :
    public ReadOnlyVector_GlobalEvaluationData
  {
    public:

      /**
       *  \brief Default Constructor.
       */
      EpetraVector_ReadOnly_GlobalEvaluationData()
        :
        isInitialized_(false)
      {
      } // end of Default Constructor

      /**
       *  \brief Copy Constructor.
       *
       *  \param[in] src The object to be copied.
       */
      EpetraVector_ReadOnly_GlobalEvaluationData(
        const EpetraVector_ReadOnly_GlobalEvaluationData& src)
        :
        isInitialized_(false)
      {
        initialize(src.importer_, src.ghostedMap_, src.ownedMap_);
      } // end of Copy Constructor

      /**
       *  \brief Initializing Constructor.
       *
       *  Initialize this object with some `Epetra` communication objects.
       *  This method must be called before an object of this type can be used.
       *
       *  \param[in] importer   Importer for doing communication from the owned
       *                        to the ghosted vector.
       *  \param[in] ghostedMap Map describing the ghosted vector.
       *  \param[in] ownedMap   Map describing the owned vector.
       */
      EpetraVector_ReadOnly_GlobalEvaluationData(
        const Teuchos::RCP<const Epetra_Import>& importer,
        const Teuchos::RCP<const Epetra_Map>&    ghostedMap,
        const Teuchos::RCP<const Epetra_Map>&    ownedMap)
        :
        isInitialized_(false)
      {
        initialize(importer, ghostedMap, ownedMap);
      } // end of Initializing Constructor

      /**
       *  \brief Choose a few GIDs and, instead of zeroing them out in the
       *         ghosted vector, set them to a specified value.
       *
       *  \note This is only useful for GIDs in the ghosted map.
       *  \note This must be called before `initialize()`.
       *  \note No attempt to synchronize these values across a processor is
       *        made, so it's up to the user to be consistent.
       *
       *  \param[in] indices A `std::vector` of global IDs.
       *  \param[in] value   The value to be assigned to the entries given by
       *                     indices.
       */
      void
      useConstantValues(
        const std::vector<int>& indices,
        double                  value);

      /**
       *  \brief Initialize this object with some `Epetra` communication
       *         objects.
       *
       *  This method must be called before an object of this type can be used.
       *
       *  \param[in] importer   Importer for doing communication from the owned
       *                        to the ghosted vector.
       *  \param[in] ghostedMap Map describing the ghosted vector.
       *  \param[in] ownedMap   Map describing the owned vector.
       */
      void
      initialize(
        const Teuchos::RCP<const Epetra_Import>& importer,
        const Teuchos::RCP<const Epetra_Map>&    ghostedMap,
        const Teuchos::RCP<const Epetra_Map>&    ownedMap);

      /**
       *  \brief Communicate the owned data to the ghosted vector.
       *
       *  For this class, this method does the halo exchange for the vector.
       *
       *  \param[in] mem Not needed for this class, but part of the
       *                 `GlobalEvaluationData` interface.
       */
      virtual void
      globalToGhost(
        int mem = 0);

      /**
       *  \brief Clear out the ghosted vector.
       */
      virtual void
      initializeData();

      /**
       *  \brief Communicate the ghosted data to the owned vector.
       *
       *  For this class this method does nothing.
       *
       *  \param[in] mem Not needed for this class, but part of the
       *                 `GlobalEvaluationData` interface.
       */
      virtual void
      ghostToGlobal(
        int mem = 0);

      /**
       *  \brief Determine if a Dirichlet adjustment is necessary.
       *
       *  For this class, there's nothing to do because it's read-only.
       *
       *  \returns False.
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
        const Teuchos::RCP<const Epetra_Vector>& ownedVector);

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
        const Teuchos::RCP<const Thyra::VectorBase<double>>& ownedVector);

      /**
       *  \brief Get the owned vector (`Thyra` version).
       *
       *  \returns The owned vector as a `const Thyra::VectorBase<double>`.
       */
      Teuchos::RCP<const Thyra::VectorBase<double>>
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
       *  \returns Whether or not the object is initialized.
       */
      bool
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
       *  \brief The importer used to communicate between the owned and ghosted
       *         vectors.
       */
      Teuchos::RCP<const Epetra_Import> importer_;

      /**
       *  \brief The ghosted vector.
       */
      Teuchos::RCP<Epetra_Vector> ghostedVector_;

      /**
       *  \brief The owned vector.
       */
      Teuchos::RCP<const Thyra::VectorBase<double>> ownedVector_;

      /**
       *  \brief A list of global IDs (which will be translated to local IDs),
       *         paired with a value to be assigned in the `ghostedVector_`.
       */
      typedef std::pair<std::vector<int>, double> FilteredPair;

      /**
       *  \brief The list of filtered pairs, used to initialize values on the
       *         `ghostedVector_`.
       */
      std::vector<FilteredPair> filteredPairs_;

  }; // end of class EpetraVector_ReadOnly_GlobalEvaluationData

} // end of namespace panzer

#endif // __Panzer_EpetraVector_ReadOnly_GlobalEvaluationData_hpp__