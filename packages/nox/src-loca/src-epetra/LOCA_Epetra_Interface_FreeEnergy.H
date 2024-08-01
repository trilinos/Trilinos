// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_EPETRA_INTERFACE_FREEENERGY_H
#define LOCA_EPETRA_INTERFACE_FREEENERGY_H

#include "NOX_Common.H"

// Forward declarations
class Epetra_Vector;

namespace LOCA {
  namespace Epetra {
    namespace Interface {

      /*!
       * \brief Used by LOCA::Epetra::Group to provide a link to the
       * external code for computing the free energy.
       */

      class FreeEnergy {

      public:

    //! Constructor
    FreeEnergy() {};

    //! Destructor
    virtual ~FreeEnergy() {};

    /*!
     * \brief Call user routine for computing the FreeEnergy
         * of a system, for use in Phase Transition Tracking alg
     */

    virtual double
        computeFreeEnergy(const Epetra_Vector& x) = 0;

      };
    } // namespace Interface
  } // namespace Epetra
} // namespace LOCA

#endif
