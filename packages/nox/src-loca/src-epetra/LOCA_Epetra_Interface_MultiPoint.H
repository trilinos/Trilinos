// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_EPETRA_INTERFACE_MULTIPOINT_H
#define LOCA_EPETRA_INTERFACE_MULTIPOINT_H

#include "LOCA.H"
#include "NOX_Common.H"
#include "LOCA_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"

#ifdef HAVE_NOX_EPETRAEXT

#include <iostream>
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_DistObject.h"
#include "Epetra_Import.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "EpetraExt_MultiMpiComm.h"
#else
#include "Epetra_SerialComm.h"
#include "EpetraExt_MultiSerialComm.h"
#endif

#include "EpetraExt_BlockVector.h"
#include "EpetraExt_BlockCrsMatrix.h"

// Forward declarations
namespace LOCA {
  namespace Epetra {
    namespace Interface {
      class Required;
    }
  }
}
namespace NOX {
  namespace Epetra {
    namespace Interface {
      class Jacobian;
    }
  }
}


namespace LOCA {
  namespace Epetra {
    namespace Interface {

      /*!
       * \brief Used by LOCA::Epetra::Group to provide a link to the
       * external code for creating and solving space-time formulations
       * of the user's problem.
       *
       * Code that takes a standard NOX-LOCA problem interface for
       * spatially-discretized problems, and creates a NOX-LOCA interface
       * for space-time problems.  This interface is called by multiple
       * replicas of the spatially discretized problem for parallelism of
       * the time domain. Each replica can also own one or more time steps.
       * This interface assembles a single space-time system including the
       * solution and residual BlockVectors and the Jacobian BlockCrsMatrix.
       * This interface relies heavily on the EpetraExt_Block* classes.
       *
       */
      class MultiPoint :
    public LOCA::Epetra::Interface::Required,
    public NOX::Epetra::Interface::Jacobian {

      public:

    /*!
      \brief Constructor

      The interface, vector, and matrix arguments are all for the spatially
      decomposed problem. The MultiMpiCOmm communicator gives this
      class the information needed to construct the space-time system.
    */
    MultiPoint(
             const Teuchos::RCP<LOCA::Epetra::Interface::Required> &iReq,
             const Teuchos::RCP<NOX::Epetra::Interface::Jacobian> &iJac,
         const Epetra_MultiVector &splitMultiVec,
         const Teuchos::RCP<Epetra_RowMatrix> &splitJac,
         const Teuchos::RCP<EpetraExt::MultiComm> &globalComm);

    //! Destructor
    virtual ~MultiPoint();

    /*!
      \brief NOX Interface function for computing the residual.

      Computes \em F given the vector \em x.
      Returns true if computation was successful.
    */
    virtual bool computeF(const Epetra_Vector& x, Epetra_Vector& F,
                  const FillType fillFlag);

    /*!
      \brief NOX Interface function for computing the Jacobian matrix.

      Computes \em Jac given the vector \em x.
      Returns true if computation was successful.
    */
    virtual bool computeJacobian(const Epetra_Vector& x,
                     Epetra_Operator& Jac);

    //! LOCA Interface function for setting the continuation parameter
    virtual void setParameters(const LOCA::ParameterVector& param);

    //! LOCA Interface function for printing the solution
    virtual void printSolution(const Epetra_Vector& x_, double conParam);

    //! Accessor for the global solution BlockVector
    virtual EpetraExt::BlockVector& getSolution();

    //! Accessor for the global BlockCrsMatrix
    virtual EpetraExt::BlockCrsMatrix& getJacobian();

    //! LOCA Interface exception handler
    virtual void throwError(const std::string& functionName,
                const std::string& errorMsg) const;

      private:

    //! Pointer to interface for usual time-dependent problem
    Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq;
    Teuchos::RCP< NOX::Epetra::Interface::Jacobian> iJac;

    //! Pointer to split (spatial) Jacobian matrix.
    Teuchos::RCP<Epetra_RowMatrix> splitJac;

    //! Pointer to the global (full XYZT) communicator.
    Teuchos::RCP<EpetraExt::MultiComm> globalComm;

    //! Split (spatial) input vector.
    Epetra_Vector splitVec;

    //! Split (spatial) residual vector.
    Epetra_Vector splitRes;

    //! Pointer to global XYZT Jacobian matrix
    EpetraExt::BlockCrsMatrix* jacobian;

    //! Pointer to global XYZT solution vector.
    EpetraExt::BlockVector* solution;

    //! Pointer to global XYZT solution overlap vector.
    EpetraExt::BlockVector* solutionOverlap;

    //! Pointer to global XYZT solution overlap vector importer.
    Epetra_Import* overlapImporter;

    //! Number of time steps computed on each time domain.
    int timeStepsOnTimeDomain;

    //! Total number of time step domains.
    int numTimeDomains;

    //! Time domain on current processor.
    int timeDomain;

    //! Continuation step number (used in printing)
    int conStep;

    /*!
      \brief Stencil for each row of global XYZT Jacobian matrix.

      Used in creating global XYZT Jacobian matrix for different
      finite difference schemes.
    */
    std::vector< std::vector<int> >* rowStencil;

    //! Set of indices into global XYZT Jacobian matrix.
    std::vector<int>* rowIndex;

    //! Pointer to split (spatial) Jacobian matrix in CRS matrix form.
    Epetra_CrsMatrix* splitJacCrs;

      };

    } // namespace Interface
  } // namespace Epetra
} // namespace LOCA

#endif // HAVE_NOX_EPETRAEXT
#endif // LOCA_EPETRA_INTERFACE_MULTIPOINT_H
