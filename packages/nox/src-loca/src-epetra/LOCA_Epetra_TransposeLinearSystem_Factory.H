// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_EPETRA_TRANSPOSELINEARSYSTEM_FACTORY_H
#define LOCA_EPETRA_TRANSPOSELINEARSYSTEM_FACTORY_H

#include "Teuchos_RCP.hpp"

// forward declarations
class Epetra_Operator;
namespace Teuchos {
  class ParameterList;
}
namespace NOX {
  namespace Epetra {
    class LinearSystem;
  }
}
namespace LOCA {
  class GlobalData;
  namespace Epetra {
    namespace TransposeLinearSystem {
      class AbstractStrategy;
    }
  }
}

namespace LOCA {

  namespace Epetra {

    namespace TransposeLinearSystem {

      //! Factory for creating transpose linear system strategy objects
      /*!
       * The parameters passed to the create() through the \c solverParams
       * argument method should specify the "Transpose Solver Method"
       * as described below, as well as any additional parameters
       * for the particular strategy.
       * <ul>
       * <li> "Transpose Solver Method" - Name of the method. Valid choices are
       *   <ul>
       *   <li> "Tranpose Preconditioner"
       *        (NOX::Epetra::TransposeLinearSystem::TransposePreconditioner)
       *        [Default]
       *   <li> "Explicit Transpose"
       *        (NOX::Epetra::TransposeLinearSystem::ExplicitTranspose)
       *   <li> "Left Preconditioning"
       *        (NOX::Epetra::TransposeLinearSystem::LeftPreconditioning)
       *   </ul>
       * </ul>
       */
      class Factory {

      public:

    //! Constructor
    Factory(const Teuchos::RCP<LOCA::GlobalData>& global_data);

    //! Destructor
    virtual ~Factory();

    //! Create transpose solver strategy
    /*!
     * \param solverParams [in] Solver parameters as described above
     * \param linsys [in] Linear system solver
     */
    Teuchos::RCP<LOCA::Epetra::TransposeLinearSystem::AbstractStrategy>
    create(
        const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
        const Teuchos::RCP<NOX::Epetra::LinearSystem>& linsys);

    //! Return strategy name given by \c solverParams
    const std::string& strategyName(Teuchos::ParameterList& solverParams) const;

      private:

    //! Private to prohibit copying
    Factory(const Factory&);

    //! Private to prohibit copying
    Factory& operator = (const Factory&);

      protected:

    //! Global data object
    Teuchos::RCP<LOCA::GlobalData> globalData;

      }; // Class Factory

    } // NameSpace TransposeLinearSystem

  } // Namespace Epetra

} // Namespace LOCA

#endif // LOCA_EPETRA_TRANSPOSELINEARSYSTEM_FACTORY_H
