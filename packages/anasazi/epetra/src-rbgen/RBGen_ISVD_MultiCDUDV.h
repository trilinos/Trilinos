// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RBGEN_ISVD_MULTICDUDV_H
#define RBGEN_ISVD_MULTICDUDV_H

#include "RBGen_ISVDMultiCD.h"
#include "RBGen_ISVDUDV.h"

namespace RBGen {

  //! IncSVD method implementing UDV storage with multiple coordinate descent passes.
  class ISVD_MultiCDUDV : public virtual ISVDUDV, public virtual ISVDMultiCD {
    public:
      //! @name Constructor/Destructor.
      //@{

      //! Default constructor.
      ISVD_MultiCDUDV();

      //! Destructor.
      virtual ~ISVD_MultiCDUDV() {};
      //@}

      //! @name Set Methods
      //@{

      //! Initialize the method with the given parameter list and snapshot set.
      void Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params,
          const Teuchos::RCP< const Epetra_MultiVector >& init,
          const Teuchos::RCP< RBGen::FileIOHandler< Epetra_Operator > >& fileio = Teuchos::null );

      //@}
  };

} // end of RBGen namespace

#endif // RBGEN_ISVD_MULTICDUDV_H
