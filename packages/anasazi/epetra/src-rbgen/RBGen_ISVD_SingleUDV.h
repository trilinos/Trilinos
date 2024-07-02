// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RBGEN_ISVD_SINGLEUDV_H
#define RBGEN_ISVD_SINGLEUDV_H

#include "RBGen_ISVDSingle.h"
#include "RBGen_ISVDUDV.h"

namespace RBGen {

  //! IncSVD method implementing UDV storage with a single coordinate descent pass.
  class ISVD_SingleUDV : public virtual ISVDUDV, public virtual ISVDSingle {
    public:

      //! @name Constructor/Destructor.
      //@{

      //! Default constructor.
      ISVD_SingleUDV();

      //! Destructor.
      virtual ~ISVD_SingleUDV() {};
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

#endif // RBGEN_ISVD_SINGLEUDV_H
