// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RBGen_ISVD_MultiCDUDV.h"

namespace RBGen {

  ISVD_MultiCDUDV::ISVD_MultiCDUDV() : IncSVDPOD(), ISVDUDV(), ISVDMultiCD() {}

  void ISVD_MultiCDUDV::Initialize( 
      const Teuchos::RCP< Teuchos::ParameterList >& params,
      const Teuchos::RCP< const Epetra_MultiVector >& init,
      const Teuchos::RCP< RBGen::FileIOHandler< Epetra_Operator > >& fileio) {
    IncSVDPOD::Initialize(params,init,fileio);
    ISVDUDV::Initialize(params,init,fileio);
    ISVDMultiCD::Initialize(params,init,fileio);
  }

} // end of RBGen namespace
