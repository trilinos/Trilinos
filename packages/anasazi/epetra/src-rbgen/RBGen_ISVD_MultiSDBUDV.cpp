// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RBGen_ISVD_MultiSDBUDV.h"

namespace RBGen {

  ISVD_MultiSDBUDV::ISVD_MultiSDBUDV() : IncSVDPOD(), ISVDUDV(), ISVDMultiSDB() {}

  void ISVD_MultiSDBUDV::Initialize( 
      const Teuchos::RCP< Teuchos::ParameterList >& params,
      const Teuchos::RCP< const Epetra_MultiVector >& init,
      const Teuchos::RCP< RBGen::FileIOHandler< Epetra_Operator > >& fileio) {
    IncSVDPOD::Initialize(params,init,fileio);
    ISVDUDV::Initialize(params,init,fileio);
    ISVDMultiSDB::Initialize(params,init,fileio);
  }

} // end of RBGen namespace
