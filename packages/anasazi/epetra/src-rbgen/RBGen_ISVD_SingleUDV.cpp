// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RBGen_ISVD_SingleUDV.h"

namespace RBGen {

  ISVD_SingleUDV::ISVD_SingleUDV() : IncSVDPOD(), ISVDUDV(), ISVDSingle() {}

  void ISVD_SingleUDV::Initialize(
      const Teuchos::RCP< Teuchos::ParameterList >& params,
      const Teuchos::RCP< const Epetra_MultiVector >& init,
      const Teuchos::RCP< RBGen::FileIOHandler< Epetra_Operator > >& fileio) {
    IncSVDPOD::Initialize(params,init,fileio);
    ISVDUDV::Initialize(params,init,fileio);
    ISVDSingle::Initialize(params,init,fileio);
  }

} // end of RBGen namespace
