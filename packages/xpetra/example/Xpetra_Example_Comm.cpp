// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_Comm.hpp>

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>

#include <Xpetra_EpetraComm.hpp>

#include "Teuchos_DefaultSerialComm.hpp"
#include "Tpetra_Core.hpp"

// This driver simply tests Teuchos2Epetra_Comm

int main(int argc, char** argv) {
  typedef int Ordinal;
  typedef double Scalar;

  using namespace Teuchos;

  oblackholestream blackhole;
  GlobalMPISession mpiSession(&argc, &argv, &blackhole);

  {
    RCP<const Comm<int> > serialTeuchosComm = rcp(new SerialComm<int>);
    RCP<const Comm<int> > teuchosComm       = rcp_implicit_cast<const SerialComm<int> >(serialTeuchosComm);
    RCP<const Epetra_Comm> epetraComm       = Teuchos2Epetra_Comm(teuchosComm);

    assert(epetraComm != Teuchos::null);
  }

  {
    RCP<const Comm<int> > teuchosComm = Tpetra::getDefaultComm();
    RCP<const Epetra_Comm> epetraComm = Teuchos2Epetra_Comm(teuchosComm);

    assert(epetraComm != Teuchos::null);
  }

  return (0);
}  // main
