// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Program to debug segfaults being reported in CDASH when
// -D KokkosClassic_DefaultNode:STRING=Tpetra::KokkosCompat::KokkosOpenMPWrapperNode 
// -D Trilinos_ENABLE_OpenMP:BOOL=ON  
// Problem appears to be in creation of Xpetra::EpetraMapT

#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Comm.hpp>

#include <Tpetra_Map.hpp>
#include <Epetra_Map.h>
#include <Xpetra_EpetraMap.hpp>
#include <Xpetra_EpetraUtils.hpp>

int main(int narg, char **arg)
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > tcomm = Tpetra::getDefaultComm();
  Teuchos::RCP<const Epetra_Comm> ecomm = Xpetra::toEpetra(tcomm);

  const int nGlobRows = 50;
  const Epetra_Map emap(nGlobRows, 0, *ecomm);
  Teuchos::RCP<const Epetra_BlockMap> ebmap = Teuchos::rcpFromRef(emap);

  typedef Xpetra::EpetraMapT<int, Tpetra::Map<>::node_type> xemap_t;
  Teuchos::RCP<const xemap_t> xmap(new xemap_t(ebmap));
  
  const Teuchos::RCP<const Teuchos::Comm<int> > &xcomm = xmap->getComm();

  std::cout << "Teuchos:  Hello from " 
            << tcomm->getRank() << " of " 
            << tcomm->getSize() << std::endl;
  std::cout << "Epetra:   Hello from " 
            << ecomm->MyPID() << " of " 
            << ecomm->NumProc() << std::endl;
  std::cout << "Xpetra:   Hello from " 
            << xcomm->getRank() << " of " 
            << xcomm->getSize() << std::endl;
}
