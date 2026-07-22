// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_PAUSE_TO_ATTACH
#define PANZER_PAUSE_TO_ATTACH

#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_FancyOStream.hpp"

#include <iostream>

#include <sys/types.h>
#include <unistd.h>

namespace panzer {

  inline void pauseToAttach(MPI_Comm mpicomm)
  {
    Teuchos::RCP<Teuchos::Comm<int> > comm =
      Teuchos::createMpiComm<int>(Teuchos::rcp(new Teuchos::OpaqueWrapper<MPI_Comm>(mpicomm)));
    Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
    out.setShowProcRank(true);
    out.setOutputToRootOnly(-1);

    comm->barrier();

    // try to get them to print out all at once
    out << "PID = " << getpid() << std::endl;

    if (comm->getRank() == 0)
      getchar();
    comm->barrier();
  }

  inline void pauseToAttach()
  {
    pauseToAttach(MPI_COMM_WORLD); // CHECK: ALLOW MPI_COMM_WORLD
  }


}

#endif
