// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ANASAZI_OUTPUT_STREAM_TRAITS_HPP
#define ANASAZI_OUTPUT_STREAM_TRAITS_HPP

/*!     \file AnasaziOutputStreamTraits.hpp
        \brief Abstract class definition for Anasazi output stream.
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RCP.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#include "AnasaziGlobalComm.hpp"
#endif

/*!  \class Anasazi::OutputStreamTraits

  \brief Output managers remove the need for the eigensolver to know any information
  about the required output.  However, a formatted output stream is needed to control
  the output during parallel computations.

  \author Mark Hoemmen and Heidi Thornquist
*/

namespace Anasazi {

template<class OperatorType>
struct OutputStreamTraits {

  // The input argument, op, is presumed to be a valid object that can be queried to
  // determine the correct output stream.
  static Teuchos::RCP<Teuchos::FancyOStream>
  getOutputStream (const OperatorType& /* op */, int rootRank = 0)
  {
    Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));

#ifdef HAVE_MPI
    // The default implementation will output on processor 0, if parallel.
    int myRank = 0;
    int numProcs = 1;
    // Initialize MPI
    int mpiStarted = 0;
    MPI_Initialized(&mpiStarted);
    if (mpiStarted)
    {
      MPI_Comm_rank(get_global_comm(), &myRank);
      MPI_Comm_size(get_global_comm(), &numProcs);
    }
    fos->setProcRankAndSize(myRank, numProcs);
    fos->setOutputToRootOnly(rootRank);
#endif

    return fos;
  }
};


} // end Anasazi namespace

#endif

// end of file AnasaziOutputStreamTraits.hpp
