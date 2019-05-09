// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
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
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
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
