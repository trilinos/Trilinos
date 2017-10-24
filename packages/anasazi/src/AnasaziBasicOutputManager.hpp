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

#ifndef ANASAZI_BASIC_OUTPUT_MANAGER_HPP
#define ANASAZI_BASIC_OUTPUT_MANAGER_HPP

/*!     \file AnasaziBasicOutputManager.hpp
        \brief Basic output manager for sending information of select verbosity levels to the appropriate output stream
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziOutputManager.hpp"
#include "Teuchos_oblackholestream.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/*!  \class Anasazi::BasicOutputManager

  \brief Anasazi's basic output manager for sending information of select verbosity levels
  to the appropriate output stream.

  \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist
*/

namespace Anasazi {

  template <class ScalarType>
  class BasicOutputManager : public OutputManager<ScalarType> {

    public:

      //! @name Constructors/Destructor
      //@{ 

      //! Default constructor
      BasicOutputManager(int vb = Anasazi::Errors, int rootRank = 0);

      //! Constructor with specified verbosity and formatted output stream
      BasicOutputManager(int vb,
                         const Teuchos::RCP<Teuchos::FancyOStream>& fos)
      : OutputManager<ScalarType>(vb, fos)
      {};

      //! Destructor.
      virtual ~BasicOutputManager() {};
      //@}


    private:

      //! @name Undefined methods
      //@{ 

      //! Copy constructor.
      BasicOutputManager( const OutputManager<ScalarType>& OM );

      //! Assignment operator.
      BasicOutputManager<ScalarType>& operator=( const OutputManager<ScalarType>& OM );

      //@}
  };

  template<class ScalarType>
  BasicOutputManager<ScalarType>::BasicOutputManager(int vb, int rootRank)
  : OutputManager<ScalarType>(vb)
  {
#ifdef HAVE_MPI
    // The OutputManger constructor will create fos_ that outputs to std::cout
    // This class will query MPI to print on processor 0, if parallel.
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
    this->fos_->setProcRankAndSize(myRank, numProcs);
    this->fos_->setOutputToRootOnly(rootRank);
#endif
  }

} // end Anasazi namespace

#endif

// end of file AnasaziOutputManager.hpp
