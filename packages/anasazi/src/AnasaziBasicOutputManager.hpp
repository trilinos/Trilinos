// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "AnasaziGlobalComm.hpp"
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
      MPI_Comm_rank(get_global_comm(), &myRank);
      MPI_Comm_size(get_global_comm(), &numProcs);
    }
    this->fos_->setProcRankAndSize(myRank, numProcs);
    this->fos_->setOutputToRootOnly(rootRank);
#endif
  }

} // end Anasazi namespace

#endif

// end of file AnasaziOutputManager.hpp
