/**
  \file   Amesos2_Status.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Thu Jan 14 08:52:04 2010
  
  \brief  Container class for status variables.
*/

#ifndef AMESOS2_STATUS_HPP
#define AMESOS2_STATUS_HPP

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace Amesos {


class Status {
public:
  /// Default constructor.
  Status()
    : isSymbolicFactorizationOK_(false)
    , isNumericFactorizationOK_(false)
    , printTiming_(false)
    , printStatus_(false)
    , computeVectorNorms_(false)
    , computeTrueResidual_(false)
      // Will use the Teuchos::VerboseObject for verbosity settings
    , verbose_(0)
    , debug_(0)

    , numSymbolicFact_(0)
    , numNumericFact_(0)
    , numSolve_(0)
    , threshold_(0.0)

    , myPID_(0)
    , root_(true)
    , numProcs_(1)
    , maxProcs_(1)
    { }


  Status(const Teuchos::RCP<const Teuchos::Comm<int> > comm)
    : isSymbolicFactorizationOK_(false)
    , isNumericFactorizationOK_(false)
    , printTiming_(false)
    , printStatus_(false)
    , computeVectorNorms_(false)
    , computeTrueResidual_(false)
      // Will use the Teuchos::VerboseObject for verbosity settings
    , verbose_(0)
    , debug_(0)

    , numSymbolicFact_(0)
    , numNumericFact_(0)
    , numSolve_(0)
    , threshold_(0.0)

    , myPID_(comm->getRank())
    , root_( myPID_ == 0 )
    , numProcs_(comm->getSize())
    , maxProcs_(1)
    { }


  /// Default destructor.
  ~Status() { };


  /**
   * \brief Set various status parameters.
   *
   * Updates internal variables from the given parameterList.
   */
  void setStatusParameters(
    const Teuchos::RCP<Teuchos::ParameterList> & parameterList );

  /// If \c true, symbolicFactorization() has been successfully called.
  bool isSymbolicFactorizationOK_;

  /// If \c true, numericFactorization() has been successfully called.
  bool isNumericFactorizationOK_;

  /// If \c true, prints timing information in the destructor.
  bool printTiming_;

  /// If \c true, print additional information in the destructor.
  bool printStatus_;

  /// If \c true, prints the norms of X and B in solve().
  bool computeVectorNorms_;

  /// If \c true, computes the true residual in solve().
  bool computeTrueResidual_;

  
  /// Toggles the output level.
  int verbose_;

  /// Sets the level of debug_ output
  int debug_;


  /// Number of symbolic factorization phases.
  int numSymbolicFact_;

  /// Number of numeric factorization phases.
  int numNumericFact_;

  /// Number of solves.
  int numSolve_;  


  double threshold_;


  /// My process ID in this MPI communicator
  int myPID_;

  /// Indicates whether this process is the root process
  bool root_;

  /// The number of processors in this MPI communicator
  int numProcs_;

  /// The maximum number of processes allowed
  int maxProcs_;
};                              // end class Amesos::Status


} // end namespace Amesos

#endif  // AMESOS2_STATUS_HPP
