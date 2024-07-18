// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Util.hpp"
#include "Teuchos_Comm.hpp"

namespace Tpetra {
namespace Details {

/// \brief Whether the two communicators are congruent.
///
/// Two communicators are <i>congruent</i> when they have the same
/// number of processes, and those processes occur in the same rank
/// order.
///
/// If both communicators are MpiComm instances, this function returns
/// <tt>true</tt> exactly when <tt>MPI_Comm_compare</tt> returns
/// <tt>MPI_IDENT</tt> (the communicators are handles for the same
/// object) or <tt>MPI_CONGRUENT</tt>.  SerialComm instances are
/// always congruent.  An MpiComm is congruent to a SerialComm if the
/// MpiComm has only one process.  This function is symmetric in its
/// arguments.
///
/// If either Comm instance is neither an MpiComm nor a SerialComm,
/// this method cannot do any better than to compare their process
/// counts.
bool
congruent (const Teuchos::Comm<int>& comm1,
           const Teuchos::Comm<int>& comm2)
{
#ifdef HAVE_MPI
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::rcpFromRef;
  using Teuchos::MpiComm;
  using Teuchos::rcp_dynamic_cast;

  RCP<const MpiComm<int> > mpiComm1 =
    rcp_dynamic_cast<const MpiComm<int> > (rcpFromRef (comm1));
  RCP<const MpiComm<int> > mpiComm2 =
    rcp_dynamic_cast<const MpiComm<int> > (rcpFromRef (comm2));

  if (mpiComm1.is_null ()) { // comm1 is not an MpiComm
    return comm1.getSize () == comm2.getSize (); // hope for the best
  } else { // comm1 is an MpiComm
    if (mpiComm2.is_null ()) { // comm2 is not an MpiComm
      return comm1.getSize () == comm2.getSize (); // hope for the best
    } else { // both comm1 and comm2 are MpiComm
      MPI_Comm rawMpiComm1 = * (mpiComm1->getRawMpiComm ());
      MPI_Comm rawMpiComm2 = * (mpiComm2->getRawMpiComm ());

      int result = MPI_UNEQUAL;
      const int err = MPI_Comm_compare (rawMpiComm1, rawMpiComm2, &result);
      TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
                                 "congruent: MPI_Comm_compare failed");
      return result == MPI_IDENT || result == MPI_CONGRUENT;
    }
  }
#else // NOT HAVE_MPI
  return comm1.getSize () == comm2.getSize (); // hope for the best
#endif // HAVE_MPI
}

std::unique_ptr<std::string>
createPrefix(const int myRank,
             const char prefix[])
{
  std::ostringstream os;
  os << "Proc " << myRank << ": " << prefix << ": ";
  return std::unique_ptr<std::string>(new std::string(os.str()));
}

std::unique_ptr<std::string>
createPrefix(const Teuchos::Comm<int>* comm,
             const char functionName[])
{
  const int myRank = comm == nullptr ? -1 : comm->getRank();
  const std::string prefix = std::string("Tpetra::") + functionName;
  return createPrefix(myRank, prefix.c_str());
}

std::unique_ptr<std::string>
createPrefix(const Teuchos::Comm<int>* comm,
             const char className[],
             const char methodName[])
{
  const int myRank = comm == nullptr ? -1 : comm->getRank();
  const std::string prefix = std::string("Tpetra::") +
    className + std::string("::") + methodName;
  return createPrefix(myRank, prefix.c_str());
}

} // namespace Details
} // namespace Tpetra
