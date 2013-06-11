// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// ************************************************************************
// @HEADER

#include <Tpetra_Util.hpp>

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

} // namespace Details
} // namespace Tpetra


