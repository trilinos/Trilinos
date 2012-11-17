// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include "MueLu_TimeMonitor.hpp"

namespace MueLu {

  // TODO: this function can be templated (T=double).
  ArrayRCP<double> ReduceMaxMinAvg(double localValue, Teuchos::Comm<int> const &comm, int rootNode) {
    ArrayRCP<double> r = ArrayRCP<double>(3, localValue);

#ifdef HAVE_MPI
    double & maxTime = r[0], & minTime = r[1], & avgTime = r[2];

    // Note: workaround because reduce() is not implemented in Teuchos::Comm
    const Teuchos::MpiComm<int> & mpiComm = dynamic_cast<const Teuchos::MpiComm<int>& >(comm);
    MPI_Comm rawMpiComm = (*mpiComm.getRawMpiComm())();
    //

    // DEBUG std::cout << comm.getRank() << ": " << localValue << std::endl;

    int ntimers=1;
    MPI_Reduce(&localValue, &maxTime, ntimers, MPI_DOUBLE, MPI_MAX, rootNode, rawMpiComm);
    MPI_Reduce(&localValue, &minTime, ntimers, MPI_DOUBLE, MPI_MIN, rootNode, rawMpiComm);
    MPI_Reduce(&localValue, &avgTime, ntimers, MPI_DOUBLE, MPI_SUM, rootNode, rawMpiComm); avgTime /= comm.getSize();
#endif // HAVE_MPI

    return r;
  }

} // namespace MueLu
