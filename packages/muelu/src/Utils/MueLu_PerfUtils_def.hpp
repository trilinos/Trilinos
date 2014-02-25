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
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NodeT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DIScalarLAIMED. IN Node EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NodeT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GlobalOrdinalODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_PERFUTILS_DEF_HPP
#define MUELU_PERFUTILS_DEF_HPP

#include <string>

#include <Xpetra_Import.hpp>

#include "MueLu_PerfUtils_decl.hpp"

#include "MueLu_Utilities.hpp"

namespace MueLu {

  template <class Scalar, class LO, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string PerfUtils<Scalar, LO, GlobalOrdinal, Node, LocalMatOps>::PrintMatrixInfo(const Matrix& A, const std::string& msgTag, RCP<const ParameterList> params) {
    std::ostringstream ss;
    ss << msgTag << " size =  " << A.getGlobalNumRows() << " x " << A.getGlobalNumCols() << ", nnz = " << A.getGlobalNumEntries() << std::endl;

    if (params.is_null())
      return ss.str();

    if (params->isParameter("printLoadBalancingInfo") && params->get<bool>("printLoadBalancingInfo")) {
      RCP<const Teuchos::Comm<int> > comm = A.getRowMap()->getComm();
      GlobalOrdinal numActiveProcesses = comm->getSize(), numProcessesWithData = 0;

      // aggregate data
      GlobalOrdinal  numMyNnz = A.getNodeNumEntries(),     minNnz,     maxNnz;
      GlobalOrdinal numMyRows = A.getNodeNumRows(),    minNumRows, maxNumRows;
      double  numMyNnz2 =  Teuchos::as<double>(numMyNnz)* numMyNnz,     sumNnz = Teuchos::as<double>(A.getGlobalNumEntries()),     sum2Nnz;
      double numMyRows2 = Teuchos::as<double>(numMyRows)*numMyRows, sumNumRows = Teuchos::as<double>(A.getGlobalNumRows()),    sum2NumRows;
      maxAll(comm,                                       numMyNnz, maxNnz);
      maxAll(comm,                                      numMyRows, maxNumRows);
      sumAll(comm, (GlobalOrdinal)((numMyRows > 0) ?         1 :          0), numProcessesWithData);
      minAll(comm, (GlobalOrdinal)(( numMyNnz > 0) ?  numMyNnz :     maxNnz), minNnz);
      minAll(comm, (GlobalOrdinal)((numMyRows > 0) ? numMyRows : maxNumRows), minNumRows);
      sumAll(comm,                                      numMyNnz2, sum2Nnz);
      sumAll(comm,                                     numMyRows2, sum2NumRows);

      double avgNumRows = sumNumRows / numProcessesWithData;
      double avgNnz     = sumNnz     / numProcessesWithData;
      double devNumRows = (numProcessesWithData != 1 ? sqrt((sum2NumRows - sumNumRows*sumNumRows/numProcessesWithData)/(numProcessesWithData-1)) : 0);
      double devNnz     = (numProcessesWithData != 1 ? sqrt((sum2Nnz     -     sumNnz*    sumNnz/numProcessesWithData)/(numProcessesWithData-1)) : 0);

      char buf[256];
      ss << msgTag << " Load balancing info:" << std::endl;
      ss << msgTag << "   # active processes: "   << numActiveProcesses << ",  # processes with data = " << numProcessesWithData << std::endl;
      sprintf(buf, "avg = %.2e,  dev = %4.1f%%,  min = %+5.1f%%,  max = %+5.1f%%", avgNumRows,
              (devNumRows/avgNumRows)*100, (minNumRows/avgNumRows-1)*100, (maxNumRows/avgNumRows-1)*100);
      ss << msgTag << "   # rows per proc   : " << buf << std::endl;
      sprintf(buf, "avg = %.2e,  dev = %4.1f%%,  min = %+5.1f%%,  max = %+5.1f%%", avgNnz,
              (devNnz/avgNnz)*100, (minNnz/avgNnz-1)*100, (maxNnz/avgNnz-1)*100);
      ss << msgTag << "   #  nnz per proc   : " << buf << std::endl;
    }

    return ss.str();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string PerfUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CommPattern(const Matrix& A, const std::string& msgTag, RCP<const ParameterList> params) {
    std::ostringstream out;

    RCP<const Teuchos::Comm<int> > comm = A.getRowMap()->getComm();
    int myRank   = comm->getRank();
    int numProcs = comm->getSize();

    out << msgTag << " " << myRank << ":";

    RCP<const Import> importer = (A.getCrsGraph() != Teuchos::null ? A.getCrsGraph()->getImporter() : Teuchos::null);
    if (importer.is_null()) {
      out << std::endl;
      return out.str();
    }

    ArrayView<const int> exportPIDs = importer->getExportPIDs();

    if (exportPIDs.size()) {
      // NodeTE: exportPIDs is sorted but not unique ( 1 1 1 2 2 3 4 4 4 )
      int neigh  = exportPIDs[0];
      GlobalOrdinal  weight = 1;
      for (int i = 1; i < exportPIDs.size(); i++) {
        if (exportPIDs[i] != exportPIDs[i-1]) {
          out << " " << neigh << "(" << weight << ")";

          neigh  = exportPIDs[i];
          weight = 1;

        } else {
          weight += 1;
        }
      }
      out << " " << neigh << "(" << weight << ")" << std::endl;
    }

    return out.str();
  }

} //namespace MueLu

#endif // MUELU_PERFUTILS_DEF_HPP
