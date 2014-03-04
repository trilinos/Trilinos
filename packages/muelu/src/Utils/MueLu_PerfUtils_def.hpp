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

#include <algorithm>
#include <string>

#include <Xpetra_Export.hpp>
#include <Xpetra_Import.hpp>

#include "MueLu_PerfUtils_decl.hpp"

#include "MueLu_Utilities.hpp"

namespace MueLu {

  template<class Type>
  void calculateStats(Type& minVal, Type& maxVal, double& avgVal, double& devVal, const RCP<const Teuchos::Comm<int> >& comm, const Type& v, RCP<ParameterList> paramList = Teuchos::null) {
    int numProcs = comm->getSize();
    if (!paramList.is_null() && paramList->isParameter("num procs"))
      numProcs = paramList->get<int>("num procs");

    Type sumVal, sum2Val;

    sumAll(comm,   v, sumVal);
    sumAll(comm, v*v, sum2Val);
    maxAll(comm,   v, maxVal);

    if (paramList.is_null() || !paramList->isParameter("avoid min zero"))
      minAll(comm, v, minVal);
    else if (paramList->get<bool>("avoid min zero") == true)
      minAll(comm, (v > 0 ? v : maxVal), minVal);

    avgVal = as<double>(sumVal) / numProcs;
    devVal = (numProcs != 1 ? sqrt((sum2Val - sumVal*avgVal)/(numProcs-1)) : 0);
  }

  template<class Type>
  std::string stringStats(const RCP<const Teuchos::Comm<int> >& comm, const Type& v, RCP<ParameterList> paramList = Teuchos::null) {
    Type minVal, maxVal;
    double avgVal, devVal;
    calculateStats<Type>(minVal, maxVal, avgVal, devVal, comm, v, paramList);

    char buf[256];
    if (avgVal && (paramList.is_null() || !paramList->isParameter("print abs") || paramList->get<bool>("print abs") == false))
      sprintf(buf, "avg = %.2e,  dev = %5.1f%%,  min = %+6.1f%%,  max = %+6.1f%%", avgVal,
              (devVal/avgVal)*100, (minVal/avgVal-1)*100, (maxVal/avgVal-1)*100);
    else
      sprintf(buf, "avg = %8.2f,  dev = %6.2f,  min = %6.1f ,  max = %6.1f", avgVal,
              devVal, as<double>(minVal), as<double>(maxVal));
    return buf;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string PerfUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::PrintMatrixInfo(const Matrix& A, const std::string& msgTag, RCP<const ParameterList> params) {
    typedef Xpetra::global_size_t global_size_t;

    std::ostringstream ss;
    ss << msgTag << " size =  " << A.getGlobalNumRows() << " x " << A.getGlobalNumCols() << ", nnz = " << A.getGlobalNumEntries() << std::endl;

    if (params.is_null())
      return ss.str();

    RCP<const Teuchos::Comm<int> > comm = A.getRowMap()->getComm();

    if (params->isParameter("printLoadBalancingInfo") && params->get<bool>("printLoadBalancingInfo")) {
      GO numProcessesWithData = 0;

      size_t numMyNnz = A.getNodeNumEntries(), numMyRows = A.getNodeNumRows();

      sumAll(comm, as<GO>((numMyRows > 0) ? 1 : 0), numProcessesWithData);

      ParameterList paramList;
      paramList.set("num procs",      numProcessesWithData);
      paramList.set("avoid min zero", true);

      ss << msgTag << " Load balancing info"    << std::endl;
      ss << msgTag << "   # active processes: " << comm->getSize() << ",  # processes with data = " << numProcessesWithData << std::endl;
      ss << msgTag << "   # rows per proc   : " << stringStats<global_size_t>(comm, numMyRows, rcpFromRef(paramList)) << std::endl;
      ss << msgTag << "   #  nnz per proc   : " << stringStats<global_size_t>(comm,  numMyNnz, rcpFromRef(paramList)) << std::endl;
    }

    if (params->isParameter("printCommInfo") && params->get<bool>("printCommInfo")) {
      RCP<const Import> importer = A.getCrsGraph()->getImporter();
      RCP<const Export> exporter = A.getCrsGraph()->getExporter();

      // Communication volume
      size_t numExport = 0, numImport = 0;
      size_t numNeigh  = 0, maxNeighMsg = 0;
      if (!importer.is_null()) {
        numExport = importer->getNumExportIDs();
        numImport = importer->getNumRemoteIDs();

        ArrayView<const int> exportPIDs = importer->getExportPIDs();

        if (exportPIDs.size()) {
          std::map<int,size_t> neighMap;
          for (int i = 0; i < exportPIDs.size(); i++)
            neighMap[exportPIDs[i]]++;

          numNeigh    = neighMap.size();
          maxNeighMsg = std::max_element(neighMap.begin(), neighMap.end())->second;
        }
      }

      ParameterList absList;
      absList.set("print abs", true);

      ss << msgTag << " Communication info"     << std::endl;
      ss << msgTag << "   # num export ids  : " << stringStats<global_size_t>(comm,   numExport)                      << std::endl;
      ss << msgTag << "   # num import ids  : " << stringStats<global_size_t>(comm,   numImport)                      << std::endl;
      ss << msgTag << "   # num neighbors   : " << stringStats<global_size_t>(comm,    numNeigh, rcpFromRef(absList)) << std::endl;
      ss << msgTag << "   # max neigh msg   : " << stringStats<global_size_t>(comm, maxNeighMsg)                      << std::endl;
    }

    return ss.str();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string PerfUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CommPattern(const Matrix& A, const std::string& msgTag, RCP<const ParameterList> params) {
    std::ostringstream out;

    RCP<const Teuchos::Comm<int> > comm = A.getRowMap()->getComm();
    int myRank = comm->getRank();

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
      GO  weight = 1;
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
