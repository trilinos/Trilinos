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

#ifdef HAVE_MPI
#include <Teuchos_CommHelpers.hpp>
#endif

#include <Xpetra_Export.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>

#include "MueLu_PerfUtils_decl.hpp"

//#include "MueLu_Utilities.hpp"

namespace MueLu {

template <class Type>
void calculateStats(Type& minVal, Type& maxVal, double& avgVal, double& devVal, int& minProc, int& maxProc, const RCP<const Teuchos::Comm<int> >& comm, int numActiveProcs, const Type& v) {
  Type sumVal, sum2Val, v2 = v * v;
  using MT    = typename Teuchos::ScalarTraits<Type>::magnitudeType;
  double zero = Teuchos::ScalarTraits<double>::zero();

  MueLu_sumAll(comm, v, sumVal);
  MueLu_sumAll(comm, v2, sum2Val);
  MueLu_minAll(comm, v, minVal);
  MueLu_maxAll(comm, v, maxVal);

  int w;
  w = (minVal == v) ? comm->getRank() : -1;
  MueLu_maxAll(comm, w, maxProc);
  w = (maxVal == v) ? comm->getRank() : -1;
  MueLu_maxAll(comm, w, minProc);

  avgVal       = (numActiveProcs > 0 ? (as<double>(Teuchos::ScalarTraits<Type>::real(sumVal)) / numActiveProcs) : zero);
  MT avgVal_MT = Teuchos::as<MT>(avgVal);
  devVal       = (numActiveProcs > 1 ? sqrt((as<double>(Teuchos::ScalarTraits<Type>::real(sum2Val - sumVal * avgVal_MT))) / (numActiveProcs - 1)) : zero);
}

template <class Type>
std::string stringStats(const RCP<const Teuchos::Comm<int> >& comm, int numActiveProcs, const Type& v, RCP<ParameterList> paramList = Teuchos::null) {
  Type minVal, maxVal;
  double avgVal, devVal;
  int minProc, maxProc;
  calculateStats<Type>(minVal, maxVal, avgVal, devVal, minProc, maxProc, comm, numActiveProcs, v);

  const double zero = Teuchos::ScalarTraits<double>::zero();
  const double one  = Teuchos::ScalarTraits<double>::one();
  std::ostringstream buf;
  buf << std::fixed;
  if ((avgVal != zero) && (paramList.is_null() || !paramList->isParameter("print abs") || paramList->get<bool>("print abs") == false)) {
    double relDev = (devVal / avgVal) * 100;
    double relMin = (as<double>(Teuchos::ScalarTraits<Type>::real(minVal)) / avgVal - one) * 100;
    double relMax = (as<double>(Teuchos::ScalarTraits<Type>::real(maxVal)) / avgVal - one) * 100;
    buf << "avg = " << std::scientific << std::setw(10) << std::setprecision(2) << avgVal << ",  "
        << "dev = " << std::fixed << std::setw(6) << std::setprecision(1) << relDev << "%,  "
        << "min = " << std::fixed << std::setw(7) << std::setprecision(1) << std::setw(7) << relMin << "%"
        << " (" << std::scientific << std::setw(10) << std::setprecision(2) << minVal << " on " << std::fixed << std::setw(4) << minProc << "),  "
        << "max = " << std::fixed << std::setw(7) << std::setprecision(1) << relMax << "%"
        << " (" << std::scientific << std::setw(10) << std::setprecision(2) << maxVal << " on " << std::fixed << std::setw(4) << maxProc << ")";
  } else {
    double relDev = (avgVal != zero ? (devVal / avgVal) * 100 : zero);
    buf << "avg = " << std::scientific << std::setw(10) << std::setprecision(2) << avgVal << ",  "
        << "dev = " << std::fixed << std::setw(6) << std::setprecision(1) << relDev << "%,  "
        << "min = " << std::scientific << std::setw(10) << std::setprecision(2) << minVal
        << "  (on " << std::fixed << std::setw(4) << minProc << "),          "
        << "max = " << std::scientific << std::setw(10) << std::setprecision(2) << maxVal
        << "  (on " << std::fixed << std::setw(4) << maxProc << ")";
  }
  return buf.str();
}

template <class Map>
bool cmp_less(typename Map::value_type& v1, typename Map::value_type& v2) {
  return v1.second < v2.second;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string PerfUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::PrintMatrixInfo(const Matrix& A, const std::string& msgTag, RCP<const ParameterList> params) {
  if (!CheckMatrix(A))
    return "";

  typedef Xpetra::global_size_t global_size_t;

  std::ostringstream ss;

  ss << msgTag << " size =  " << A.getGlobalNumRows() << " x " << A.getGlobalNumCols();
  if (A.haveGlobalConstants())
    ss << ", nnz = " << A.getGlobalNumEntries();
  ss << std::endl;

  if (params.is_null())
    return ss.str();

  bool printLoadBalanceInfo = false, printCommInfo = false, printEntryStats = false;
  if (params->isParameter("printLoadBalancingInfo") && params->get<bool>("printLoadBalancingInfo"))
    printLoadBalanceInfo = true;
  if (params->isParameter("printCommInfo") && params->get<bool>("printCommInfo"))
    printCommInfo = true;
  if (params->isParameter("printEntryStats") && params->get<bool>("printEntryStats"))
    printEntryStats = true;

  if (!printLoadBalanceInfo && !printCommInfo && !printEntryStats)
    return ss.str();

  RCP<const Import> importer = A.getCrsGraph()->getImporter();
  RCP<const Export> exporter = A.getCrsGraph()->getExporter();

  size_t numMyNnz = A.getLocalNumEntries(), numMyRows = A.getLocalNumRows();

  // Create communicator only for active processes
  RCP<const Teuchos::Comm<int> > origComm = A.getRowMap()->getComm();
  bool activeProc                         = true;
  int numProc                             = origComm->getSize();
  int numActiveProcs                      = 0;
#ifdef HAVE_MPI
  RCP<const Teuchos::MpiComm<int> > mpiComm = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(origComm);
  MPI_Comm rawComm                          = (*mpiComm->getRawMpiComm())();

  std::vector<size_t> numRowsPerProc(numProc);
  Teuchos::gatherAll(*origComm, 1, &numMyRows, numProc, &numRowsPerProc[0]);

  int root      = 0;
  bool rootFlag = true;
  for (int i = 0; i < numProc; i++) {
    if (numRowsPerProc[i]) {
      ++numActiveProcs;
      if (rootFlag) {
        root     = i;
        rootFlag = false;
      }
    }
  }

  if (numMyRows == 0) {
    activeProc = false;
    numMyNnz   = 0;
  }  // Reset numMyNnz to avoid adding it up in reduceAll
#else
  if (numMyRows == 0) {
    // FIXME JJH 10-May-2017 Is there any case in serial where numMyRows would be zero?
    //  Reset numMyNnz to avoid adding it up in reduceAll
    numActiveProcs = 0;
    activeProc     = false;
    numMyNnz       = 0;
  } else {
    numActiveProcs = 1;
  }
#endif

  std::string outstr;
  ParameterList absList;
  absList.set("print abs", true);

  RCP<const Matrix> rcpA            = rcpFromRef(A);
  RCP<const CrsMatrixWrap> crsWrapA = rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(rcpA);
  RCP<const CrsMatrix> crsA;
  if (!crsWrapA.is_null())
    crsA = crsWrapA->getCrsMatrix();
  if (printEntryStats && !crsA.is_null()) {
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType magnitudeType;
    typedef Teuchos::ScalarTraits<magnitudeType> MTS;
    ArrayRCP<const size_t> rowptr_RCP;
    ArrayRCP<const LocalOrdinal> colind_RCP;
    ArrayRCP<const Scalar> vals_RCP;
    ArrayRCP<size_t> offsets_RCP;
    ArrayView<const size_t> rowptr;
    ArrayView<const Scalar> vals;
    ArrayView<size_t> offsets;

    crsA->getAllValues(rowptr_RCP, colind_RCP, vals_RCP);
    crsA->getLocalDiagOffsets(offsets_RCP);
    rowptr  = rowptr_RCP();
    vals    = vals_RCP();
    offsets = offsets_RCP();

    Scalar val, minVal, maxVal;
    magnitudeType absVal, minAbsVal, maxAbsVal;
    {
      minVal    = STS::rmax();
      maxVal    = STS::rmin();
      minAbsVal = MTS::rmax();
      maxAbsVal = MTS::zero();

      for (int i = 0; i < offsets.size(); i++) {
        val = vals[rowptr[i] + offsets[i]];
        if (STS::real(val) < STS::real(minVal))
          minVal = val;
        if (STS::real(val) > STS::real(maxVal))
          maxVal = val;
        absVal    = STS::magnitude(val);
        minAbsVal = std::min(minAbsVal, absVal);
        maxAbsVal = std::max(maxAbsVal, absVal);
      }

      ss << msgTag << " diag min       : " << stringStats<Scalar>(origComm, numActiveProcs, minVal) << std::endl;
      ss << msgTag << " diag max       : " << stringStats<Scalar>(origComm, numActiveProcs, maxVal) << std::endl;
      ss << msgTag << " abs(diag) min  : " << stringStats<Scalar>(origComm, numActiveProcs, minAbsVal) << std::endl;
      ss << msgTag << " abs(diag) max  : " << stringStats<Scalar>(origComm, numActiveProcs, maxAbsVal) << std::endl;
    }

    {
      minVal    = STS::rmax();
      maxVal    = STS::rmin();
      minAbsVal = MTS::rmax();
      maxAbsVal = MTS::zero();

      for (int i = 0; i < vals.size(); i++) {
        val = vals[i];
        if (STS::real(val) < STS::real(minVal))
          minVal = val;
        if (STS::real(val) > STS::real(maxVal))
          maxVal = val;
        absVal    = STS::magnitude(val);
        minAbsVal = std::min(minAbsVal, absVal);
        maxAbsVal = std::max(maxAbsVal, absVal);
      }

      ss << msgTag << " entry min      : " << stringStats<Scalar>(origComm, numActiveProcs, minVal) << std::endl;
      ss << msgTag << " entry max      : " << stringStats<Scalar>(origComm, numActiveProcs, maxVal) << std::endl;
      ss << msgTag << " abs(entry) min : " << stringStats<Scalar>(origComm, numActiveProcs, minAbsVal) << std::endl;
      ss << msgTag << " abs(entry) max : " << stringStats<Scalar>(origComm, numActiveProcs, maxAbsVal) << std::endl;
    }
  }

  if (printLoadBalanceInfo) {
    ss << msgTag << " Load balancing info" << std::endl;
    ss << msgTag << "   # active processes: " << numActiveProcs << "/" << numProc << std::endl;
    ss << msgTag << "   # rows per proc   : " << stringStats<global_size_t>(origComm, numActiveProcs, numMyRows) << std::endl;
    ss << msgTag << "   #  nnz per proc   : " << stringStats<global_size_t>(origComm, numActiveProcs, numMyNnz) << std::endl;
  }

  if (printCommInfo && numActiveProcs != 1) {
    typedef std::map<int, size_t> map_type;
    map_type neighMap;
    if (!importer.is_null()) {
      ArrayView<const int> exportPIDs = importer->getExportPIDs();
      if (exportPIDs.size())
        for (int i = 0; i < exportPIDs.size(); i++)
          neighMap[exportPIDs[i]]++;
    }

    // Communication volume
    size_t numExportSend = 0;
    size_t numImportSend = 0;
    size_t numMsgs       = 0;
    size_t minMsg        = 0;
    size_t maxMsg        = 0;

    if (activeProc) {
      numExportSend               = (!exporter.is_null() ? exporter->getNumExportIDs() : 0);
      numImportSend               = (!importer.is_null() ? importer->getNumExportIDs() : 0);
      numMsgs                     = neighMap.size();
      map_type::const_iterator it = std::min_element(neighMap.begin(), neighMap.end(), cmp_less<map_type>);
      minMsg                      = (it != neighMap.end() ? it->second : 0);
      it                          = std::max_element(neighMap.begin(), neighMap.end(), cmp_less<map_type>);
      maxMsg                      = (it != neighMap.end() ? it->second : 0);
    }

    ss << msgTag << " Communication info" << std::endl;
    ss << msgTag << "   # num export send : " << stringStats<global_size_t>(origComm, numActiveProcs, numExportSend) << std::endl;
    ss << msgTag << "   # num import send : " << stringStats<global_size_t>(origComm, numActiveProcs, numImportSend) << std::endl;
    ss << msgTag << "   # num msgs        : " << stringStats<global_size_t>(origComm, numActiveProcs, numMsgs, rcpFromRef(absList)) << std::endl;
    ss << msgTag << "   # min msg size    : " << stringStats<global_size_t>(origComm, numActiveProcs, minMsg) << std::endl;
    ss << msgTag << "   # max msg size    : " << stringStats<global_size_t>(origComm, numActiveProcs, maxMsg) << std::endl;
  }

  outstr = ss.str();

#ifdef HAVE_MPI
  int strLength = outstr.size();
  MPI_Bcast(&strLength, 1, MPI_INT, root, rawComm);
  if (origComm->getRank() != root)
    outstr.resize(strLength);
  MPI_Bcast(&outstr[0], strLength, MPI_CHAR, root, rawComm);
#endif

  return outstr;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string PerfUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::PrintImporterInfo(RCP<const Import> importer, const std::string& msgTag) {
  typedef Xpetra::global_size_t global_size_t;

  std::ostringstream ss;

  // Create communicator only for active processes
  RCP<const Teuchos::Comm<int> > origComm = importer->getSourceMap()->getComm();
  bool activeProc                         = true;
  int numActiveProcs                      = origComm->getSize();
#ifdef HAVE_MPI
  RCP<const Teuchos::MpiComm<int> > mpiComm = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(origComm);
  MPI_Comm rawComm                          = (*mpiComm->getRawMpiComm())();
  int root                                  = 0;
#endif

  std::string outstr;
  ParameterList absList;
  absList.set("print abs", true);

  typedef std::map<int, size_t> map_type;
  map_type neighMap;
  ArrayView<const int> exportPIDs = importer->getExportPIDs();
  if (exportPIDs.size())
    for (int i = 0; i < exportPIDs.size(); i++)
      neighMap[exportPIDs[i]]++;

  // Communication volume
  size_t numImportSend = 0;
  size_t numMsgs       = 0;
  size_t minMsg        = 0;
  size_t maxMsg        = 0;

  if (activeProc) {
    numImportSend               = importer->getNumExportIDs();
    numMsgs                     = neighMap.size();
    map_type::const_iterator it = std::min_element(neighMap.begin(), neighMap.end(), cmp_less<map_type>);
    minMsg                      = (it != neighMap.end() ? it->second : 0);
    it                          = std::max_element(neighMap.begin(), neighMap.end(), cmp_less<map_type>);
    maxMsg                      = (it != neighMap.end() ? it->second : 0);
  }

  ss << msgTag << " Communication info" << std::endl;
  ss << msgTag << "   # num import send : " << stringStats<global_size_t>(origComm, numActiveProcs, numImportSend) << std::endl;
  ss << msgTag << "   # num msgs        : " << stringStats<global_size_t>(origComm, numActiveProcs, numMsgs, rcpFromRef(absList)) << std::endl;
  ss << msgTag << "   # min msg size    : " << stringStats<global_size_t>(origComm, numActiveProcs, minMsg) << std::endl;
  ss << msgTag << "   # max msg size    : " << stringStats<global_size_t>(origComm, numActiveProcs, maxMsg) << std::endl;

  outstr = ss.str();

#ifdef HAVE_MPI
  int strLength = outstr.size();
  MPI_Bcast(&strLength, 1, MPI_INT, root, rawComm);
  if (origComm->getRank() != root)
    outstr.resize(strLength);
  MPI_Bcast(&outstr[0], strLength, MPI_CHAR, root, rawComm);
#endif

  return outstr;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string PerfUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CommPattern(const Matrix& A, const std::string& msgTag, RCP<const ParameterList> /* params */) {
  if (!CheckMatrix(A))
    return "";

  std::ostringstream out;

  RCP<const Teuchos::Comm<int> > comm = A.getRowMap()->getComm();
  int myRank                          = comm->getRank();

  out << msgTag << " " << myRank << ":";

  RCP<const Import> importer = (A.getCrsGraph() != Teuchos::null ? A.getCrsGraph()->getImporter() : Teuchos::null);
  if (importer.is_null()) {
    out << std::endl;
    return out.str();
  }

  ArrayView<const int> exportPIDs = importer->getExportPIDs();

  if (exportPIDs.size()) {
    // NodeTE: exportPIDs is sorted but not unique ( 1 1 1 2 2 3 4 4 4 )
    int neigh = exportPIDs[0];
    GO weight = 1;
    for (int i = 1; i < exportPIDs.size(); i++) {
      if (exportPIDs[i] != exportPIDs[i - 1]) {
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

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool PerfUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CheckMatrix(const Matrix& A) {
  // We can only print statistics for matrices that have a crs graph. A
  // potential issue is regarding Xpetra::TpetraBlockCrsMatrix which has no
  // CrsGraph.  It is held as a private data member by Xpetra::CrsMatrix,
  // which itself is an Xpetra::Matrix. So we check directly whether the
  // request for the graph throws.
  bool hasCrsGraph = true;
  try {
    A.getCrsGraph();

  } catch (...) {
    hasCrsGraph = false;
  }

  return hasCrsGraph;
}

}  // namespace MueLu

#endif  // MUELU_PERFUTILS_DEF_HPP
