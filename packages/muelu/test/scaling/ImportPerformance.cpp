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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <unistd.h>

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Xpetra
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_IO.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>

#include <MueLu.hpp>

#include <MueLu_BaseClass.hpp>
#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION
#include <MueLu_ExplicitInstantiation.hpp>
#endif
#include <MueLu_Level.hpp>
#include <MueLu_MutuallyExclusiveTime.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Utilities.hpp>

#ifdef HAVE_MUELU_EPETRA
#include <EpetraExt_MMHelpers.h>
#include <EpetraExt_RowMatrixOut.h>
#endif

#include <TpetraExt_MatrixMatrix.hpp>

#include <MueLu_CreateXpetraPreconditioner.hpp>

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::TimeMonitor;

// =========================================================================
// =========================================================================
// =========================================================================
// Provide a "resize" operation for double*'s.
inline void resize_doubles(int nold, int nnew, double*& d) {
  if (nnew > nold) {
    double* tmp = new double[nnew];
    for (int i = 0; i < nold; i++)
      tmp[i] = d[i];
    delete[] d;
    d = tmp;
  }
}

// =========================================================================
// =========================================================================
// =========================================================================
#if defined(HAVE_MUELU_EPETRA)

extern void MakeColMapAndReindexSort(int& NumRemoteColGIDs, int*& RemoteColindices,
                                     std::vector<int>& RemotePermuteIDs, std::vector<int>& RemoteOwningPIDs);

extern void MakeColMapAndReindexSort(int& NumRemoteColGIDs, long long*& RemoteColindices,
                                     std::vector<int>& RemotePermuteIDs, std::vector<int>& RemoteOwningPIDs);

void build_remote_pids(int MyPID, const std::vector<int>& ColMapOwningPIDs, std::vector<int>& RemotePIDs) {
  // Presume the column map has Aztec ordering
  int N = (int)ColMapOwningPIDs.size();
  int first_idx;
  for (first_idx = 0; first_idx < N; first_idx++)
    if (ColMapOwningPIDs[first_idx] != MyPID)
      break;

  /*   printf("[%d] ColMapOwningPIDs(%d) =",MyPID,(int)ColMapOwningPIDs.size());
   for(int i=0;i<(int)ColMapOwningPIDs.size(); i++)
     printf("%d ",ColMapOwningPIDs[i]);
     printf("\n");*/

  // Make sure there are some non-local unknowns
  if (first_idx == N) {
    printf("[%d] No remotes\n", MyPID);
    return;
  }

  RemotePIDs.resize(ColMapOwningPIDs.size() - first_idx);
  for (int i = first_idx; i < N; i++)
    RemotePIDs[i - first_idx] = ColMapOwningPIDs[i];

  /*   printf("[%d] RemotePIDs(%d) =",MyPID,(int)RemotePIDs.size());
   for(int i=0;i<(int)RemotePIDs.size(); i++)
     printf("%d ",RemotePIDs[i]);
     printf("\n");*/
}

Epetra_Map* convert_lightweightmap_to_map(const EpetraExt::LightweightMap& A, const Epetra_Comm& Comm) {
  Epetra_Map* Aout = 0;
  if (A.GlobalIndicesInt()) {
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    Aout = new Epetra_Map(-1, A.NumMyElements(), A.MyGlobalElements(), 0, Comm);
#endif
  } else if (A.GlobalIndicesLongLong()) {
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    Aout = new Epetra_Map((long long)-1, A.NumMyElements(), A.MyGlobalElements64(), (long long)0, Comm);
#endif
  }

  return Aout;
}

Epetra_CrsMatrix* convert_lightweightcrsmatrix_to_crsmatrix(const EpetraExt::LightweightCrsMatrix& A) {
  auto tm                 = TimeMonitor::getNewTimer("OptimizedTransfer: Convert: MapConstructor");
  const Epetra_Comm& Comm = A.DomainMap_.Comm();

  // Build Maps
  Epetra_Map *RowMap, *ColMap;
  if (A.use_lw)
    RowMap = convert_lightweightmap_to_map(*A.RowMapLW_, Comm);
  else
    throw std::runtime_error("Only works in LW mode");
  ColMap                 = convert_lightweightmap_to_map(A.ColMap_, Comm);
  Epetra_CrsMatrix* Aout = new Epetra_CrsMatrix(Copy, *RowMap, *ColMap, 0);
  int N                  = RowMap->NumMyElements();
  int nnz                = A.colind_.size();

  tm = Teuchos::null;
  // Copy pointers over
  auto tm2                            = TimeMonitor::getNewTimer("OptimizedTransfer: Convert: Data Copy");
  Epetra_IntSerialDenseVector& rowptr = Aout->ExpertExtractIndexOffset();
  Epetra_IntSerialDenseVector& colind = Aout->ExpertExtractIndices();
  double*& vals                       = Aout->ExpertExtractValues();
  rowptr.Resize(N + 1);
  colind.Resize(nnz);
  resize_doubles(0, nnz, vals);

  for (int i = 0; i < N + 1; i++)
    rowptr[i] = A.rowptr_[i];

  for (int i = 0; i < nnz; i++) {
    colind[i] = A.colind_[i];
    vals[i]   = A.vals_[i];
  }
  tm2      = Teuchos::null;
  auto tm3 = TimeMonitor::getNewTimer("OptimizedTransfer: Convert: BuildRemote");

  // Get RemotePIDs
  std::vector<int> RemotePIDs_;
  build_remote_pids(Comm.MyPID(), A.ColMapOwningPIDs_, RemotePIDs_);

  tm3      = Teuchos::null;
  auto tm4 = TimeMonitor::getNewTimer("OptimizedTransfer: Convert: BuildImport");

  // Importer build
  const int* ExportLIDs   = A.ExportLIDs_.size() ? &A.ExportLIDs_[0] : 0;
  const int* ExportPIDs   = A.ExportPIDs_.size() ? &A.ExportPIDs_[0] : 0;
  const int* RemotePIDs   = RemotePIDs_.size() ? &RemotePIDs_[0] : 0;
  Epetra_Import* Importer = new Epetra_Import(*ColMap, A.DomainMap_, RemotePIDs_.size(), RemotePIDs, A.ExportLIDs_.size(), ExportLIDs, ExportPIDs);

  tm4      = Teuchos::null;
  auto tm5 = TimeMonitor::getNewTimer("OptimizedTransfer: Convert: ESFC");

  // ESFC
  Aout->ExpertStaticFillComplete(A.DomainMap_, *RowMap, Importer, 0);

  // Cleanup
  if (A.use_lw) delete RowMap;
  delete ColMap;

  return Aout;
}

#endif

// =========================================================================
// =========================================================================
// =========================================================================
#if defined(HAVE_MUELU_EPETRA)
bool epetra_check_importer_correctness(const Epetra_Import& A, const Epetra_Import& B) {
  int MyPID       = A.SourceMap().Comm().MyPID();
  bool is_correct = true;

  // Same
  if (A.NumSameIDs() != B.NumSameIDs()) {
    printf("[%d] NumSameIDs %d vs. %d\n", MyPID, A.NumSameIDs(), B.NumSameIDs());
    is_correct = false;
  }
  // Permutes
  if (A.NumPermuteIDs() != B.NumPermuteIDs()) {
    printf("[%d] NumPermuteIDs %d vs. %d\n", MyPID, A.NumPermuteIDs(), B.NumPermuteIDs());
    is_correct = false;
  } else {
    int N               = A.NumPermuteIDs();
    bool error_detected = false;
    for (int i = 0; !error_detected && i < N; i++)
      error_detected = (A.PermuteFromLIDs()[i] != B.PermuteFromLIDs()[i]) || (A.PermuteToLIDs()[i] != B.PermuteToLIDs()[i]);

    if (error_detected) {
      printf("[%d] A Permutes = ", MyPID);
      for (int i = 0; i < N; i++)
        printf("%d->%d ", A.PermuteFromLIDs()[i], A.PermuteToLIDs()[i]);
      printf("\n[%d] B Permutes = ", MyPID);
      for (int i = 0; i < N; i++)
        printf("%d->%d ", B.PermuteFromLIDs()[i], B.PermuteToLIDs()[i]);
      printf("\n");
      is_correct = false;
    }
  }

  // Remotes
  if (A.NumRemoteIDs() != B.NumRemoteIDs()) {
    printf("[%d] NumRemoteIDs %d vs. %d\n", MyPID, A.NumRemoteIDs(), B.NumRemoteIDs());
    is_correct = false;
  } else {
    int N               = A.NumRemoteIDs();
    bool error_detected = false;
    for (int i = 0; !error_detected && i < N; i++)
      error_detected = A.RemoteLIDs()[i] != B.RemoteLIDs()[i];

    if (error_detected) {
      printf("[%d] A RemoteLIDs = ", MyPID);
      for (int i = 0; i < N; i++)
        printf("%d ", A.RemoteLIDs()[i]);
      printf("\n[%d] B RemoteLIDs = ", MyPID);
      for (int i = 0; i < N; i++)
        printf("%d ", B.RemoteLIDs()[i]);
      printf("\n");
      is_correct = false;
    }
  }

  // Exports
  if (A.NumExportIDs() != B.NumExportIDs()) {
    printf("[%d] NumExportIDs %d vs. %d\n", MyPID, A.NumExportIDs(), B.NumExportIDs());
    is_correct = false;
  } else {
    int N               = A.NumExportIDs();
    bool error_detected = false;
    for (int i = 0; !error_detected && i < N; i++)
      error_detected = (A.ExportLIDs()[i] != B.ExportLIDs()[i]) || (A.ExportPIDs()[i] != B.ExportPIDs()[i]);

    if (error_detected) {
      printf("[%d] A Exports(%d) = ", MyPID, A.NumExportIDs());
      for (int i = 0; i < N; i++)
        printf("%d(%d)->%d ", A.ExportLIDs()[i], A.SourceMap().GID(A.ExportLIDs()[i]), A.ExportPIDs()[i]);
      printf("\n[%d] B Exports(%d) = ", MyPID, B.NumExportIDs());
      for (int i = 0; i < N; i++)
        printf("%d(%d)->%d ", B.ExportLIDs()[i], B.SourceMap().GID(A.ExportLIDs()[i]), B.ExportPIDs()[i]);
      printf("\n");
      is_correct = false;
    }
  }

  // Message Counts
  if (A.NumSend() != B.NumSend()) {
    printf("[%d] NumSend %d vs. %d\n", MyPID, A.NumSend(), B.NumSend());
    is_correct = false;
  }
  if (A.NumRecv() != B.NumRecv()) {
    printf("[%d] NumRecv %d vs. %d\n", MyPID, A.NumRecv(), B.NumRecv());
    is_correct = false;
  }

#ifdef HAVE_MPI
  const Epetra_MpiDistributor& Ad = *dynamic_cast<Epetra_MpiDistributor*>(&A.Distributor());
  const Epetra_MpiDistributor& Bd = *dynamic_cast<Epetra_MpiDistributor*>(&B.Distributor());

  if (Ad.MaxSendLength() != Bd.MaxSendLength()) {
    printf("[%d] Distor.MaxSendLength %d vs. %d\n", MyPID, Ad.MaxSendLength(), Bd.MaxSendLength());
    is_correct = false;
  }
  if (Ad.TotalReceiveLength() != Bd.TotalReceiveLength()) {
    printf("[%d] Distor.TotalReceiveLength %d vs. %d\n", MyPID, Ad.TotalReceiveLength(), Bd.TotalReceiveLength());
    is_correct = false;
  }

  if (Ad.NumSends() != Bd.NumSends()) {
    printf("[%d] Distor.NumSends %d vs. %d\n", MyPID, Ad.NumSends(), Bd.NumSends());
    is_correct = false;
  } else {
    int N               = Ad.NumSends();
    bool error_detected = false;
    for (int i = 0; !error_detected && i < N; i++)
      error_detected = (Ad.ProcsTo()[i] != Bd.ProcsTo()[i]) || (Ad.LengthsTo()[i] != Bd.LengthsTo()[i]);

    if (error_detected) {
      printf("[%d] Ad Sends = ", MyPID);
      for (int i = 0; i < N; i++)
        printf("%d->%d ", Ad.LengthsTo()[i], Ad.ProcsTo()[i]);
      printf("\n[%d] Bd Sends = ", MyPID);
      for (int i = 0; i < N; i++)
        printf("%d->%d ", Bd.LengthsTo()[i], Bd.ProcsTo()[i]);
      printf("\n");
      is_correct = false;
    }
  }

  if (Ad.NumReceives() != Bd.NumReceives()) {
    printf("[%d] Distor.NumReceives %d vs. %d\n", MyPID, Ad.NumReceives(), Bd.NumReceives());
    is_correct = false;
  } else {
    int N               = Ad.NumReceives();
    bool error_detected = false;
    for (int i = 0; !error_detected && i < N; i++)
      error_detected = (Ad.ProcsFrom()[i] != Bd.ProcsFrom()[i]) || (Ad.LengthsFrom()[i] != Bd.LengthsFrom()[i]);

    if (error_detected) {
      printf("[%d] Ad Receives = ", MyPID);
      for (int i = 0; i < N; i++)
        printf("%d->%d ", Ad.LengthsFrom()[i], Ad.ProcsFrom()[i]);
      printf("\n[%d] Bd Receives = ", MyPID);
      for (int i = 0; i < N; i++)
        printf("%d->%d ", Bd.LengthsFrom()[i], Bd.ProcsFrom()[i]);
      printf("\n");
      is_correct = false;
    }
  }
#endif

  return is_correct;
}
#endif  // if defined(HAVE_MUELU_EPETRA)

// =========================================================================
// =========================================================================
// =========================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TestTransfer(Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A, Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > P) {
#include <MueLu_UseShortNames.hpp>
  Xpetra::UnderlyingLib lib = A->getRowMap()->lib();

  if (lib == Xpetra::UseTpetra) {
    typedef Tpetra::CrsMatrixStruct<SC, LO, GO, NO> crs_matrix_struct_type;
    typedef Tpetra::CrsMatrix<SC, LO, GO, NO> crs_matrix_type;
    typedef Tpetra::Import<LO, GO, NO> import_type;

    RCP<const crs_matrix_type> Au = Utilities::Op2TpetraCrs(A);
    RCP<const crs_matrix_type> Pu = Utilities::Op2TpetraCrs(P);
    if (Au->getComm()->getSize() == 1) return;

    // ==================
    // Optimized Transfer
    // ==================
    auto tm = TimeMonitor::getNewTimer("OptimizedTransfer: Import");
    crs_matrix_struct_type Pview;
    Tpetra::MMdetails::import_and_extract_views(*Pu, Au->getColMap(), Pview, Au->getGraph()->getImporter(), false, "ImportPerf: ");

    Au->getComm()->barrier();

    // ==================
    // Naive Transfer
    // ==================
    // Use the columnmap from Aopt and build an importer ex nihilo

    tm       = Teuchos::null;
    auto tm2 = TimeMonitor::getNewTimer("NaiveTransfer: BuildImport");
    import_type NaiveImport(Pview.importMatrix->getColMap(), Pu->getDomainMap());
    Au->getComm()->barrier();
  } else if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_MUELU_EPETRA)
    RCP<const Epetra_CrsMatrix> Au = Utilities::Op2EpetraCrs(A);
    RCP<const Epetra_CrsMatrix> Pu = Utilities::Op2EpetraCrs(P);
    if (Au->Comm().NumProc() == 1) return;

    // ==================
    // Optimized Transfer
    // ==================
    // Build the LightweightCrsMatrix

    auto tm3 = TimeMonitor::getNewTimer("OptimizedTransfer: Import");
    EpetraExt::CrsMatrixStruct Pview;
    bool SortGhosts = true;

    if (Au->RowMap().GlobalIndicesInt()) {
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
      EpetraExt::import_only<int>(*Pu, Au->ColMap(), Pview, Au->Importer(), SortGhosts, "ImportPerf: ");
#endif
    } else if (Au->RowMap().GlobalIndicesInt()) {
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
      EpetraExt::import_only<long long>(*Pu, Au->ColMap(), Pview, Au->Importer(), SortGhosts, "ImportPerf: ");
#endif
    }
    tm3                    = Teuchos::null;
    auto tm4               = TimeMonitor::getNewTimer("OptimizedTransfer: Convert");
    Epetra_CrsMatrix* Aopt = convert_lightweightcrsmatrix_to_crsmatrix(*Pview.importMatrix);

    Au->Comm().Barrier();
    // ==================
    // Naive Transfer
    // ==================
    // Use the columnmap from Aopt and build an importer ex nihilo
    tm4                           = Teuchos::null;
    auto tm5                      = TimeMonitor::getNewTimer("NaiveTransfer: BuildImport");
    const Epetra_Map& NaiveColMap = Aopt->ColMap();
    Epetra_Import NaiveImport(NaiveColMap, Pu->DomainMap());

    Au->Comm().Barrier();

    // Check importer for correctness
    fflush(stdout);
    const Epetra_Import* OptImport = Aopt->Importer();
    bool is_correct                = epetra_check_importer_correctness(NaiveImport, *OptImport);
    fflush(stdout);
    int is_OK_local = is_correct, is_OK_global;
    Au->Comm().MinAll(&is_OK_local, &is_OK_global, 1);
    if (!is_OK_global) throw std::runtime_error("Importer correctness test failed.");

    // Cleanup
    delete Aopt;
#endif  // defined(HAVE_MUELU_EPETRA)
  }
}

// =========================================================================
// =========================================================================
// =========================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor& clp, Xpetra::UnderlyingLib& lib, int argc, char* argv[]) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // =========================================================================
  // Convenient definitions
  // =========================================================================
  typedef Teuchos::ScalarTraits<SC> STS;
  SC one = STS::one();
  typedef typename STS::coordinateType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  // =========================================================================
  // Parameters initialization
  // =========================================================================
  GO nx = 100, ny = 100, nz = 100;
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D");  // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);                                       // manage parameters of Xpetra

  std::string xmlFileName = "import.xml";
  clp.setOption("xml", &xmlFileName, "read parameters from a file");
  bool printTimings = true;
  clp.setOption("timings", "notimings", &printTimings, "print timings to screen");
  std::string timingsFormat = "table-fixed";
  clp.setOption("time-format", &timingsFormat, "timings format (table-fixed | table-scientific | yaml)");
  int numImports = 100;
  clp.setOption("numImport", &numImports, "#times to test");
  int MM_TAFC_OptCoreCnt = 3000;
  clp.setOption("MM_TAFC_OptimizationCoreCount", &MM_TAFC_OptCoreCnt, "Num Cores above which Optimized MatrixMatrix transferAndFillComplete is used");

  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
  }

  ParameterList paramList;
  Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<ParameterList>(&paramList), *comm);

  ParameterList& mmlist = paramList.sublist("matrixmatrix: kernel params", false);
  int commandcc         = mmlist.get("MM_TAFC_OptimizationCoreCount", MM_TAFC_OptCoreCnt);
  commandcc             = paramList.get("MM_TAFC_OptimizationCoreCount", commandcc);
  paramList.remove("MM_TAFC_OptimizationCoreCount", false);
  mmlist.set("MM_TAFC_OptimizationCoreCount", commandcc);

  bool isDriver = paramList.isSublist("Run1");
  if (isDriver) {
    // update galeriParameters with the values from the XML file
    ParameterList& realParams = galeriParameters.GetParameterList();

    for (ParameterList::ConstIterator it = realParams.begin(); it != realParams.end(); it++) {
      const std::string& name = realParams.name(it);
      if (paramList.isParameter(name))
        realParams.setEntry(name, paramList.getEntry(name));
    }
  }

  // Retrieve matrix parameters (they may have been changed on the command line)
  // [for instance, if we changed matrix type from 2D to 3D we need to update nz]
  ParameterList galeriList = galeriParameters.GetParameterList();

  // =========================================================================
  // Problem construction
  // =========================================================================
  std::ostringstream galeriStream;
  comm->barrier();

  RCP<Matrix> A;
  RCP<const Map> map;
  RCP<RealValuedMultiVector> coordinates;
  typedef typename RealValuedMultiVector::scalar_type Real;
  RCP<MultiVector> nullspace;
  std::string matrixType = galeriParameters.GetMatrixType();

  RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: S - Global Time")));
  {
    auto tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 1 - Matrix Build")));

    galeriStream << "========================================================\n"
                 << xpetraParameters << galeriParameters;

    // Galeri will attempt to create a square-as-possible distribution of subdomains di, e.g.,
    //                                 d1  d2  d3
    //                                 d4  d5  d6
    //                                 d7  d8  d9
    //                                 d10 d11 d12
    // A perfect distribution is only possible when the #processors is a perfect square.
    // This *will* result in "strip" distribution if the #processors is a prime number or if the factors are very different in
    // size. For example, np=14 will give a 7-by-2 distribution.
    // If you don't want Galeri to do this, specify mx or my on the galeriList.

    // Create map and coordinates
    // In the future, we hope to be able to first create a Galeri problem, and then request map and coordinates from it
    // At the moment, however, things are fragile as we hope that the Problem uses same map and coordinates inside
    if (matrixType == "Laplace1D") {
      map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian1D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<Real, LO, GO, Map, RealValuedMultiVector>("1D", map, galeriList);

    } else if (matrixType == "Laplace2D" || matrixType == "Star2D" ||
               matrixType == "BigStar2D" || matrixType == "Elasticity2D") {
      map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian2D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<Real, LO, GO, Map, RealValuedMultiVector>("2D", map, galeriList);

    } else if (matrixType == "Laplace3D" || matrixType == "Brick3D" || matrixType == "Elasticity3D") {
      map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian3D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<Real, LO, GO, Map, RealValuedMultiVector>("3D", map, galeriList);
    }

    // Expand map to do multiple DOF per node for block problems
    if (matrixType == "Elasticity2D")
      map = Xpetra::MapFactory<LO, GO, Node>::Build(map, 2);
    if (matrixType == "Elasticity3D")
      map = Xpetra::MapFactory<LO, GO, Node>::Build(map, 3);

    galeriStream << "Processor subdomains in x direction: " << galeriList.get<GO>("mx") << std::endl
                 << "Processor subdomains in y direction: " << galeriList.get<GO>("my") << std::endl
                 << "Processor subdomains in z direction: " << galeriList.get<GO>("mz") << std::endl
                 << "========================================================" << std::endl;

    if (matrixType == "Elasticity2D" || matrixType == "Elasticity3D") {
      // Our default test case for elasticity: all boundaries of a square/cube have Neumann b.c. except left which has Dirichlet
      galeriList.set("right boundary", "Neumann");
      galeriList.set("bottom boundary", "Neumann");
      galeriList.set("top boundary", "Neumann");
      galeriList.set("front boundary", "Neumann");
      galeriList.set("back boundary", "Neumann");
    }

    RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(galeriParameters.GetMatrixType(), map, galeriList);
    A = Pr->BuildMatrix();

    if (matrixType == "Elasticity2D" ||
        matrixType == "Elasticity3D") {
      nullspace = Pr->BuildNullspace();
      A->SetFixedBlockSize((galeriParameters.GetMatrixType() == "Elasticity2D") ? 2 : 3);
    }

    comm->barrier();
    tm = Teuchos::null;
  }

  galeriStream << "Galeri complete.\n========================================================" << std::endl;

  int numReruns = 1;
  if (paramList.isParameter("number of reruns"))
    numReruns = paramList.get<int>("number of reruns");

  const bool mustAlreadyExist = true;
  for (int rerunCount = 1; rerunCount <= numReruns; rerunCount++) {
    ParameterList mueluList, runList;

    bool stop = false;
    if (isDriver) {
      runList   = paramList.sublist("Run1", mustAlreadyExist);
      mueluList = runList.sublist("MueLu", mustAlreadyExist);
    } else {
      mueluList = paramList;
      stop      = true;
    }

    if (nullspace.is_null()) {
      int blkSize = 1;
      if (mueluList.isSublist("Matrix")) {
        // Factory style parameter list
        const Teuchos::ParameterList& operatorList = paramList.sublist("Matrix");
        if (operatorList.isParameter("PDE equations"))
          blkSize = operatorList.get<int>("PDE equations");

      } else if (paramList.isParameter("number of equations")) {
        // Easy style parameter list
        blkSize = paramList.get<int>("number of equations");
      }

      nullspace = MultiVectorFactory::Build(map, blkSize);
      for (int i = 0; i < blkSize; i++) {
        RCP<const Map> domainMap = A->getDomainMap();
        GO indexBase             = domainMap->getIndexBase();

        ArrayRCP<SC> nsData = nullspace->getDataNonConst(i);
        for (int j = 0; j < nsData.size(); j++) {
          GO GID = domainMap->getGlobalElement(j) - indexBase;

          if ((GID - i) % blkSize == 0)
            nsData[j] = one;
        }
      }
    }

    int runCount = 1;
    do {
      int savedOut    = -1;
      FILE* openedOut = NULL;
      if (isDriver) {
        if (runList.isParameter("filename")) {
          // Redirect all output into a filename We have to redirect all output,
          // including printf's, therefore we cannot simply replace C++ cout
          // buffers, and have to use heavy machinary (dup2)
          std::string filename = runList.get<std::string>("filename");
          if (numReruns > 1)
            filename += "_run" + MueLu::toString(rerunCount);
          filename += (lib == Xpetra::UseEpetra ? ".epetra" : ".tpetra");

          savedOut  = dup(STDOUT_FILENO);
          openedOut = fopen(filename.c_str(), "w");
          dup2(fileno(openedOut), STDOUT_FILENO);
        }
      }

      // Instead of checking each time for rank, create a rank 0 stream
      RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      Teuchos::FancyOStream& out       = *fancy;
      out.setOutputToRootOnly(0);

      out << galeriStream.str();

      // =========================================================================
      // Preconditioner construction
      // =========================================================================
      comm->barrier();

      RCP<Hierarchy> H;
      {
        auto MueLuSU_D2 = TimeMonitor(*TimeMonitor::getNewTimer("Driver: 2 - MueLu Setup"));

        A->SetMaxEigenvalueEstimate(-one);
        if (lib == Xpetra::UseEpetra) {
          mueluList.set("use kokkos refactor", false);
        }
        Teuchos::ParameterList& userParamList = mueluList.sublist("user data");
        userParamList.set<RCP<RealValuedMultiVector> >("Coordinates", coordinates);
        H = MueLu::CreateXpetraPreconditioner(A, mueluList);
        comm->barrier();
      }

      // =========================================================================
      // Grab useful pieces
      // =========================================================================

      RCP<Matrix> P;
      {
        auto MueLuES_D3 = TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3 - Extract Stuff"));
        H->GetLevel(1)->Get("P", P);
        comm->barrier();
      }

      RCP<Teuchos::ParameterList> opt_list = rcp(new Teuchos::ParameterList);
      opt_list->set("Timer Label", "OptTAFC");

      for (int i = 0; i < numImports; i++) {
        // =========================================================================
        // Optimized transfer & fill complete loop for P_1
        // =========================================================================
        auto D4 = TimeMonitor(*TimeMonitor::getNewTimer("Driver: 4 - TransferAndFillComplete"));
        TestTransfer(A, P);
        comm->barrier();
      }

      globalTimeMonitor = Teuchos::null;

      // =========================================================================
      // Timing stuff
      // =========================================================================
      if (printTimings) {
        RCP<ParameterList> reportParams = rcp(new ParameterList);
        if (timingsFormat == "yaml") {
          reportParams->set("Report format", "YAML");  // "Table" or "YAML"
          reportParams->set("YAML style", "compact");  // "spacious" or "compact"
        }
        reportParams->set("How to merge timer sets", "Union");
        reportParams->set("alwaysWriteLocal", false);
        reportParams->set("writeGlobalStats", true);
        reportParams->set("writeZeroTimers", false);
        // FIXME: no "ignoreZeroTimers"

        const std::string filter = "";

        std::ios_base::fmtflags ff(out.flags());
        if (timingsFormat == "table-fixed")
          out << std::fixed;
        else
          out << std::scientific;
        TimeMonitor::report(comm.ptr(), out, filter, reportParams);
        out << std::setiosflags(ff);
      }

      TimeMonitor::clearCounters();

      if (isDriver) {
        if (openedOut != NULL) {
          TEUCHOS_ASSERT(savedOut >= 0);
          dup2(savedOut, STDOUT_FILENO);
          fclose(openedOut);
          openedOut = NULL;
        }
        try {
          runList   = paramList.sublist("Run" + MueLu::toString(++runCount), mustAlreadyExist);
          mueluList = runList.sublist("MueLu", mustAlreadyExist);
        } catch (Teuchos::Exceptions::InvalidParameterName&) {
          stop = true;
        }
      }

    } while (!stop);
  }

  comm->barrier();

  return EXIT_SUCCESS;
}

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char* argv[]) {
  auto val = Automatic_Test_ETI(argc, argv);
  return val;
}
