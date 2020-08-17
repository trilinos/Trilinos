/*
// @HEADER
// ***********************************************************************
//
//               Epetra: Linear Algebra Services Package
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
// ************************************************************************
// @HEADER
*/

// Creates vectors with different maps; tests results of export into them.
// Bug7758 indicates that incorrect results of export are produced when
// the source map is NOT a superset of the target map.
// Same tests as tpetra/core/test/MultiVector/Bug7758.cpp.

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"

#ifdef EPETRA_MPI
#include <mpi.h>
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

//////////////////////////////////////////////////////////////////////////////
int DefaultToDefaultEpetra(const Epetra_Comm &comm)
{
  int me = comm.MyPID();
  int np = comm.NumProc();
  int ierr = 0;

  using vector_t = Epetra_Vector;
  using map_t = Epetra_Map;

  const int nGlobalEntries = 8 * np;
  const double tgtScalar = 100. * (me+1);
  const double srcScalar = 2.;
  std::vector<int> myEntries(nGlobalEntries); 

  // Default one-to-one linear block map in Trilinos

  map_t defaultMap(nGlobalEntries, 0, comm);

  std::cout << me << " DEFAULT MAP" << std::endl;
  defaultMap.Print(std::cout);

  // Create vectors; see what the result is with CombineMode=ADD

  vector_t defaultVecTgt(defaultMap);
  defaultVecTgt.PutScalar(tgtScalar);

  vector_t defaultVecSrc(defaultMap);
  defaultVecSrc.PutScalar(srcScalar);

  // Export Default-to-default:  should be a copy of src to tgt

  Epetra_Export defaultToDefault(defaultMap, defaultMap);
  defaultVecTgt.Export(defaultVecSrc, defaultToDefault, Add);

  std::cout << me << " DEFAULT TO DEFAULT " << std::endl;
  defaultVecTgt.Print(std::cout);

  // Check result; all vector entries should be srcScalar
  for (int i = 0; i < defaultVecTgt.MyLength(); i++)
    if (defaultVecTgt[i] != srcScalar) ierr++;

  if (ierr > 0) 
    std::cout << "TEST FAILED:  DEFAULT-TO-DEFAULT-EPETRA TEST HAD " << ierr 
              << " FAILURES ON RANK " << me << std::endl;

  int gerr;
  comm.SumAll(&ierr, &gerr, 1);

  return gerr;
}

//////////////////////////////////////////////////////////////////////////////
int CyclicToDefaultEpetra(const Epetra_Comm &comm)
{
  int me = comm.MyPID();
  int np = comm.NumProc();
  int ierr = 0;

  using vector_t = Epetra_Vector;
  using map_t = Epetra_Map;

  const int nGlobalEntries = 8 * np;
  const double tgtScalar = 100. * (me+1);
  const double srcScalar = 2.;
  std::vector<int> myEntries(nGlobalEntries); 

  // Default one-to-one linear block map in Trilinos

  map_t defaultMap(nGlobalEntries, 0, comm);

  std::cout << me << " DEFAULT MAP" << std::endl;
  defaultMap.Print(std::cout);

  // One-to-one cyclic map:  deal out entries like cards

  int nMyEntries = 0;
  for (int i = 0; i < nGlobalEntries; i++) {
    if (i % np == me) {
      myEntries[nMyEntries++] = i;
    }
  }

  int dummy = -1;
  map_t cyclicMap(dummy, nMyEntries, &myEntries[0], 0, comm);

  std::cout << me << " CYCLIC MAP" << std::endl;
  cyclicMap.Print(std::cout);

  // Create vectors; see what the result is with CombineMode=ADD

  vector_t defaultVecTgt(defaultMap);
  defaultVecTgt.PutScalar(tgtScalar);

  vector_t cyclicVecSrc(cyclicMap);
  cyclicVecSrc.PutScalar(srcScalar);

  // Export Cyclic-to-default

  Epetra_Export cyclicToDefault(cyclicMap, defaultMap);
  defaultVecTgt.Export(cyclicVecSrc, cyclicToDefault, Add);

  std::cout << me << " CYCLIC TO DEFAULT " << std::endl;
  defaultVecTgt.Print(std::cout);

  // Check result; all vector entries should be srcScalar

  for (int i = 0; i < defaultVecTgt.MyLength(); i++)
    if (defaultVecTgt[i] != srcScalar) ierr++;
  if (ierr > 0) 
    std::cout << "TEST FAILED:  CYCLIC-TO-DEFAULT-EPETRA TEST HAD " << ierr 
              << " FAILURES ON RANK " << me << std::endl;

  int gerr;
  comm.SumAll(&ierr, &gerr, 1);

  return gerr;
}

//////////////////////////////////////////////////////////////////////////////
int OverlapToDefaultEpetra(const Epetra_Comm &comm)
{
  int me = comm.MyPID();
  int np = comm.NumProc();
  int ierr = 0;

  using vector_t = Epetra_Vector;
  using map_t = Epetra_Map;

  if (np > 1) {  // Need more than one proc to avoid duplicate entries in maps

    const int nGlobalEntries = 8 * np;
    const double tgtScalar = 100. * (me+1);
    const double srcScalar = 2.;
    std::vector<int> myEntries(nGlobalEntries); 

    // Default one-to-one linear block map in Trilinos

    map_t defaultMap(nGlobalEntries, 0, comm);

    std::cout << me << " DEFAULT MAP" << std::endl;
    defaultMap.Print(std::cout);

    // Overlap map; some entries are stored on two procs
    int nMyEntries = 0;
    for (int i = 0; i < defaultMap.NumMyElements()/2; i++) {
      myEntries[nMyEntries++] = defaultMap.GID(i);
    }
    for (int i = 0; i < defaultMap.NumMyElements(); i++) {
      myEntries[nMyEntries++] =
        (defaultMap.MaxMyGID() + 1 + i) % nGlobalEntries;
    }
  
    int dummy = -1;
    map_t overlapMap(dummy, nMyEntries, &myEntries[0], 0, comm);
  
    std::cout << me << " OVERLAP MAP" << std::endl;
    overlapMap.Print(std::cout);

    // Create vectors; see what the result is with CombineMode=ADD

    vector_t defaultVecTgt(defaultMap);
    defaultVecTgt.PutScalar(tgtScalar);

    vector_t overlapVecSrc(overlapMap);
    overlapVecSrc.PutScalar(srcScalar);

    // Export Overlap-to-default

    Epetra_Export overlapToDefault(overlapMap, defaultMap);
    defaultVecTgt.Export(overlapVecSrc, overlapToDefault, Add);

    std::cout << me << " OVERLAP TO DEFAULT " << std::endl;
    defaultVecTgt.Print(std::cout);

    for (int i = 0; i < defaultVecTgt.MyLength()/2; i++)
      if (defaultVecTgt[i] != 2 * srcScalar) ierr++;  // overlapped
    for (int i = defaultVecTgt.MyLength()/2;
             i < defaultVecTgt.MyLength(); i++)
      if (defaultVecTgt[i] != srcScalar) ierr++;  // not overlapped
    if (ierr > 0) 
      std::cout << "TEST FAILED:  OVERLAP-TO-DEFAULT-EPETRA TEST HAD " << ierr 
                << " FAILURES ON RANK " << me << std::endl;
  }

  int gerr;
  comm.SumAll(&ierr, &gerr, 1);

  return gerr;
}

//////////////////////////////////////////////////////////////////////////////
int SupersetToDefaultEpetra(const Epetra_Comm &comm)
{
  int me = comm.MyPID();
  int np = comm.NumProc();
  int ierr = 0;

  using vector_t = Epetra_Vector;
  using map_t = Epetra_Map;

  if (np > 1) {  // Need more than one proc to avoid duplicate entries in maps

    const int nGlobalEntries = 8 * np;
    const double tgtScalar = 100. * (me+1);
    const double srcScalar = 2.;
    std::vector<int> myEntries(nGlobalEntries); 

    // Default one-to-one linear block map in Trilinos

    map_t defaultMap(nGlobalEntries, 0, comm);

    std::cout << me << " DEFAULT MAP" << std::endl;
    defaultMap.Print(std::cout);

    // Superset map; some entries are stored on two procs
    int nMyEntries = 0;
    for (int i = 0; i < defaultMap.NumMyElements(); i++) {
      myEntries[nMyEntries++] = defaultMap.GID(i);
    }
    for (int i = 0; i < defaultMap.NumMyElements()/2; i++) {
      myEntries[nMyEntries++] =
        (defaultMap.MaxMyGID() + 1 + i) % nGlobalEntries;
    }
  
    int dummy = -1;
    map_t supersetMap(dummy, nMyEntries, &myEntries[0], 0, comm);
  
    std::cout << me << " SUPERSET MAP" << std::endl;
    supersetMap.Print(std::cout);

    // Create vectors; see what the result is with CombineMode=ADD

    vector_t defaultVecTgt(defaultMap);
    defaultVecTgt.PutScalar(tgtScalar);

    vector_t supersetVecSrc(supersetMap);
    supersetVecSrc.PutScalar(srcScalar);

    // Export Superset-to-default

    Epetra_Export supersetToDefault(supersetMap, defaultMap);
    defaultVecTgt.Export(supersetVecSrc, supersetToDefault, Add);

    std::cout << me << " SUPERSET TO DEFAULT " << std::endl;
    defaultVecTgt.Print(std::cout);

    for (int i = 0; i < defaultVecTgt.MyLength()/2; i++)
      if (defaultVecTgt[i] != 2 * srcScalar) ierr++;  // overlapped
    for (int i = defaultVecTgt.MyLength()/2;
             i < defaultVecTgt.MyLength(); i++)
      if (defaultVecTgt[i] != srcScalar) ierr++;  // not overlapped
    if (ierr > 0) 
      std::cout << "TEST FAILED:  SUPERSET-TO-DEFAULT-EPETRA TEST HAD " << ierr 
                << " FAILURES ON RANK " << me << std::endl;
  }

  int gerr;
  comm.SumAll(&ierr, &gerr, 1);

  return gerr;
}

//////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv) 
{
#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int gerr = 0;

  gerr += DefaultToDefaultEpetra(Comm);
  gerr += CyclicToDefaultEpetra(Comm);
  gerr += OverlapToDefaultEpetra(Comm);
  gerr += SupersetToDefaultEpetra(Comm);

  if (Comm.MyPID() == 0) {
    if (gerr > 0) std::cout << "TEST FAILED" << std::endl;
    else std::cout << "TEST PASSED" << std::endl;
  }

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return gerr;

}
