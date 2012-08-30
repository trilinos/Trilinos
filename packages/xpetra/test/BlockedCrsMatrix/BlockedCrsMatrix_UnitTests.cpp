// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
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
/*
 * BlockedCrsMatrix_UnitTests.cpp
 *
 *  Created on: Aug 22, 2011
 *      Author: wiesner
 */

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#  include "mpi.h"
#endif
#  include "Epetra_SerialComm.h"

#include <Xpetra_ConfigDefs.hpp>

#ifdef HAVE_XPETRA_EPETRAEXT
// EpetraExt
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#endif

#include <Xpetra_DefaultPlatform.hpp>
#include <Teuchos_as.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#ifdef HAVE_XPETRA_TPETRA
#include <Xpetra_TpetraCrsMatrix.hpp>
#endif
#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraCrsMatrix.hpp>
#endif
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_Exceptions.hpp>

#include <MueLu_Utilities.hpp> //TODO: Xpetra tests should not use MueLu

namespace {

  using Teuchos::Array;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::arcp;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::Tuple;
  using Teuchos::tuple;
  using std::sort;
  using std::find;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::Comm;

  using Xpetra::DefaultPlatform;
  using Xpetra::Matrix;
  using Xpetra::CrsMatrix;
#ifdef HAVE_XPETRA_TPETRA
  using Xpetra::TpetraCrsMatrix; //TMP
#endif
  using Xpetra::Map;

  using Xpetra::viewLabel_t;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

  //////////////////////////////////////////////////////////////////////////
  // EPETRA helper functions
  Teuchos::RCP<Epetra_Map> SplitMap(const Epetra_Map& Amap,
                                    const Epetra_Map& Agiven)
  {
    const Epetra_Comm& Comm = Amap.Comm();
    const Epetra_Map&  Ag = Agiven;

    int count=0;
    std::vector<int> myaugids(Amap.NumMyElements());
    for (int i=0; i<Amap.NumMyElements(); ++i)
      {
        const int gid = Amap.GID(i);
        if (Ag.MyGID(gid)) continue;
        myaugids[count] = gid;
        ++count;
      }
    myaugids.resize(count);
    int gcount;
    Comm.SumAll(&count,&gcount,1);
    Teuchos::RCP<Epetra_Map> Aunknown = Teuchos::rcp(new Epetra_Map(gcount,count,&myaugids[0],0,Comm));

    return Aunknown;
  }

  Teuchos::RCP<Epetra_Map> CreateMap(const std::set<int>& gids, const Epetra_Comm& comm)
  {
    std::vector<int> mapvec;
    mapvec.reserve(gids.size());
    mapvec.assign(gids.begin(), gids.end());
    Teuchos::RCP<Epetra_Map> map =
      Teuchos::rcp(new Epetra_Map(-1,
                                  mapvec.size(),
                                  &mapvec[0],
                                  0,
                                  comm));
    mapvec.clear();
    return map;
  }

  Teuchos::RCP<Epetra_Map> MergeMaps(const std::vector<Teuchos::RCP<const Epetra_Map> >& maps)
  {
    if (maps.size()==0)
      cout << "no maps to merge" << endl;
    for (unsigned i=0; i<maps.size(); ++i)
      {
        if (maps[i]==Teuchos::null)
          cout << "can not merge extractor with null maps" << endl;
        if (not maps[i]->UniqueGIDs())
          cout << "map " << i <<  " not unique" << endl;
      }
    std::set<int> mapentries;
    for (unsigned i=0; i<maps.size(); ++i)
      {
        const Epetra_Map& map = *maps[i];
        std::copy(map.MyGlobalElements(),
                  map.MyGlobalElements()+map.NumMyElements(),
                  std::inserter(mapentries,mapentries.begin()));
      }
    return CreateMap(mapentries, maps[0]->Comm());
  }


  bool SplitMatrix2x2(Teuchos::RCP<const Epetra_CrsMatrix> A,
                      const Epetra_Map& A11rowmap,
                      const Epetra_Map& A22rowmap,
                      Teuchos::RCP<Epetra_CrsMatrix>& A11,
                      Teuchos::RCP<Epetra_CrsMatrix>& A12,
                      Teuchos::RCP<Epetra_CrsMatrix>& A21,
                      Teuchos::RCP<Epetra_CrsMatrix>& A22)
  {
    if (A==Teuchos::null)
      {
        cout << "ERROR: SplitMatrix2x2: A==null on entry" << endl;
        return false;
      }

    const Epetra_Comm& Comm   = A->Comm();
    const Epetra_Map&  A22map = A22rowmap;
    const Epetra_Map&  A11map = A11rowmap;

    //----------------------------- create a parallel redundant map of A22map
    std::map<int,int> a22gmap;
    {
      std::vector<int> a22global(A22map.NumGlobalElements());
      int count=0;
      for (int proc=0; proc<Comm.NumProc(); ++proc)
        {
          int length = 0;
          if (proc==Comm.MyPID())
            {
              for (int i=0; i<A22map.NumMyElements(); ++i)
                {
                  a22global[count+length] = A22map.GID(i);
                  ++length;
                }
            }
          Comm.Broadcast(&length,1,proc);
          Comm.Broadcast(&a22global[count],length,proc);
          count += length;
        }
      if (count != A22map.NumGlobalElements())
        {
          cout << "ERROR SplitMatrix2x2: mismatch in dimensions" << endl;
          return false;
        }

      // create the map
      for (int i=0; i<count; ++i)
        a22gmap[a22global[i]] = 1;
      a22global.clear();
    }

    //--------------------------------------------------- create matrix A22
    A22 = Teuchos::rcp(new Epetra_CrsMatrix(Copy,A22map,100));
    {
      std::vector<int>    a22gcindices(100);
      std::vector<double> a22values(100);
      for (int i=0; i<A->NumMyRows(); ++i)
        {
          const int grid = A->GRID(i);
          if (A22map.MyGID(grid)==false)
            continue;
          int     numentries;
          double* values;
          int*    cindices;
          int err = A->ExtractMyRowView(i,numentries,values,cindices);
          if (err)
            {
              cout << "ERROR: SplitMatrix2x2: A->ExtractMyRowView returned " << err << endl;
              return false;
            }

          if (numentries>(int)a22gcindices.size())
            {
              a22gcindices.resize(numentries);
              a22values.resize(numentries);
            }
          int count=0;
          for (int j=0; j<numentries; ++j)
            {
              const int gcid = A->ColMap().GID(cindices[j]);
              // see whether we have gcid in a22gmap
              std::map<int,int>::iterator curr = a22gmap.find(gcid);
              if (curr==a22gmap.end()) continue;
              //cout << gcid << " ";
              a22gcindices[count] = gcid;
              a22values[count]    = values[j];
              ++count;
            }
          //cout << endl; fflush(stdout);
          // add this filtered row to A22
          err = A22->InsertGlobalValues(grid,count,&a22values[0],&a22gcindices[0]);
          if (err<0)
            {
              cout << "ERROR: SplitMatrix2x2: A->InsertGlobalValues returned " << err << endl;
              return false;
            }

        } //for (int i=0; i<A->NumMyRows(); ++i)
      a22gcindices.clear();
      a22values.clear();
    }
    A22->FillComplete();
    A22->OptimizeStorage();

    //----------------------------------------------------- create matrix A11
    A11 = Teuchos::rcp(new Epetra_CrsMatrix(Copy,A11map,100));
    {
      std::vector<int>    a11gcindices(100);
      std::vector<double> a11values(100);
      for (int i=0; i<A->NumMyRows(); ++i)
        {
          const int grid = A->GRID(i);
          if (A11map.MyGID(grid)==false) continue;
          int     numentries;
          double* values;
          int*    cindices;
          int err = A->ExtractMyRowView(i,numentries,values,cindices);
          if (err)
            {
              cout << "ERROR: SplitMatrix2x2: A->ExtractMyRowView returned " << err << endl;
              return false;
            }

          if (numentries>(int)a11gcindices.size())
            {
              a11gcindices.resize(numentries);
              a11values.resize(numentries);
            }
          int count=0;
          for (int j=0; j<numentries; ++j)
            {
              const int gcid = A->ColMap().GID(cindices[j]);
              // see whether we have gcid as part of a22gmap
              std::map<int,int>::iterator curr = a22gmap.find(gcid);
              if (curr!=a22gmap.end()) continue;
              a11gcindices[count] = gcid;
              a11values[count] = values[j];
              ++count;
            }
          err = A11->InsertGlobalValues(grid,count,&a11values[0],&a11gcindices[0]);
          if (err<0)
            {
              cout << "ERROR: SplitMatrix2x2: A->InsertGlobalValues returned " << err << endl;
              return false;
            }

        } // for (int i=0; i<A->NumMyRows(); ++i)
      a11gcindices.clear();
      a11values.clear();
    }
    A11->FillComplete();
    A11->OptimizeStorage();

    //---------------------------------------------------- create matrix A12
    A12 = Teuchos::rcp(new Epetra_CrsMatrix(Copy,A11map,100));
    {
      std::vector<int>    a12gcindices(100);
      std::vector<double> a12values(100);
      for (int i=0; i<A->NumMyRows(); ++i)
        {
          const int grid = A->GRID(i);
          if (A11map.MyGID(grid)==false) continue;
          int     numentries;
          double* values;
          int*    cindices;
          int err = A->ExtractMyRowView(i,numentries,values,cindices);
          if (err)
            {
              cout << "ERROR: SplitMatrix2x2: A->ExtractMyRowView returned " << err << endl;
              return false;
            }

          if (numentries>(int)a12gcindices.size())
            {
              a12gcindices.resize(numentries);
              a12values.resize(numentries);
            }
          int count=0;
          for (int j=0; j<numentries; ++j)
            {
              const int gcid = A->ColMap().GID(cindices[j]);
              // see whether we have gcid as part of a22gmap
              std::map<int,int>::iterator curr = a22gmap.find(gcid);
              if (curr==a22gmap.end()) continue;
              a12gcindices[count] = gcid;
              a12values[count] = values[j];
              ++count;
            }
          err = A12->InsertGlobalValues(grid,count,&a12values[0],&a12gcindices[0]);
          if (err<0)
            {
              cout << "ERROR: SplitMatrix2x2: A->InsertGlobalValues returned " << err << endl;
              return false;
            }

        } // for (int i=0; i<A->NumMyRows(); ++i)
      a12values.clear();
      a12gcindices.clear();
    }
    A12->FillComplete(A22map,A11map);
    A12->OptimizeStorage();

    //----------------------------------------------------------- create A21
    A21 = Teuchos::rcp(new Epetra_CrsMatrix(Copy,A22map,100));
    {
      std::vector<int>    a21gcindices(100);
      std::vector<double> a21values(100);
      for (int i=0; i<A->NumMyRows(); ++i)
        {
          const int grid = A->GRID(i);
          if (A22map.MyGID(grid)==false) continue;
          int     numentries;
          double* values;
          int*    cindices;
          int err = A->ExtractMyRowView(i,numentries,values,cindices);
          if (err)
            {
              cout << "ERROR: SplitMatrix2x2: A->ExtractMyRowView returned " << err << endl;
              return false;
            }

          if (numentries>(int)a21gcindices.size())
            {
              a21gcindices.resize(numentries);
              a21values.resize(numentries);
            }
          int count=0;
          for (int j=0; j<numentries; ++j)
            {
              const int gcid = A->ColMap().GID(cindices[j]);
              // see whether we have gcid as part of a22gmap
              std::map<int,int>::iterator curr = a22gmap.find(gcid);
              if (curr!=a22gmap.end()) continue;
              a21gcindices[count] = gcid;
              a21values[count] = values[j];
              ++count;
            }
          err = A21->InsertGlobalValues(grid,count,&a21values[0],&a21gcindices[0]);
          if (err<0)
            {
              cout << "ERROR: SplitMatrix2x2: A->InsertGlobalValues returned " << err << endl;
              return false;
            }

        } // for (int i=0; i<A->NumMyRows(); ++i)
      a21values.clear();
      a21gcindices.clear();
    }
    A21->FillComplete(A11map,A22map);
    A21->OptimizeStorage();

    //-------------------------------------------------------------- tidy up
    a22gmap.clear();
    return true;
  }

  bool SplitMatrix2x1(Teuchos::RCP<const Epetra_CrsMatrix> A,
                      const Epetra_Map& A1rowmap,
                      const Epetra_Map& A2rowmap,
                      const Epetra_Map& domainmap,
                      Teuchos::RCP<Epetra_CrsMatrix>& A1,
                      Teuchos::RCP<Epetra_CrsMatrix>& A2)
  {
    if (A==Teuchos::null)
      {
        cout << "ERROR: SplitMatrix2x2: A==null on entry" << endl;
        return false;
      }

    const Epetra_Map&  A1map = A1rowmap;
    const Epetra_Map&  A2map = A2rowmap;

    //----------------------------- create a parallel redundant map of A22map
    /*std::map<int,int> a22gmap;
      {
      std::vector<int> a22global(A22map.NumGlobalElements());
      int count=0;
      for (int proc=0; proc<Comm.NumProc(); ++proc)
      {
      int length = 0;
      if (proc==Comm.MyPID())
      {
      for (int i=0; i<A22map.NumMyElements(); ++i)
      {
      a22global[count+length] = A22map.GID(i);
      ++length;
      }
      }
      Comm.Broadcast(&length,1,proc);
      Comm.Broadcast(&a22global[count],length,proc);
      count += length;
      }
      if (count != A22map.NumGlobalElements())
      {
      cout << "ERROR SplitMatrix2x2: mismatch in dimensions" << endl;
      return false;
      }

      // create the map
      for (int i=0; i<count; ++i)
      a22gmap[a22global[i]] = 1;
      a22global.clear();
      }*/

    //--------------------------------------------------- create matrix A22
    A2 = Teuchos::rcp(new Epetra_CrsMatrix(Copy,A2map,100));
    {
      std::vector<int>    a22gcindices(100);
      std::vector<double> a22values(100);
      for (int i=0; i<A->NumMyRows(); ++i)
        {
          const int grid = A->GRID(i);
          if (A2map.MyGID(grid)==false)
            continue;
          int     numentries;
          double* values;
          int*    cindices;
          int err = A->ExtractMyRowView(i,numentries,values,cindices);
          if (err)
            {
              cout << "ERROR: SplitMatrix2x1: A->ExtractMyRowView returned " << err << endl;
              return false;
            }

          if (numentries>(int)a22gcindices.size())
            {
              a22gcindices.resize(numentries);
              a22values.resize(numentries);
            }
          int count=0;
          for (int j=0; j<numentries; ++j)
            {
              const int gcid = A->ColMap().GID(cindices[j]);
              a22gcindices[count] = gcid;
              a22values[count]    = values[j];
              ++count;
            }
          //cout << endl; fflush(stdout);
          // add this filtered row to A22
          err = A2->InsertGlobalValues(grid,count,&a22values[0],&a22gcindices[0]);
          if (err<0)
            {
              cout << "ERROR: SplitMatrix2x2: A->InsertGlobalValues returned " << err << endl;
              return false;
            }
          //cout << "PROC " << A->Comm().MyPID() << " InsertGlobalValues row " << grid << " entries " << count << endl;

        } //for (int i=0; i<A->NumMyRows(); ++i)
      a22gcindices.clear();
      a22values.clear();
    }
    A2->FillComplete(domainmap,A2map);
    A2->OptimizeStorage();

    //----------------------------------------------------- create matrix A11
    A1 = Teuchos::rcp(new Epetra_CrsMatrix(Copy,A1map,100));
    {
      std::vector<int>    a11gcindices(100);
      std::vector<double> a11values(100);
      for (int i=0; i<A->NumMyRows(); ++i)
        {
          const int grid = A->GRID(i);
          if (A1map.MyGID(grid)==false) continue;
          int     numentries;
          double* values;
          int*    cindices;
          int err = A->ExtractMyRowView(i,numentries,values,cindices);
          if (err)
            {
              cout << "ERROR: SplitMatrix2x2: A->ExtractMyRowView returned " << err << endl;
              return false;
            }

          if (numentries>(int)a11gcindices.size())
            {
              a11gcindices.resize(numentries);
              a11values.resize(numentries);
            }
          int count=0;
          for (int j=0; j<numentries; ++j)
            {
              const int gcid = A->ColMap().GID(cindices[j]);
              a11gcindices[count] = gcid;
              a11values[count] = values[j];
              ++count;
            }
          err = A1->InsertGlobalValues(grid,count,&a11values[0],&a11gcindices[0]);
          if (err<0)
            {
              cout << "ERROR: SplitMatrix2x2: A->InsertGlobalValues returned " << err << endl;
              return false;
            }

        } // for (int i=0; i<A->NumMyRows(); ++i)
      a11gcindices.clear();
      a11values.clear();
    }
    A1->FillComplete(domainmap,A1map);
    A1->OptimizeStorage();

    return true;
  }

  RCP<const Comm<int> > getDefaultComm()
  {
    if (testMpi) {
      return DefaultPlatform::getDefaultPlatform().getComm();
    }
    return rcp(new Teuchos::SerialComm<int>());
  }

  /////////////////////////////////////////////////////

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
                  "test-mpi", "test-serial", &testMpi,
                  "Test MPI (if available) or force test of serial.  In a serial build,"
                  " this option is ignored and a serial comm is always used." );
    clp.setOption(
                  "error-tol-slack", &errorTolSlack,
                  "Slack off of machine epsilon used to check test results" );
  }

  //
  // UNIT TESTS
  //


  /// simple test routine for the apply function of BlockedCrsMatrix
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockedCrsMatrix, EpetraApply, Scalar, LO, GO, Node )
  {
#ifdef HAVE_XPETRA_EPETRAEXT

#ifdef HAVE_MPI

    RCP<Epetra_Comm> Comm;
    if(testMpi)
#ifdef HAVE_MPI
      Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
      Comm = Teuchos::rcp(new Epetra_SerialComm);
#endif
    else
      Comm = Teuchos::rcp(new Epetra_SerialComm);

    // 1) load all matrices
    Epetra_Map pointmap(5148,0,*Comm);  // 5148  2

    // generate local maps for loading matrices
    std::vector<int> velgidvec; // global strided maps
    std::vector<int> pregidvec;
    std::vector<int> fullgidvec; // full global map
    for (int i=0; i<pointmap.NumMyElements(); i++)
      {
        // loop over all local ids in pointmap

        // get corresponding global id
        int gid = pointmap.GID(i);

        // store global strided gids
        velgidvec.push_back(3*gid);
        velgidvec.push_back(3*gid+1);
        pregidvec.push_back(3*gid+2);

        // gid for full map
        fullgidvec.push_back(3*gid);
        fullgidvec.push_back(3*gid+1);
        fullgidvec.push_back(3*gid+2);
      }

    // generate strided maps
    Teuchos::RCP<const Epetra_Map> velmap = Teuchos::rcp(new const Epetra_Map(-1, velgidvec.size(), &velgidvec[0], 0, *Comm));
    Teuchos::RCP<const Epetra_Map> premap = Teuchos::rcp(new const Epetra_Map(-1, pregidvec.size(), &pregidvec[0], 0, *Comm));

    // generate full map
    const Teuchos::RCP<const Epetra_Map> fullmap = Teuchos::rcp(new const Epetra_Map(-1, fullgidvec.size(), &fullgidvec[0], 0, *Comm));

    // read in matrices
    Epetra_CrsMatrix* ptrA = 0;
    Epetra_Vector*    ptrx = 0;
    Epetra_Vector*    ptrf = 0;
    EpetraExt::MatrixMarketFileToCrsMatrix("nsjac_test.mm",*fullmap,*fullmap,*fullmap,ptrA);
    EpetraExt::MatrixMarketFileToVector("nsrhs_test.mm",*fullmap,ptrf);
    EpetraExt::MatrixMarketFileToVector("nslhs_test.mm",*fullmap,ptrx);

    Teuchos::RCP<Epetra_CrsMatrix> fullA = Teuchos::rcp(ptrA);
    Teuchos::RCP<Epetra_Vector> b = Teuchos::rcp(ptrf);
    Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp(ptrx);

    // split fullA into A11,..., A22
    Teuchos::RCP<Epetra_CrsMatrix> A11;
    Teuchos::RCP<Epetra_CrsMatrix> A12;
    Teuchos::RCP<Epetra_CrsMatrix> A21;
    Teuchos::RCP<Epetra_CrsMatrix> A22;

    TEST_EQUALITY(SplitMatrix2x2(fullA,*velmap,*premap,A11,A12,A21,A22),true);

    // build Xpetra objects from Epetra_CrsMatrix objects
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xfuA = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(fullA));
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xA11 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A11));
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xA12 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A12));
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xA21 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A21));
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xA22 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A22));

    // build map extractor
    Teuchos::RCP<Xpetra::EpetraMap> xfullmap = Teuchos::rcp(new Xpetra::EpetraMap(fullmap));
    Teuchos::RCP<Xpetra::EpetraMap> xvelmap  = Teuchos::rcp(new Xpetra::EpetraMap(velmap ));
    Teuchos::RCP<Xpetra::EpetraMap> xpremap  = Teuchos::rcp(new Xpetra::EpetraMap(premap ));

    std::vector<Teuchos::RCP<const Xpetra::Map<LO,GO,Node> > > xmaps;
    xmaps.push_back(xvelmap);
    xmaps.push_back(xpremap);

    Teuchos::RCP<const Xpetra::MapExtractor<Scalar,LO,GO,Node> > map_extractor = Xpetra::MapExtractorFactory<Scalar,LO,GO>::Build(xfullmap,xmaps);

    // build blocked operator
    Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> > bOp = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar,LO,GO>(map_extractor,map_extractor,10));
    bOp->setMatrix(0,0,xA11);
    bOp->setMatrix(0,1,xA12);
    bOp->setMatrix(1,0,xA21);
    bOp->setMatrix(1,1,xA22);

    bOp->fillComplete();

    // build vector
    Teuchos::RCP<Xpetra::EpetraVector> xx  = Teuchos::rcp(new Xpetra::EpetraVector(x));
    Teuchos::RCP<Xpetra::Vector<Scalar,LO,GO> > result =  Xpetra::VectorFactory<Scalar,LO,GO>::Build(xfullmap ,true);

    // matrix vector product
    bOp->apply(*xx,*result,Teuchos::NO_TRANS);

    // check results
    Teuchos::RCP<Xpetra::Vector<Scalar,LO,GO,Node> > result2 =  Xpetra::VectorFactory<Scalar,LO,GO>::Build(xfullmap ,true);
    xfuA->apply(*xx,*result2);

    result2->update(1.0,*result,-1.0);

    //cout << "norm of difference " << result2->norm2() << endl;

    TEUCHOS_TEST_COMPARE(result2->norm2(), <, 1e-16, out, success);
#endif // MPI

#endif // XPETRA_EPETRA

  }

  /// simple test for matrix-matrix multiplication for two 2x2 blocked matrices
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockedCrsMatrix, EpetraMatrixMatrixMult, Scalar, LO, GO, Node ) //TODO: add template parameter <Node,...>
  {
#ifdef HAVE_XPETRA_EPETRAEXT

    RCP<Epetra_Comm> Comm;
    if(testMpi)
#ifdef HAVE_MPI
      Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
      Comm = Teuchos::rcp(new Epetra_SerialComm);
#endif
    else
      Comm = Teuchos::rcp(new Epetra_SerialComm);

    // 1) load all matrices
    Epetra_Map pointmap(5148,0,*Comm);  // 5148  2

    // generate local maps for loading matrices
    std::vector<int> velgidvec; // global strided maps
    std::vector<int> pregidvec;
    std::vector<int> fullgidvec; // full global map
    for (int i=0; i<pointmap.NumMyElements(); i++)
      {
        // loop over all local ids in pointmap

        // get corresponding global id
        int gid = pointmap.GID(i);

        // store global strided gids
        velgidvec.push_back(3*gid);
        velgidvec.push_back(3*gid+1);
        pregidvec.push_back(3*gid+2);

        // gid for full map
        fullgidvec.push_back(3*gid);
        fullgidvec.push_back(3*gid+1);
        fullgidvec.push_back(3*gid+2);
      }

    // generate strided maps
    Teuchos::RCP<const Epetra_Map> velmap = Teuchos::rcp(new const Epetra_Map(-1, velgidvec.size(), &velgidvec[0], 0, *Comm));
    Teuchos::RCP<const Epetra_Map> premap = Teuchos::rcp(new const Epetra_Map(-1, pregidvec.size(), &pregidvec[0], 0, *Comm));

    // generate full map
    const Teuchos::RCP<const Epetra_Map> fullmap = Teuchos::rcp(new const Epetra_Map(-1, fullgidvec.size(), &fullgidvec[0], 0, *Comm));

    // read in matrices
    Epetra_CrsMatrix* ptrA = 0;
    Epetra_Vector*    ptrx = 0;
    Epetra_Vector*    ptrf = 0;
    EpetraExt::MatrixMarketFileToCrsMatrix("nsjac_test.mm",*fullmap,*fullmap,*fullmap,ptrA);
    EpetraExt::MatrixMarketFileToVector("nsrhs_test.mm",*fullmap,ptrf);
    EpetraExt::MatrixMarketFileToVector("nslhs_test.mm",*fullmap,ptrx);

    Teuchos::RCP<Epetra_CrsMatrix> fullA = Teuchos::rcp(ptrA);
    Teuchos::RCP<Epetra_Vector> b = Teuchos::rcp(ptrf);
    Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp(ptrx);

    // split fullA into A11,..., A22
    Teuchos::RCP<Epetra_CrsMatrix> A11;
    Teuchos::RCP<Epetra_CrsMatrix> A12;
    Teuchos::RCP<Epetra_CrsMatrix> A21;
    Teuchos::RCP<Epetra_CrsMatrix> A22;

    TEST_EQUALITY(SplitMatrix2x2(fullA,*velmap,*premap,A11,A12,A21,A22),true);

    //////////////////////////////////////
    // build 1st block matrix

    // build Xpetra objects from Epetra_CrsMatrix objects
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xfuA = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(fullA)); //TODO: should use ALL the template parameters
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xA11 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A11));
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xA12 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A12));
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xA21 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A21));
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xA22 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A22));

    // build map extractor
    Teuchos::RCP<Xpetra::EpetraMap> xfullmap = Teuchos::rcp(new Xpetra::EpetraMap(fullmap));
    Teuchos::RCP<Xpetra::EpetraMap> xvelmap  = Teuchos::rcp(new Xpetra::EpetraMap(velmap ));
    Teuchos::RCP<Xpetra::EpetraMap> xpremap  = Teuchos::rcp(new Xpetra::EpetraMap(premap ));

    std::vector<Teuchos::RCP<const Xpetra::Map<LO,GO> > > xmaps;
    xmaps.push_back(xvelmap);
    xmaps.push_back(xpremap);

    Teuchos::RCP<const Xpetra::MapExtractor<Scalar,LO,GO> > map_extractor = Xpetra::MapExtractorFactory<Scalar,LO,GO>::Build(xfullmap,xmaps);

    // build blocked operator
    Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO> > bOp = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar,LO,GO>(map_extractor,map_extractor,10));
    bOp->setMatrix(0,0,xA11);
    bOp->setMatrix(0,1,xA12);
    bOp->setMatrix(1,0,xA21);
    bOp->setMatrix(1,1,xA22);

    bOp->fillComplete();
    //////////////////////////////


    //////////////////////////////////////
    // build 2nd block matrix

    // split fullA into A11,..., A22
    Teuchos::RCP<Epetra_CrsMatrix> fullA_2 = Teuchos::rcp(new Epetra_CrsMatrix(*fullA));
    Teuchos::RCP<Epetra_CrsMatrix> A11_2;
    Teuchos::RCP<Epetra_CrsMatrix> A12_2;
    Teuchos::RCP<Epetra_CrsMatrix> A21_2;
    Teuchos::RCP<Epetra_CrsMatrix> A22_2;

    TEST_EQUALITY(SplitMatrix2x2(fullA_2,*velmap,*premap,A11_2,A12_2,A21_2,A22_2),true);

    // build Xpetra objects from Epetra_CrsMatrix objects
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xfuA_2 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(fullA_2));
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xA11_2 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A11_2));
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xA12_2 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A12_2));
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xA21_2 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A21_2));
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xA22_2 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A22_2));

    // build map extractor
    Teuchos::RCP<Xpetra::EpetraMap> xfullmap_2 = Teuchos::rcp(new Xpetra::EpetraMap(fullmap));
    Teuchos::RCP<Xpetra::EpetraMap> xvelmap_2  = Teuchos::rcp(new Xpetra::EpetraMap(velmap ));
    Teuchos::RCP<Xpetra::EpetraMap> xpremap_2  = Teuchos::rcp(new Xpetra::EpetraMap(premap ));

    std::vector<Teuchos::RCP<const Xpetra::Map<LO,GO> > > xmaps_2;
    xmaps_2.push_back(xvelmap_2);
    xmaps_2.push_back(xpremap_2);

    Teuchos::RCP<const Xpetra::MapExtractor<Scalar,LO,GO> > map_extractor_2 = Xpetra::MapExtractorFactory<Scalar,LO,GO>::Build(xfullmap_2,xmaps_2);

    // build blocked operator
    Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO> > bOp_2 = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar,LO,GO>(map_extractor_2,map_extractor_2,10));
    bOp_2->setMatrix(0,0,xA11_2);
    bOp_2->setMatrix(0,1,xA12_2);
    bOp_2->setMatrix(1,0,xA21_2);
    bOp_2->setMatrix(1,1,xA22_2);

    bOp_2->fillComplete();
    //////////////////////////////

    RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO> > bOpbOp_2 = MueLu::Utils<Scalar,LO,GO>::TwoMatrixMultiplyBlock(bOp,false,bOp_2,false);

    //////////////////////////////
    // matrix-matrix multiplication of standard matrices
    RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xfuAop   = rcp(new Xpetra::CrsMatrix<Scalar,LO,GO>(xfuA));
    RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xfuAop_2 = rcp(new Xpetra::CrsMatrix<Scalar,LO,GO>(xfuA_2));
    Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO> > fuAfuA_2 = MueLu::Utils<Scalar,LO,GO>::TwoMatrixMultiply(xfuAop,false,xfuAop_2,false);

    /////////////////////////////
    // build vector
    Teuchos::RCP<Xpetra::EpetraVector> xx  = Teuchos::rcp(new Xpetra::EpetraVector(x));
    Teuchos::RCP<Xpetra::Vector<Scalar,LO,GO> > result =  Xpetra::VectorFactory<Scalar,LO,GO>::Build(xfullmap ,true);

    // matrix vector product
    bOpbOp_2->apply(*xx,*result,Teuchos::NO_TRANS);

    // check results
    Teuchos::RCP<Xpetra::Vector<Scalar,LO,GO> > result2 =  Xpetra::VectorFactory<Scalar,LO,GO>::Build(xfullmap ,true);
    fuAfuA_2->apply(*xx,*result2);

    result2->update(1.0,*result,-1.0);

    //cout << "norm of difference " << result2->norm2() << endl;

    TEUCHOS_TEST_COMPARE(result2->norm2(), <, 1e-16, out, success);

    //
    TEUCHOS_TEST_EQUALITY(bOpbOp_2->getGlobalNumRows(), fuAfuA_2->getGlobalNumRows(), out, success );
    TEUCHOS_TEST_EQUALITY(bOpbOp_2->getGlobalNumCols(), fuAfuA_2->getGlobalNumCols(), out, success );

    //RCP<Xpetra::EpetraCrsMatrix> rd = MueLu::Gallery::Random<Scalar, LO, GO, Xpetra::EpetraMap, Xpetra::EpetraCrsMatrix>(xvelmap,xpremap);

    //rd->describe(out,Teuchos::VERB_EXTREME);
#endif

  }

  /// simple test for matrix-matrix multiplication for a 2x2 blocked matrix with a 2x1 blocked matrix
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockedCrsMatrix, EpetraMatrixMatrixMult2x1, Scalar, LO, GO, Node)
  {
#ifdef HAVE_XPETRA_EPETRAEXT
    RCP<const Comm<int> > comm = getDefaultComm();

    // build maps
    RCP<Epetra_Map> rowmap1     = Teuchos::rcp(new Epetra_Map(24,0,*Xpetra::toEpetra(comm)));
    RCP<Epetra_Map> rowmap2     = Teuchos::rcp(new Epetra_Map(12,24,*Xpetra::toEpetra(comm)));
    RCP<Epetra_Map> dommap1     = Teuchos::rcp(new Epetra_Map(8,0,*Xpetra::toEpetra(comm)));
    RCP<Epetra_Map> dommap2     = Teuchos::rcp(new Epetra_Map(4,8,*Xpetra::toEpetra(comm)));

    std::vector<RCP<const Epetra_Map> > rowmaps;
    rowmaps.push_back(rowmap1); rowmaps.push_back(rowmap2);
    std::vector<RCP<const Epetra_Map> > dommaps;
    dommaps.push_back(dommap1); dommaps.push_back(dommap2);

    RCP<Epetra_Map> fullrowmap = MergeMaps(rowmaps);
    RCP<Epetra_Map> fulldommap = MergeMaps(dommaps);

    // read in matrices in matrix market format
    Epetra_CrsMatrix* ptrA   = 0;
    Epetra_CrsMatrix* ptrP   = 0;
    EpetraExt::MatrixMarketFileToCrsMatrix("A.mat",*fullrowmap,*fullrowmap,*fullrowmap,ptrA);
    EpetraExt::MatrixMarketFileToCrsMatrix("P.mat",*fullrowmap,*fullrowmap,*fulldommap,ptrP);
    Teuchos::RCP<Epetra_CrsMatrix> epA = Teuchos::rcp(ptrA);
    Teuchos::RCP<Epetra_CrsMatrix> epP = Teuchos::rcp(ptrP);

    // split fullA into A11,..., A22
    Teuchos::RCP<Epetra_CrsMatrix> epA11;
    Teuchos::RCP<Epetra_CrsMatrix> epA12;
    Teuchos::RCP<Epetra_CrsMatrix> epA21;
    Teuchos::RCP<Epetra_CrsMatrix> epA22;

    SplitMatrix2x2(epA,*rowmap1,*rowmap2,epA11,epA12,epA21,epA22);

    Teuchos::RCP<Epetra_CrsMatrix> epP1;
    Teuchos::RCP<Epetra_CrsMatrix> epP2;

    SplitMatrix2x1(epP,*rowmap1,*rowmap2,*fulldommap,epP1,epP2);

    ////////////////// transform Epetra stuff to Xpetra

    // build Xpetra objects from Epetra_CrsMatrix objects
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xA11 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epA11));
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xA12 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epA12));
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xA21 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epA21));
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xA22 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epA22));

    // build Xpetra objects from Epetra_CrsMatrix objects
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xP1 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epP1));
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xP2 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epP2));

    Teuchos::RCP<Xpetra::EpetraMap> xfullrowmap = Teuchos::rcp(new Xpetra::EpetraMap(fullrowmap));
    Teuchos::RCP<Xpetra::EpetraMap> xfulldommap = Teuchos::rcp(new Xpetra::EpetraMap(fulldommap));
    Teuchos::RCP<Xpetra::EpetraMap> xrowmap1  = Teuchos::rcp(new Xpetra::EpetraMap(rowmap1 ));
    Teuchos::RCP<Xpetra::EpetraMap> xrowmap2  = Teuchos::rcp(new Xpetra::EpetraMap(rowmap2 ));
    Teuchos::RCP<Xpetra::EpetraMap> xdommap1  = Teuchos::rcp(new Xpetra::EpetraMap(dommap1 ));
    Teuchos::RCP<Xpetra::EpetraMap> xdommap2  = Teuchos::rcp(new Xpetra::EpetraMap(dommap2 ));

    // build map extractor objects
    std::vector<Teuchos::RCP<const Xpetra::Map<LO,GO> > > xrowmaps;
    xrowmaps.push_back(xrowmap1);
    xrowmaps.push_back(xrowmap2);

    Teuchos::RCP<const Xpetra::MapExtractor<Scalar,LO,GO> > map_extractor = Xpetra::MapExtractorFactory<Scalar,LO,GO>::Build(xfullrowmap,xrowmaps);

    std::vector<Teuchos::RCP<const Xpetra::Map<LO,GO> > > xdommaps;
    xdommaps.push_back(xfulldommap);

    Teuchos::RCP<const Xpetra::MapExtractor<Scalar,LO,GO> > map_domextractor = Xpetra::MapExtractorFactory<Scalar,LO,GO>::Build(xfulldommap,xdommaps);

    // build blocked operators

    // build 2x2 blocked operator
    Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO> > bA = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar,LO,GO>(map_extractor,map_extractor,10));
    bA->setMatrix(0,0,xA11);
    bA->setMatrix(0,1,xA12);
    bA->setMatrix(1,0,xA21);
    bA->setMatrix(1,1,xA22);
    bA->fillComplete();

    // build 2x1 blocked operator
    Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO> > bP = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar,LO,GO>(map_extractor,map_domextractor,10));
    bP->setMatrix(0,0,xP1);
    bP->setMatrix(1,0,xP2);
    bP->fillComplete();

    RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO> > bAbP = MueLu::Utils<Scalar,LO,GO>::TwoMatrixMultiplyBlock(bA,false,bP,false);

    TEUCHOS_TEST_EQUALITY(bAbP->Rows(), 2, out, success );
    TEUCHOS_TEST_EQUALITY(bAbP->Cols(), 1, out, success );

    RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > bAbPmerged = bAbP->Merge();

    // tests
    Teuchos::RCP<Xpetra::Vector<Scalar,LO,GO> > onevector =  Xpetra::VectorFactory<Scalar,LO,GO>::Build(bAbP->getDomainMap() ,true);
    Teuchos::RCP<Xpetra::Vector<Scalar,LO,GO> > resvector =  Xpetra::VectorFactory<Scalar,LO,GO>::Build(bAbP->getRangeMap() ,true);
    Teuchos::RCP<Xpetra::Vector<Scalar,LO,GO> > resvector2=  Xpetra::VectorFactory<Scalar,LO,GO>::Build(bAbPmerged->getRangeMap() ,true);
    onevector->putScalar(1.0);
    bAbP->apply(*onevector,*resvector);
    bAbPmerged->apply(*onevector,*resvector2);

    resvector2->update(1.0,*resvector,-1.0);
    TEUCHOS_TEST_COMPARE(resvector2->norm2(), <, 1e-16, out, success);
#endif
  }

  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP_ORDINAL( SC, LO, GO, Node )                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockedCrsMatrix, EpetraApply, SC, LO, GO, Node ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockedCrsMatrix, EpetraMatrixMatrixMult, SC, LO, GO, Node ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockedCrsMatrix, EpetraMatrixMatrixMult2x1, SC, LO, GO, Node )

  typedef Kokkos::DefaultNode::DefaultNodeType DefaultNodeType;
  
  UNIT_TEST_GROUP_ORDINAL(double, int, int, DefaultNodeType)

}
