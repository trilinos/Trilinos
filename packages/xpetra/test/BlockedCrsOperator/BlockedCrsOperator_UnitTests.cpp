/*
 * BlockedCrsOperator_UnitTests.cpp
 *
 *  Created on: Aug 22, 2011
 *      Author: wiesner
 */


#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include "Epetra_MpiComm.h"
#include "mpi.h"
#include "Epetra_SerialComm.h"

// EpetraExt
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_RowMatrixOut.h"

#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <Teuchos_as.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_Operator.hpp>
#include <Xpetra_CrsOperator.hpp>
#ifdef HAVE_XPETRA_TPETRA
#include <Xpetra_TpetraCrsMatrix.hpp>
#endif
#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraCrsMatrix.hpp>
#endif
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_BlockedCrsOperator.hpp>

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
  using Xpetra::Operator;
  using Xpetra::CrsOperator;
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


  bool SplitMatrix2x2(Teuchos::RCP<const Epetra_CrsMatrix> A,
                      Teuchos::RCP<const Epetra_Map>& A11rowmap,
                      Teuchos::RCP<const Epetra_Map>& A22rowmap,
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

    if (A11rowmap==Teuchos::null && A22rowmap != Teuchos::null)
      A11rowmap = SplitMap(A->RowMap(),*A22rowmap);
    else if (A11rowmap != Teuchos::null && A22rowmap != Teuchos::null);
    else if (A11rowmap != Teuchos::null && A22rowmap == Teuchos::null)
      A22rowmap = SplitMap(A->RowMap(),*A11rowmap);
    else
    {
      cout << "SplitMatrix2x2: Both A11rowmap and A22rowmap == null on entry" << endl;
      return false;
    }

    const Epetra_Comm& Comm   = A->Comm();
    const Epetra_Map&  A22map = *(A22rowmap.get());
    const Epetra_Map&  A11map = *(A11rowmap.get());

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
  //////////////////////////////////////////////////////////////////////////

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignord and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  /*RCP<const Comm<int> > getDefaultComm()
  {
    if (testMpi) {
      return DefaultPlatform::getDefaultPlatform().getComm();
    }
    return rcp(new Teuchos::SerialComm<int>());
  }*/

  //
  // UNIT TESTS
  //


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockedCrsOperator, EpetraApply, Scalar, LO, GO, Node ) //TODO: add template parameter <Node,...>
  {
#ifdef HAVE_XPETRA_EPETRA

    RCP<Epetra_Comm> Comm;
    if(testMpi)
      Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
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

    TEST_EQUALITY(SplitMatrix2x2(fullA,velmap,premap,A11,A12,A21,A22),true);

    // build Xpetra objects from Epetra_CrsMatrix objects
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xfuA = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(fullA));
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
    Teuchos::RCP<Xpetra::BlockedCrsOperator<Scalar,LO,GO> > bOp = Teuchos::rcp(new Xpetra::BlockedCrsOperator<Scalar,LO,GO>(map_extractor,map_extractor,10));
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
    Teuchos::RCP<Xpetra::Vector<Scalar,LO,GO> > result2 =  Xpetra::VectorFactory<Scalar,LO,GO>::Build(xfullmap ,true);
    xfuA->apply(*xx,*result2);

    result2->update(1.0,*result,-1.0);

    //cout << "norm of difference " << result2->norm2() << endl;

    TEUCHOS_TEST_COMPARE(result2->norm2(), <, 1e-16, out, success);


#endif

  }

  //
  // INSTANTIATIONS
  //

#   define UNIT_TEST_GROUP_ORDINAL( SC, LO, GO, Node )                       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockedCrsOperator, EpetraApply, SC, LO, GO, Node )

  typedef Kokkos::DefaultNode::DefaultNodeType DefaultNodeType;
  UNIT_TEST_GROUP_ORDINAL(double, int, int, DefaultNodeType)

}

