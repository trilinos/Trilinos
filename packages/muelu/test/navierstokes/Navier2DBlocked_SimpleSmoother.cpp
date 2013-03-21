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

#include <unistd.h>
#include <iostream>
#include <fstream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ArrayRCP.hpp>

// Epetra
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_VectorIn.h>
#include <EpetraExt_MultiVectorIn.h>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_StridedEpetraMap.hpp>

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Memory.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_CoupledAggregationFactory.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_PreDropFunctionConstVal.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_EpetraOperator.hpp"
#include "MueLu_SubBlockAFactory.hpp"
#include "MueLu_SimpleSmoother.hpp"
#include "MueLu_SchurComplementFactory.hpp"
#include "MueLu_Utilities.hpp"

//TODO is it really needed?
#include "MueLu_HierarchyHelpers.hpp"


#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

#include <Epetra_LinearProblem.h>
#include <AztecOO.h>

// helper routines
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

/*!
 *  2d Navier Stokes example (for Epetra)
 *
 *  using block matrices
 */


int main(int argc, char *argv[]) {
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  //
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);
  *out << MueLu::MemUtils::PrintMemoryUsage() << std::endl;

  // Timing
  Teuchos::Time myTime("global");
  Teuchos::TimeMonitor MM(myTime);

#ifndef HAVE_TEUCHOS_LONG_LONG_INT
  *out << "Warning: scaling test was not compiled with long long int support" << std::endl;
#endif

  // read in input parameters

  // default parameters
  LO SIMPLE_nSweeps = 100;
  Scalar SIMPLE_omega = 0.02;
  LO SC_nSweeps = 1;
  Scalar SC_omega = 1.0;
  LO PRED_nSweeps = 3;
  Scalar PRED_omega = 1.0;

  int SC_bUseDirectSolver = 1;

  // Note: use --help to list available options.
  Teuchos::CommandLineProcessor clp(false);
  clp.setOption("SIMPLE_sweeps",&SIMPLE_nSweeps,"number of sweeps with SIMPLE smoother");
  clp.setOption("SIMPLE_omega", &SIMPLE_omega,  "scaling factor for SIMPLE smoother");
  clp.setOption("Predict_sweeps", &PRED_nSweeps,  "number of sweeps for SIMPLE internal velocity prediction smoother (GaussSeidel)");
  clp.setOption("Predict_omega", &PRED_omega,  "damping parameter for SIMPLE internal velocity prediction smoother (GaussSeidel)");
  clp.setOption("SchurComp_sweeps",    &SC_nSweeps,"number of sweeps for SIMPLE internal SchurComp solver/smoother (GaussSeidel)");
  clp.setOption("SchurComp_omega",     &SC_omega,  "damping parameter for SIMPLE internal SchurComp solver/smoother (GaussSeidel)");
  clp.setOption("SchurComp_solver",    &SC_bUseDirectSolver,  "if 1: use direct solver for SchurComp equation, otherwise use GaussSeidel smoother (=default)");

  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_ERROR:
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  int globalNumDofs = 8898;  // used for the maps
  //int nDofsPerNode = 3;      // used for generating the fine level null-space

  // build strided maps
  // striding information: 2 velocity dofs and 1 pressure dof = 3 dofs per node
  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(2);
  stridingInfo.push_back(1);

  /////////////////////////////////////// build strided maps
  // build strided maps:
  // xstridedfullmap: full map (velocity and pressure dof gids), continous
  // xstridedvelmap: only velocity dof gid maps (i.e. 0,1,3,4,6,7...)
  // xstridedpremap: only pressure dof gid maps (i.e. 2,5,8,...)
  Teuchos::RCP<Xpetra::StridedEpetraMap> xstridedfullmap = Teuchos::rcp_dynamic_cast<Xpetra::StridedEpetraMap>(StridedMapFactory::Build(Xpetra::UseEpetra,globalNumDofs,0,stridingInfo,comm,-1));
  Teuchos::RCP<Xpetra::StridedEpetraMap> xstridedvelmap  = Teuchos::rcp_dynamic_cast<Xpetra::StridedEpetraMap>(StridedMapFactory::Build(Xpetra::UseEpetra,globalNumDofs,0,stridingInfo,comm,0));
  Teuchos::RCP<Xpetra::StridedEpetraMap> xstridedpremap  = Teuchos::rcp_dynamic_cast<Xpetra::StridedEpetraMap>(StridedMapFactory::Build(Xpetra::UseEpetra,globalNumDofs,0,stridingInfo,comm,1));

  /////////////////////////////////////// transform Xpetra::Map objects to Epetra
  // this is needed for our splitting routine
  const Teuchos::RCP<const Epetra_Map> fullmap = Teuchos::rcpFromRef(xstridedfullmap->getEpetra_Map());
  Teuchos::RCP<const Epetra_Map> velmap = Teuchos::rcpFromRef(xstridedvelmap->getEpetra_Map());
  Teuchos::RCP<const Epetra_Map> premap = Teuchos::rcpFromRef(xstridedpremap->getEpetra_Map());

  /////////////////////////////////////// import problem matrix and RHS from files (-> Epetra)

  // read in problem
  Epetra_CrsMatrix * ptrA = 0;
  Epetra_Vector * ptrf = 0;
  Epetra_MultiVector* ptrNS = 0;

  *out << "Reading matrix market file" << std::endl;

  EpetraExt::MatrixMarketFileToCrsMatrix("A5932_re1000.txt",*fullmap,*fullmap,*fullmap,ptrA);
  EpetraExt::MatrixMarketFileToVector("b5932_re1000.txt",*fullmap,ptrf);

  RCP<Epetra_CrsMatrix> epA = Teuchos::rcp(ptrA);
  RCP<Epetra_Vector> epv = Teuchos::rcp(ptrf);
  RCP<Epetra_MultiVector> epNS = Teuchos::rcp(ptrNS);


  /////////////////////////////////////// split system into 2x2 block system

  *out << "Split matrix into 2x2 block matrix" << std::endl;

  // split fullA into A11,..., A22
  Teuchos::RCP<Epetra_CrsMatrix> A11;
  Teuchos::RCP<Epetra_CrsMatrix> A12;
  Teuchos::RCP<Epetra_CrsMatrix> A21;
  Teuchos::RCP<Epetra_CrsMatrix> A22;

  if(SplitMatrix2x2(epA,*velmap,*premap,A11,A12,A21,A22)==false)
    *out << "Problem with splitting matrix"<< std::endl;

  /////////////////////////////////////// transform Epetra objects to Xpetra (needed for MueLu)

  // build Xpetra objects from Epetra_CrsMatrix objects
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xA11 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A11));
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xA12 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A12));
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xA21 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A21));
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xA22 = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A22));

  /////////////////////////////////////// generate MapExtractor object

  std::vector<Teuchos::RCP<const Xpetra::Map<LO,GO,Node> > > xmaps;

  xmaps.push_back(xstridedvelmap);
  xmaps.push_back(xstridedpremap);

  Teuchos::RCP<const Xpetra::MapExtractor<Scalar,LO,GO,Node> > map_extractor = Xpetra::MapExtractorFactory<Scalar,LO,GO>::Build(xstridedfullmap,xmaps);

  /////////////////////////////////////// build blocked transfer operator
  // using the map extractor
  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> > bOp = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar,LO,GO>(map_extractor,map_extractor,10));
  bOp->setMatrix(0,0,xA11);
  bOp->setMatrix(0,1,xA12);
  bOp->setMatrix(1,0,xA21);
  bOp->setMatrix(1,1,xA22);

  bOp->fillComplete();
  //////////////////////////////////////////////////////// finest Level
  RCP<MueLu::Level> Finest = rcp(new Level());
  Finest->setDefaultVerbLevel(Teuchos::VERB_NONE);
  Finest->Set("A",Teuchos::rcp_dynamic_cast<Matrix>(bOp));


  ///////////////////////////////////
  // Test Braess Sarazin Smoother as a solver

  *out << "Test: Creating SIMPLE Smoother" << std::endl;
  *out << "Test: Omega for SIMPLE = " << SIMPLE_omega << std::endl;
  *out << "Test: Number of sweeps for SIMPLE = " << SIMPLE_nSweeps << std::endl;
  *out << "Test: Omega for Schur Complement solver= " << SC_omega << std::endl;
  *out << "Test: Number of Schur Complement solver= " << SC_nSweeps << std::endl;
  *out << "Test: Setting up Braess Sarazin Smoother" << std::endl;

  // define SIMPLE Smoother with SIMPLE_nSweeps and SIMPLE_omega as scaling factor
  // AFact_ = Teuchos::null (= default) for the 2x2 blocked operator
  RCP<SimpleSmoother> SimpleSm = rcp( new SimpleSmoother(SIMPLE_nSweeps,SIMPLE_omega) );

  RCP<SmootherFactory>   smootherFact          = rcp( new SmootherFactory(SimpleSm) );

  // define smoother for velocity prediction
  RCP<SubBlockAFactory> A00Fact = Teuchos::rcp(new SubBlockAFactory(MueLu::NoFactory::getRCP(), 0, 0));
  RCP<SmootherPrototype> smoProtoPredict = Teuchos::null;
  std::string ifpackPredictType;
  Teuchos::ParameterList ifpackPredictList;
  ifpackPredictList.set("relaxation: sweeps", PRED_nSweeps );
  ifpackPredictList.set("relaxation: damping factor", PRED_omega );
  ifpackPredictType = "RELAXATION";
  ifpackPredictList.set("relaxation: type", "Gauss-Seidel");
  smoProtoPredict = rcp( new TrilinosSmoother(ifpackPredictType, ifpackPredictList, 0, A00Fact) );
  RCP<SmootherFactory> SmooPredictFact = rcp( new SmootherFactory(smoProtoPredict) );
  // define temporary FactoryManager that is used as input for BraessSarazin smoother
  RCP<FactoryManager> MPredict = rcp(new FactoryManager());
  MPredict->SetFactory("A",                 A00Fact);         // SchurComplement operator for correction step (defined as "A")
  MPredict->SetFactory("Smoother",          SmooPredictFact);    // solver/smoother for correction step
  MPredict->SetFactory("PreSmoother",               SmooPredictFact);
  MPredict->SetFactory("PostSmoother",              SmooPredictFact);
  MPredict->SetIgnoreUserData(true);               // always use data from factories defined in factory manager
  SimpleSm->SetVelocityPredictionFactoryManager(MPredict);    // set temporary factory manager in BraessSarazin smoother


  // define SchurComplement Factory
  // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
  // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
  // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
  RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
  SFact->SetParameter("omega", Teuchos::ParameterEntry(1.0)); // for Simple, omega is always 1.0 in the SchurComplement
  SFact->SetFactory("A",MueLu::NoFactory::getRCP());

  // define smoother/solver for BraessSarazin
  RCP<SmootherPrototype> smoProtoSC = Teuchos::null;
  if(SC_bUseDirectSolver != 1) {
    //Smoother Factory, using SFact as a factory for A
    std::string ifpackSCType;
    Teuchos::ParameterList ifpackSCList;
    ifpackSCList.set("relaxation: sweeps", SC_nSweeps );
    ifpackSCList.set("relaxation: damping factor", SC_omega );
    ifpackSCType = "RELAXATION";
    ifpackSCList.set("relaxation: type", "Gauss-Seidel");
    smoProtoSC     = rcp( new TrilinosSmoother(ifpackSCType, ifpackSCList, 0, SFact) );
  }
  else {
    Teuchos::ParameterList ifpackDSList;
    std::string ifpackDSType;
    smoProtoSC     = rcp( new DirectSolver(ifpackDSType,ifpackDSList) ); smoProtoSC->SetFactory("A", SFact);
  }

  RCP<SmootherFactory> SmooSCFact = rcp( new SmootherFactory(smoProtoSC) );

  // define temporary FactoryManager that is used as input for BraessSarazin smoother
  RCP<FactoryManager> MB = rcp(new FactoryManager());
  MB->SetFactory("A",                 SFact);         // SchurComplement operator for correction step (defined as "A")
  MB->SetFactory("Smoother",          SmooSCFact);    // solver/smoother for correction step
  MB->SetFactory("PreSmoother",               SmooSCFact);
  MB->SetFactory("PostSmoother",              SmooSCFact);
  MB->SetIgnoreUserData(true);               // always use data from factories defined in factory manager
  SimpleSm->SetSchurCompFactoryManager(MB);    // set temporary factory manager in BraessSarazin smoother

  // setup main factory manager
  RCP<FactoryManager> M = rcp(new FactoryManager());
  M->SetFactory("A",               MueLu::NoFactory::getRCP()); // this is the 2x2 blocked operator
  M->SetFactory("Smoother",        smootherFact);               // BraessSarazin block smoother
  M->SetFactory("PreSmoother",     smootherFact);
  M->SetFactory("PostSmoother",    smootherFact);

  MueLu::SetFactoryManager SFMCoarse(Finest, M);
  Finest->Request(MueLu::TopSmootherFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>(M, "Smoother"));

  // call setup (= extract blocks and extract diagonal of F)
  SimpleSm->Setup(*Finest);

  RCP<MultiVector> xtest = MultiVectorFactory::Build(xstridedfullmap,1);
  xtest->putScalar( (SC) 0.0);

  RCP<Vector> xR = Teuchos::rcp(new Xpetra::EpetraVector(epv));
  // calculate initial (absolute) residual
  Teuchos::Array<ST::magnitudeType> norms(1);

  xR->norm2(norms);
  *out << "Test: ||x_0|| = " << norms[0] << std::endl;
  *out << "Test: Applying Simple Smoother" << std::endl;
  *out << "Test: START DATA" << std::endl;
  *out << "iterations\tVelocity_residual\tPressure_residual" << std::endl;
  SimpleSm->Apply(*xtest,*xR);
  xtest->norm2(norms);
  *out << "Test: ||x_1|| = " << norms[0] << std::endl;

  Teuchos::Array<Teuchos::ScalarTraits<double>::magnitudeType> test = MueLu::Utils<double, int, int>::ResidualNorm(*bOp, *xtest, *xR);
  *out << "residual norm: " << test[0] << std::endl;

  return EXIT_SUCCESS;
}



