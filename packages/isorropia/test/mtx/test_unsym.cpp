//@HEADER
// ************************************************************************
//
//               Isorropia: Partitioning and Load Balancing Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
//
// ************************************************************************
//@HEADER

// Read in a matrix market file.  Use isorropia to do graph
// partitioning.  Compute the graph metrics both before and after partitioning.
//
// This tests Isorropia::Epetra::Partitioner followed by
// redistribution with a Isorropia::Epetra::Redistributor.
//
// For graph partitioning:
// The GRAPH_SYMMETRIZE option is assumed to be set to default "TRANSPOSE".
// This test has to be modified for other options.
//
// If run with --f={filename} a matrix market file other than west0067.mtx
// will be processed.
//

#include <Isorropia_ConfigDefs.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#include <ispatest_lbeval_utils.hpp>
#include <ispatest_epetra_utils.hpp>

#ifdef HAVE_EPETRA
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#ifdef HAVE_EPETRAEXT
#include <EpetraExt_CrsMatrixIn.h>
#endif
#endif

#include <Teuchos_CommandLineProcessor.hpp>

int main(int argc, char** argv)
{

  int rc=0, fail = 0;
#ifdef HAVE_EPETRAEXT
  int localProc = 0;

  int numProcs;
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  const Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  const Epetra_SerialComm Comm;
  numProcs = 1;
#endif

  Teuchos::CommandLineProcessor clp(false,true);

  // --f=fileName provides a different matrix market file for input
  // --run-all will continue to run all tests even if there is a failure

  std::string *inputFile = new std::string("west0067.mtx");
  bool runAll = false; 

  clp.setOption( "f", inputFile,
		"Name of input matrix market file");
  clp.setOption( "run-all", "abort", &runAll,
		"Don't abort if one test fails, run all of them.");

  Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return =
    clp.parse(argc,argv);

  if( parse_return == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
  {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
  }
  if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL )
  {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 1;
  }

  const char *fname = inputFile->c_str();

  // Read in the matrix market file and distribute its rows across the
  // processes.
  //
  // This reader uses the default Epetra_Map for number of rows for the
  // RowMap() and for the RangeMap().  For non-square matrices it uses
  // the default Epetra_Map for the number of columns for the DomainMap(),
  // otherwise it uses the RowMap().
  //
  // The maps can be specified with other versions of MMFtoCrsMatrix().

  Epetra_CrsMatrix *matrixPtr;
  rc = EpetraExt::MatrixMarketFileToCrsMatrix(fname, Comm, matrixPtr);
  if (rc < 0)
  {
    if (localProc==0) std::cerr << "error reading input file"<< std::cout;
    return(1);
  }

  bool square = (matrixPtr->NumGlobalRows() == matrixPtr->NumGlobalCols());

  // Run some partitioning tests
  //   Test graph partitioning with square unsymmetric matrices

  Teuchos::RCP<Epetra_CrsMatrix> testm = Teuchos::rcp(matrixPtr);
  Teuchos::RCP<const Epetra_RowMatrix> rm = testm;
  Isorropia::Epetra::CostDescriber costs;

  double myShareBefore = 1.0 / numProcs;
  double myShare = myShareBefore;

  double balance1, balance2, cutn1, cutn2, cutl1, cutl2;
  double cutWgt1, cutWgt2;
  int numCuts1, numCuts2;
  rc = ispatest::compute_graph_metrics(*rm, costs, myShare, 
                balance1, numCuts1, cutWgt1, cutn1, cutl1);

  if (square)
  {
    Teuchos::ParameterList params;
    params.set("partitioning method", "graph");
    Isorropia::Epetra::Partitioner iso_part (rm, params);

    Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner =
                                                Teuchos::rcpFromRef (iso_part);
    Isorropia::Epetra::Redistributor redist(partitioner);

    Teuchos::RCP<Epetra_CrsMatrix> bal_matrix;
    try
    {
        bal_matrix = redist.redistribute(*testm);
    }
    catch (std::exception& exc)
    {
        std::cout << "Isorropia exception " << exc.what() << std::endl;
        return  1;
    }
    fail = 0; // Bad way to test, If Zoltan did not crash then success !!
              // But all metrics could get worse from the initial metrics
              // Leaving metrics testing to zdrive
    rc = ispatest::compute_graph_metrics(*bal_matrix, costs, myShare,
                balance2, numCuts2, cutWgt2, cutn2, cutl2);
    /*if (balance1 >= balance2)
    {
        fail = 0;
    }
    else 
    {
        fail = 1;
    }*/
    if (localProc == 0)
    {
        cout << "CutN before and after :" << cutn1 << " " << cutn2 << endl;
        cout << "CutL before and after :" << cutl1 << " " << cutl2 << endl;
        cout << "Balance before and after :" << balance1 << " " << balance2
                << endl;
        cout << "NumCuts before and after :" << numCuts1 << " " << numCuts2
                << endl;
    }
  }
  else
  {
      // We do not support this case now , TODO : Add it when we support
      // GRAPH_SYMMETRIZE as bipartite graph
      if (localProc == 0)
      {
        std::cout << "Test not run because i/p is not square" << std::endl;
      }
      fail = 0; 
  }


#else
  fail = 0;
  if (localProc == 0)
  {
    std::cout << "Test not run because it requires EPETRA_EXT" << std::endl;
  }
#endif

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return fail;
}
