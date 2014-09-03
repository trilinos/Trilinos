/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

//
// INTRO: PLEASE READ!!!
//
// This example will show how to build and use a Vanka smoother for
// matrices of the type
//
//   A = [  F   B_1^T ]
//       [ B_2    -C  ]
//
// For the example shown below, we take B_1 = B_2 and C = 0. This
// might correspond to, say, (Navier-)Stokes or some such
// thing. That's not important right now.
//
// To do Vanka, we want to think of grabbing a row of B. Then we look
// at the nonzero columns of B and we grab the corresponding
// rows/columns of F. This is our "Vanka block" for this row of
// B. There is one Vanka block for every row of B. The Vanka method is
// then doing block Gauss-Seidel with the Vanka blocks.
//
//
//
// IMPORTANT: If your matrix comes from an FEM code, for example, you
// may run into troubles. If you happen to be using stable elements,
// such as Taylor-Hood elements for (Navier-)Stokes, your
// discretization software should correctly compute that C =
// 0. However, the graph of this matrix will not be empty! It will
// store a bunch of zeros corresponding to the element structure. This
// means that when you blindly grab a row of B, it will have
// connections that correspond to columns of F (which you want) and
// connections to columns of C (which you don't (necessarily)
// want). Again to use (Navier-)Stokes terminology, your pressure will
// be connected to a bunch of velocities (which you want) and a bunch
// of other pressures (which you don't (necessarily) want).
//
// At this time, and probably fairly permanently,
// Ifpack2::Details::UserPartitioner does nothing to deal with
// this. (How can it? It just finds graph entries, and these graph
// entries are there...) Thus if you are encountering this phenomenon
// and are displeased by it, your best options are either to eliminate
// the zero-graph completely or to write a custom partitioner that
// computes the same partition that Ifpack2::Details::UserPartitioner
// computes and remove the rows corresponding to unwanted connections
// after the fact. Also as an alternative, you could write a custom
// partition that really has nothing to do with the graph of the
// matrix that does whatever you want, however you want to do it. But
// if you do that by following the Ifpack2::Partitioner interface, you
// probably have to lie and at least pretend to care about the
// graph. But that just seems silly! How could we ever hope to
// partition rows of a matrix in any other way but somehow based on
// some graph information?!
//
//
//
// Anyway, that's my story and I'm sticking to it. Good luck.
//


// Teuchos
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#else
#include "Teuchos_DefaultSerialComm.hpp"
#endif

// Tpetra
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Map.hpp"

// Ifpack2
#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_AdditiveSchwarz_def.hpp"
#include "Ifpack2_BlockRelaxation_decl.hpp"
#include "Ifpack2_BlockRelaxation_def.hpp"
#include "Ifpack2_DenseContainer_decl.hpp"
#include "Ifpack2_DenseContainer_def.hpp"
#include "Ifpack2_Details_UserPartitioner.hpp"
#include "Ifpack2_OverlappingPartitioner.hpp"
#include "Ifpack2_Parameters.hpp"
#include "Ifpack2_Partitioner.hpp"
#include "Ifpack2_Preconditioner.hpp"

// C++ stuff
#include <iostream>


int
main (int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpisess (&argc, &argv);

  //Define communicator:
  Tpetra::DefaultPlatform::DefaultPlatformType& platform =
    Tpetra::DefaultPlatform::getDefaultPlatform ();
  Teuchos::RCP<const Teuchos::Comm<int> > comm = platform.getComm ();

  // Fancy ostream, just in case
  Teuchos::RCP<Teuchos::FancyOStream> out =
    Teuchos::fancyOStream (Teuchos::rcpFromRef (std::cout));
  out->setOutputToRootOnly (0);

  // General stuff
  typedef Tpetra::Vector<>::scalar_type ST;
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<LO>::global_ordinal_type GO;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;

  // Matrix stuff
  typedef Tpetra::Map<LO,GO,Node> MAP;
  typedef Tpetra::CrsMatrix<ST,LO,GO,Node> CRS;
  typedef Tpetra::Vector<ST,LO,GO,Node> VEC;

  // Ifpack2 stuff
  typedef Ifpack2::Preconditioner< ST, LO, GO, Node > PrecType;
  typedef Ifpack2::DenseContainer< CRS, ST > ContainerType;
  typedef Ifpack2::BlockRelaxation< CRS, ContainerType > BlockRelax;
  typedef Ifpack2::AdditiveSchwarz< CRS > TheSchwarz;

  // using stuff
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::rcp;
  using Teuchos::RCP;


  // Define the problem size. Let's define a local problem for
  // simplicity.
  size_t numNodeRowsF = 6;
  size_t numNodeRowsB = numNodeRowsF / 2;
  size_t numNodeRowsTotal = numNodeRowsF + numNodeRowsB;


  //
  // This will create a matrix of the form
  //
  //   A = [ F  B^T ]
  //       [ B   0  ].
  //
  // Locally, A will have numNodeRowsF, B will have numNodeRowsB, etc.
  //


  // Create a map.
  //
  // The first numNodeRowsF local elements in the map will be the
  // local rows of F. The remaining local elements will be local rows
  // of B.
  RCP< const MAP > myMap = rcp( new MAP( Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                                         numNodeRowsTotal,
                                         0,
                                         comm ) );


  // Start building a CrsMatrix using this Map
  ArrayView<const GO> myGIDs = myMap->getNodeElementList();

  // Create a CrsMatrix using the map
  RCP< CRS > A = rcp(new CRS(myMap,4));

  LO bLID;

  // Add rows to "[F B^T]" first
  for (size_t row = 0; row < numNodeRowsF; ++row){
    if (row == 0) { //left boundary

      bLID = numNodeRowsF;

      A->insertGlobalValues(myGIDs[row],
                            Teuchos::tuple<GO>(myGIDs[row],myGIDs[row+1],myGIDs[bLID]),
                            Teuchos::tuple<ST>(2.0,-1.0,1.0));

    }
    else if ((size_t) row == numNodeRowsF - 1) { //right boundary

      bLID = numNodeRowsTotal - 1;

      A->insertGlobalValues(myGIDs[row],
                            Teuchos::tuple<GO>(myGIDs[row-1],myGIDs[row],myGIDs[bLID]),
                            Teuchos::tuple<ST>(-1.0,2.0,1.0));
    }
    else { //interior

      bLID = numNodeRowsF + ( row % numNodeRowsB );

      A->insertGlobalValues(myGIDs[row],
                            Teuchos::tuple<GO>(myGIDs[row-1],myGIDs[row],myGIDs[row+1],myGIDs[bLID]),
                            Teuchos::tuple<ST>(-1.0,2.0,-1.0,1.0));
    }
  }

  // Now add rows to "[B 0]"
  for (size_t row = numNodeRowsF; row < numNodeRowsTotal; ++row){
    A->insertGlobalValues(myGIDs[row],
                          Teuchos::tuple<GO>(myGIDs[row-numNodeRowsF],myGIDs[row-numNodeRowsF+numNodeRowsB]),
                          Teuchos::tuple<ST>(1.0,1.0));
  }

  // Complete the fill
  A->fillComplete();


  //
  // At this point, the matrix should be ready to go. F is
  // tridiagonal, B has two entries per row. Thus, our Vanka blocks
  // will be 3x3 matrices -- a 2x2 F_\ell block and a 1x2 B_\ell
  // block.
  //


  // Set up ifpack2 parameters
  Teuchos::ParameterList MyList;
  // We are defining a custom nonoverlapping partition (see below).
  MyList.set("partitioner: type","user");
  // There is one block for each local row of B.
  MyList.set("partitioner: local parts",(int) numNodeRowsB);
  // We want to grab the rows/columns of F that correspond to the
  // nonzero entries in a row of B, so this is a lenght-1 graph
  // connection.
  MyList.set("partitioner: overlap"    ,(int) 1);
  // "Vanka" smoothing is block Gauss-Seidel
  MyList.set("relaxation: type"        ,"Gauss-Seidel");
  MyList.set("relaxation: sweeps"      ,(int) 1);
  MyList.set("relaxation: damping factor",1.0);
  MyList.set("relaxation: zero starting solution", false);
  // Use the Schwarz
  MyList.set("schwarz: overlap level", (int) 0);


  //
  // The last thing we need to do here is to set the nonoverlapping
  // partition on this proc.
  //
  // Every row has to be assigned to a nonoverlapping partition. If a
  // row is not a "seed row" for Vanka, it is given to nonoverlapping
  // partition invalid(). Rows that have valid nonoverlapping
  // partition numbers then become the "seed rows" for Vanka.
  //
  // In this case, we want the seed rows to be the rows of B, and we
  // want only one seed row per Vanka block.
  //
  ArrayRCP<LO> blockSeeds( numNodeRowsTotal,
                           Teuchos::OrdinalTraits<LO>::invalid() );

  for (size_t rowOfB = numNodeRowsF; rowOfB < numNodeRowsTotal; ++rowOfB ){
    blockSeeds[rowOfB] = rowOfB - numNodeRowsF; // This starts the partition indexing
                                                // at 0. I don't know if that is required
                                                // or not...
  }


  MyList.set("partitioner: map", blockSeeds);


  // Create the preconditioner
  RCP<TheSchwarz> prec = rcp(new TheSchwarz(A));

  // Set the parameters
  prec->setParameters(MyList);

  //
  // Now things get a little weird. With the way preconditioner
  // classes now work, we need to create our inner solver separately.
  //

  RCP<PrecType> innerPrec = rcp(new BlockRelax(A));

  // Set the BlockRelaxation parameters
  innerPrec->setParameters(MyList);

  // Set BlockRelaxation as the inner preconditioner for Additive Schwarz
  prec->setInnerPreconditioner(innerPrec);


  //
  // Now that our AdditiveSchwarz preconditioner has everything it
  // needs, it can be initialized and computed.
  //


  prec->initialize();
  prec->compute();


  // Define RHS / Initial Guess
  VEC X( myMap , true );
  VEC B( myMap , false );
  B.putScalar((ST) 2.0);

  //
  // This sets X = 0; B = 2;
  //


  // Apply the preconditioner
  prec->apply(B,X);

  Array<ST> Answer(numNodeRowsTotal);
  Answer[0] = 1.0;
  Answer[1] = 1.0;
  Answer[2] = 5.0 / 4;
  Answer[3] = 1.0;
  Answer[4] = 1.0;
  Answer[5] = 3.0 / 4;
  Answer[6] = 0.0;
  Answer[7] = 1.0;
  Answer[8] = 3.0 / 2;

  // Get data out of X
  Teuchos::ArrayRCP<const ST> ComputedSol = X.getData();
  // I know the ArrayRCP is "bad" here, but there are only 9
  // entries... I think we'll all be ok.

  // Compare data
  for (size_t ii = 0; ii < numNodeRowsTotal; ++ii){
    if ( Answer[ii] - ComputedSol[ii] > 1e-14 ){
      *out << "The test failed... Sad face.\n";
      return -1;
    }
  }


  // Success!
  *out << "The test succeeded!\n";
  return 0;


}//int main
