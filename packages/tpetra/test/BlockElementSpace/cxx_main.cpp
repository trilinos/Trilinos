// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#define ORDINALTYPE int

#include <iostream>
#include <iomanip>
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_BlockElementSpace.hpp"
#include "Tpetra_Version.hpp"

#ifdef TPETRA_MPI
#include <mpi.h>
#include "Tpetra_MpiPlatform.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#endif // TPETRA_MPI

int main(int argc, char* argv[]) {

	const ORDINALTYPE INDEXBASE = 0;
	const ORDINALTYPE NUMELEMENTS = 5;
	const ORDINALTYPE ELEMENTSIZE = 2;
  
	bool verbose = false;
	bool debug = false;
	if(argc > 1) {
		if(argv[1][0] == '-' && argv[1][1] == 'v')
			verbose = true;
		if(argv[1][0] == '-' && argv[1][1] == 'd') {
			debug = true;
			verbose = true;
		}
	}

  int rank = 0; // assume we are on serial
  int size = 1; // if MPI, will be reset later
  
  // initialize MPI if needed
#ifdef TPETRA_MPI
  size = -1;
  rank = -1;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(verbose) cout << "MPI Startup: Image " << rank << " of " << size << " is alive." << endl;
  MPI_Barrier(MPI_COMM_WORLD);
#endif // TPETRA_MPI
  
  // change verbose to only be true on Image 0
  verbose = (verbose && (rank == 0));
  
  // start the testing
	if(verbose) {
    cout << "\n****************************************\n" 
         << "Starting BlockElementSpaceTest..." << endl
         << Tpetra::Tpetra_Version() << endl
         << "****************************************\n";
  }
  
  // Platform
  if(verbose) cout << "Creating platform" << endl;
#ifdef TPETRA_MPI
  Tpetra::MpiPlatform<ORDINALTYPE, ORDINALTYPE> platform(MPI_COMM_WORLD);
#else
  Tpetra::SerialPlatform<ORDINALTYPE, ORDINALTYPE> platform;
#endif // TPETRA_MPI


  // ElementSpace
	// commented out lines are for creating alternate es objects.
  if(verbose) cout << "Creating es, constructor1" << endl;
  Tpetra::ElementSpace<ORDINALTYPE> es(NUMELEMENTS, INDEXBASE, platform);
  //if(verbose) cout << "Creating es, constructor2" << endl;
  //Tpetra::ElementSpace<ORDINALTYPE> es(-1, NUMELEMENTS, INDEXBASE, platform);
  //if(verbose) cout << "Creating es, constructor3" << endl;
  //ORDINALTYPE gidList[NUMELEMENTS] = {1,4,7,8,9};//,15,22,54,55,58};
  //Tpetra::ElementSpace<ORDINALTYPE> es(-1, NUMELEMENTS, gidList, INDEXBASE, platform);
  //if(debug) cout << es;

  //BlockElementSpace
  if(verbose) cout << "Creating bes1(fixed-size)...";
  Tpetra::BlockElementSpace<ORDINALTYPE> bes1(es, ELEMENTSIZE);
  if(verbose) cout << "Successful." << endl;
  if(debug) cout << bes1;

	if(verbose) cout << "Creating bes2(variable-sized)...";
	ORDINALTYPE* esizelist = new ORDINALTYPE[NUMELEMENTS];
	esizelist[0] = 1;
	esizelist[1] = 1;
	esizelist[2] = 2;
	esizelist[3] = 3;
	esizelist[4] = 5;
	Tpetra::BlockElementSpace<ORDINALTYPE> bes2(es, esizelist);
	delete[] esizelist;
	esizelist = 0;
  if(verbose) cout << "Successful." << endl;
	if(debug) cout << bes2;

  if(verbose) cout <<"Creating bes3(copy constructor)...";
  Tpetra::BlockElementSpace<ORDINALTYPE> bes3(bes1);
  if(verbose) cout << "Successful." << endl;
  if(debug) cout << bes3;

  if(verbose) cout << "Checking isSameAs...";
  Tpetra::BlockElementSpace<ORDINALTYPE> bes4(es, ELEMENTSIZE+1);
  assert(bes1.isSameAs(bes3) == true);
  assert(bes1.isSameAs(bes4) == false);
  if(verbose) cout << "Successful." << endl;

  if(verbose) cout << "Checking assignment operator...";
  assert(bes1.isSameAs(bes4) == false);
  bes4 = bes1;
  assert(bes1.isSameAs(bes4) == true);
  if(verbose) cout << "Successful." << endl;

	if(verbose) cout << "Creating compatible ElementSpace" << endl;
	Tpetra::ElementSpace<ORDINALTYPE> const* bes2es = bes2.generateCompatibleElementSpace();
	if(debug) cout << (*bes2es);
	delete bes2es;

#ifdef TPETRA_MPI
  MPI_Finalize();
#endif // TPETRA_MPI
  
	if(verbose) cout << "BlockElementSpace testing successful." << endl;
  return(0); 
}
