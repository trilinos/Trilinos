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

/*Paul
03-August-2002 BES tester. Initial writeup.
18-Oct-2002 Modified.
22-Oct-2002 Changed to test out BES/BESData friend wrinkles.
21-Jan-2003 Updated for .hpp
*/

#define ORDINALTYPE int

#include <iostream>
#include <iomanip>
#include "Tpetra_SerialPlatform.hpp" 
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_BlockElementSpace.hpp"

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
  
  // Platform
  if(verbose) cout << "Creating platform" << endl;
  Tpetra::SerialPlatform<ORDINALTYPE, ORDINALTYPE> platform;

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
  if(verbose) cout << "Creating bes, constructor1" << endl;
  Tpetra::BlockElementSpace<ORDINALTYPE> bes(es, ELEMENTSIZE);
  if(debug) cout << bes;

	if(verbose) cout << "Creating bes, constructor2" << endl;
	ORDINALTYPE* esizelist = new ORDINALTYPE[NUMELEMENTS];
	esizelist[0] = 1;
	esizelist[1] = 1;
	esizelist[2] = 2;
	esizelist[3] = 3;
	esizelist[4] = 5;
	Tpetra::BlockElementSpace<ORDINALTYPE> bes2(es, esizelist);
	delete[] esizelist;
	esizelist = 0;
	if(debug) cout << bes2;

	if(verbose) cout << "Creating compatible ElementSpace" << endl;
	Tpetra::ElementSpace<ORDINALTYPE> const* bes2es = bes2.generateCompatibleElementSpace();
	if(debug) cout << (*bes2es);
	delete bes2es;

	if(verbose) cout << "BlockElementSpace testing successful." << endl;
  return(0); 
}
