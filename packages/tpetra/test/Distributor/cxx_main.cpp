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

#define PACKETTYPE float
#define ORDINALTYPE int

#include "Tpetra_SerialDistributor.hpp"
#include "Tpetra_Version.hpp"

int main(int argc, char* argv[]) {
	// initialize verbose & debug flags
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

  if(verbose)
		cout << Tpetra::Tpetra_Version() << endl << endl;

  if(verbose) cout << "Creating SerialDistributor object...";
  Tpetra::SerialDistributor<PACKETTYPE, ORDINALTYPE> distributor;
  //if(debug) cout <<distributor.label() << endl;
	if(verbose) cout << "Successful." << endl;

	//  void createFromSends(const OrdinalType& numExportIDs, const OrdinalType* exportImageIDs,
	//											 const bool& deterministic, OrdinalType& numRemoteIDs ) 

	ORDINALTYPE nEIDs = 2;
	ORDINALTYPE* eIIDs = 0;
	bool determ = false;
	ORDINALTYPE nRIDs = 2;
	//distributor.createFromSends(nEIDs, eIIDs, false, nRIDs);
  
	if(verbose) cout << "Distributor test successful." << endl;
  
  return(0);
}
