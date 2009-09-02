//@HEADER
// ************************************************************************
// 
//               Tpetra: Linear Algebra Services Package 
//                 Copyright (2009) Sandia Corporation
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
// ************************************************************************
//@HEADER

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <string>
#include <cmath>
#include <vector>
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#ifdef Tpetra_MPI
#include "mpi.h"
#include "Teuchos_MpiComm.hpp"
#endif
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#include "Tpetra_Version.hpp"
#include "Teuchos_RCP.hpp"
#include "Kokkos_DefaultNode.hpp"

//#define getNumMyElements getNumMyEntries
//#define getNumGlobalElements getNumGlobalEntries
//#define getMyGlobalElements getMyGlobalEntries
#define isNodeGlobalElement isMyGlobalIndex
#define getNodeElementList getElementList
#define getNodeNumElements getLocalNumElements
// prototype
template<class ST,class MT,class OT>
void power_method(Tpetra::CrsMatrix<ST,OT>& A, ST & lambda, OT niters, MT tolerance, bool verbose);


int main(int argc, char *argv[]) {
	typedef double ST;
	typedef Teuchos::ScalarTraits<ST>::magnitudeType MT;
	typedef int OT;
	
	
	
#ifdef Tpetra_MPI
	
	// Initialize MPI
	
	MPI_Init(&argc,&argv);
	
	Teuchos::RCP<Teuchos::Comm<OT> > comm = Teuchos::rcp(new Teuchos::MpiComm<OT>( RCP<OpaqueWrapper<MPI_Comm> &MPI_COMM_WORLD ));
	
#else
	
	Teuchos::RCP<const Teuchos::Comm<OT> > comm = Teuchos::rcp(new Teuchos::SerialComm<OT>());
	
#endif
	
	OT myRank = comm->getRank();
	OT numProc = comm->getSize();
	bool verbose = (myRank==0);
	
	if (verbose)
		std::cout << Tpetra::version() << std::endl << std::endl;
	
	std::cout << comm << std::endl;
	
	// Get the number of local equations from the command line
	if (argc!=2)
	{
		if (verbose) 
			std::cout << "Usage: " << argv[0] << " number_of_equations" << std::endl;
		std::exit(1);
	}
	OT numGlobalElements = std::atoi(argv[1]);
	
	if (numGlobalElements < numProc)
	{
		if (verbose)
			std::cout << "numGlobalBlocks = " << numGlobalElements 
			<< " cannot be < number of processors = " << numProc << std::endl;
		std::exit(1);
	}
	
	// Construct a Map that puts approximately the same number of 
	// equations on each processor.
	
	Tpetra::Map<OT> map(numGlobalElements, 0, comm);
	
	// Get update list and number of local equations from newly created map.
	
	size_t numMyElements = map.getNodeNumElements();
	
	Teuchos::ArrayView<const OT> myGlobalElements = map.getNodeElementList();
	
	// Create an OTeger vector NumNz that is used to build the Petra Matrix.
	// NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
	// on this processor
	
    std::vector<size_t> NumNz(numMyElements);
	
	// We are building a tridiagonal matrix where each row has (-1 2 -1)
	// So we need 2 off-diagonal terms (except for the first and last equation)
	
	for (size_t i=0; i<numMyElements; i++)
		if (myGlobalElements[i]==0 || myGlobalElements[i] == numGlobalElements-1)
			NumNz[i] = 2;
		else
			NumNz[i] = 3;
	
	// Create a Tpetra::Matrix
	
	Tpetra::CrsMatrix<ST,OT> A(map, Teuchos::ArrayView<size_t>(&NumNz[0], numMyElements));
	
	// Add  rows one-at-a-time
	// Need some vectors to help
	// Off diagonal values will always be -1
	
	
	std::vector<ST> values(2);
	values[0] = -1.0; values[1] = -1.0;
	std::vector<OT> indices(2);
	ST two = 2.0;
	OT numEntries;
	
	for (size_t i=0; i<numMyElements; i++)
    {
		if (myGlobalElements[i]==0)
		{
			indices[0] = 1;
			numEntries = 1;
		}
		else if (myGlobalElements[i] == numGlobalElements-1)
		{
			indices[0] = numGlobalElements-2;
			numEntries = 1;
		}
		else
		{
			indices[0] = myGlobalElements[i]-1;
			indices[1] = myGlobalElements[i]+1;
			numEntries = 2;
		}
		A.insertGlobalValues(myGlobalElements[i], Teuchos::ArrayView<OT>(&indices[0], numEntries), 
									Teuchos::ArrayView<ST>(&values[0], numEntries));
		// Put in the diagonal entry
		A.insertGlobalValues(myGlobalElements[i], Teuchos::ArrayView<const OT>(&myGlobalElements[i], 1), Teuchos::ArrayView<const ST>(&two, 1));
    }
	
	// Finish up
	A.fillComplete();
	
	// Create vectors for Power method
	
	
	// variable needed for iteration
	ST lambda = 0.0;
	OT niters = numGlobalElements*10;
	ST tolerance = 1.0e-2;
	
	// Iterate
	power_method<ST,MT,OT>(A, lambda, niters, tolerance, verbose);
	
	// Increase diagonal dominance
	if (verbose) 
		std::cout << "\nIncreasing magnitude of first diagonal term, solving again\n\n"
		<< std::endl;
	
	if (A.getRowMap().isNodeGlobalElement(0)) {
		size_t numvals = A.getNumEntriesForGlobalRow(0);
		std::vector<ST> Rowvals(numvals);
		std::vector<OT> Rowinds(numvals);
		A.getGlobalRowCopy(0, Teuchos::ArrayView<OT>(&Rowinds[0], numvals), Teuchos::ArrayView<ST>(&Rowvals[0], numvals), numvals); // Get A[0,0]
		for (size_t i=0; i<numvals; i++) if (Rowinds[i] == 0) Rowvals[i] *= 10.0;
		
		A.replaceGlobalValues(0, Teuchos::ArrayView<OT>(&Rowinds[0], numvals), Teuchos::ArrayView<ST>(&Rowvals[0], numvals));
	}
	
	// Iterate (again)
	lambda = 0.0;
	power_method<ST,MT,OT>(A, lambda, niters, tolerance, verbose);	
	
	// Release all objects
#ifdef Tpetra_MPI
	MPI_Finalize() ;
#endif
	
	/* end main
	 */
}
template<class ST,class MT,class OT>
void power_method(Tpetra::CrsMatrix<ST,OT> & A, ST &lambda, OT niters, MT tolerance, bool verbose) {  
	Kokkos::DefaultNode::DefaultNodeType node;
	Tpetra::Vector<ST,OT> q(node, A.getRowMap());
	Tpetra::Vector<ST,OT> z(node, A.getRowMap());
	Tpetra::Vector<ST,OT> resid(node, A.getRowMap());
	
	// Fill z with random Numbers
	z.randomize();
	
	// variable needed for iteration
	MT normz, residual;
	
	ST one = Teuchos::ScalarTraits<ST>::one();
	ST zero = Teuchos::ScalarTraits<ST>::zero();

	
	for (OT iter = 0; iter < niters; iter++)
    {
		normz = z.norm2(); // Compute 2-norm of z
		q.scale(one/normz, z);
		A.apply(q, z); // Compute z = A*q
		lambda = q.dot(z); // Approximate maximum eigenvalue
		if (iter%100==0 || iter+1==niters)
		{
			resid.update(one, z, -lambda, q, zero); // Compute A*q - lambda*q
			residual = resid.norm2();
			if (verbose) std::cout << "Iter = " << iter << "  Lambda = " << lambda 
			    << "  Residual of A*q - lambda*q = " 
			    << residual << std::endl;
		} 
		if (residual < tolerance) {
			break;
		}
    }
	return;
}
