#include "Tpetra_RowMatrixTransposer.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Teuchos_VerboseObject.hpp"


int main(int argc, char* argv[]){
	Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
	Teuchos::oblackholestream blackhole;
	Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
	typedef double Scalar;
	typedef Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
	typedef int Ordinal;
	using Tpetra::global_size_t;
	
	Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
	
	size_t myRank = comm->getRank();
	size_t numProc = comm->getSize();
	bool verbose = (myRank==0);
	
	std::cout << *comm;
	
	const global_size_t numGlobalElements = 4;
	if (numGlobalElements < numProc) {
		if (verbose) {
			std::cout << "numGlobalBlocks = " << numGlobalElements 
			<< " cannot be less than the number of processors = " << numProc << std::endl;
		}
		return -1;
	}
	
	// Construct a Map that puts approximately the same number of equations on each processor.
	
	Teuchos::RCP<const Tpetra::Map<Ordinal> > map = Tpetra::createUniformContigMap<Ordinal,Ordinal>(numGlobalElements, comm);
	
	// Get update list and number of local equations from newly created map.
	
	const size_t numMyElements = map->getNodeNumElements();
	
	Teuchos::ArrayView<const Ordinal> myGlobalElements = map->getNodeElementList();
	
	// Create an OTeger vector NumNz that is used to build the Petra Matrix.
	// NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
	// on this processor
	
	Teuchos::ArrayRCP<size_t> NumNz = Teuchos::arcp<size_t>(numMyElements);

	// We are building a tridiagonal matrix where each row has (-1 2 -1)
	// So we need 2 off-diagonal terms (except for the first and last equation)

	for (size_t i=0; i < numMyElements; ++i) {
		if (myGlobalElements[i] == 0 || static_cast<global_size_t>(myGlobalElements[i]) == numGlobalElements-1) {
		// boundary
			NumNz[i] = 2;
		}
		else {
			NumNz[i] = 3;
		}
	}

	// Create a Tpetra::Matrix using the Map, with a static allocation dictated by NumNz
	Tpetra::CrsMatrix<Scalar,Ordinal>  A (map, NumNz, Tpetra::StaticProfile);
	Tpetra::CrsMatrix<Scalar,Ordinal>  AT(map, NumNz, Tpetra::StaticProfile);
	Teuchos::RCP< Tpetra::CrsMatrix<Scalar,Ordinal> > TestMatrix = Teuchos::null;

	// We are done with NumNZ
	NumNz = Teuchos::null;

	// Add  rows one-at-a-time
	// Off diagonal values will always be -1
	const Scalar two    = static_cast<Scalar>( 2.0);
	const Scalar negOne = static_cast<Scalar>(-1.0);
	const Scalar three = static_cast<Scalar>(3.0);
	for (size_t i=0; i<numMyElements; i++) {
		if (myGlobalElements[i] == 0) {
			A.insertGlobalValues( myGlobalElements[i],
			Teuchos::tuple<Ordinal>( myGlobalElements[i], myGlobalElements[i]+1 ),
			Teuchos::tuple<Scalar> ( two, negOne ) );
		}
		else if (static_cast<global_size_t>(myGlobalElements[i]) == numGlobalElements-1) {
			A.insertGlobalValues( myGlobalElements[i],
			Teuchos::tuple<Ordinal>( myGlobalElements[i]-1, myGlobalElements[i] ),
			Teuchos::tuple<Scalar> ( negOne, two ) );
		}
		else {
			A.insertGlobalValues( myGlobalElements[i],
			Teuchos::tuple<Ordinal>( myGlobalElements[i]-1, myGlobalElements[i], myGlobalElements[i]+1 ),
			Teuchos::tuple<Scalar> ( three, two, negOne ) );
		}
	}

	/*
	for (size_t i=0; i<numMyElements; i++) {
		if (myGlobalElements[i] == 0) {
			AT->insertGlobalValues( myGlobalElements[i],
			Teuchos::tuple<Ordinal>( myGlobalElements[i], myGlobalElements[i]+1 ),
			Teuchos::tuple<Scalar> ( two, three ) );
		}
		else if (static_cast<global_size_t>(myGlobalElements[i]) == numGlobalElements-1) {
			AT->insertGlobalValues( myGlobalElements[i],
			Teuchos::tuple<Ordinal>( myGlobalElements[i]-1, myGlobalElements[i] ),
			Teuchos::tuple<Scalar> ( negOne, two ) );
		}
		else if(static_cast<global_size_t>(myGlobalElements[i])==1){
			AT->insertGlobalValues( myGlobalElements[i],
			Teuchos::tuple<Ordinal>( myGlobalElements[i]-1, myGlobalElements[i], myGlobalElements[i]+1 ),
			Teuchos::tuple<Scalar> ( negOne, two, three ) );
		}
		else if(static_cast<global_size_t>(myGlobalElements[i])==2){
			AT->insertGlobalValues( myGlobalElements[i],
			Teuchos::tuple<Ordinal>( myGlobalElements[i]-1, myGlobalElements[i], myGlobalElements[i]+1 ),
			Teuchos::tuple<Scalar> ( negOne, two, negOne ) );
		}
	}*/
	
	// Finish up
	A.fillComplete(Tpetra::DoOptimizeStorage);
	AT.fillComplete(Tpetra::DoOptimizeStorage);
	
	//	A->describe(*out, Teuchos::VERB_EXTREME); 
//		AT->describe(*out, Teuchos::VERB_EXTREME); 
	Tpetra::RowMatrixTransposer<Scalar, Ordinal> transposer = Tpetra::RowMatrixTransposer<Scalar, Ordinal>(A);
	/*Teuchos::RCP<Tpetra::Map<Ordinal> > tMap = Teuchos::rcp(new Tpetra::Map<Ordinal>(4, 0, comm));*/
	transposer.createTranspose(Tpetra::DoOptimizeStorage, TestMatrix/*, tMap*/);

//	if (verbose) {
	TestMatrix->describe(*out, Teuchos::VERB_EXTREME); 
	//}
	

	return 0;
}

