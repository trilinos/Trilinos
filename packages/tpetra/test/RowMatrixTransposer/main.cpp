#include "Tpetra_DefaultPlatform.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Tpetra_MatrixMatrix.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include <cmath>



template<class CrsMatrix_t> double getNorm(CrsMatrix_t& matrix){
  double mySum = 0;
  Teuchos::Array<int> inds(matrix.getNodeMaxNumRowEntries());
  Teuchos::Array<double> vals(matrix.getNodeMaxNumRowEntries());
  for(int i =0; ((size_t)i)<matrix.getNodeNumRows(); ++i){
    size_t numRowEnts = matrix.getNumEntriesInLocalRow(i);
    Teuchos::ArrayView<const int> indsView = inds();
    Teuchos::ArrayView<const double> valsView = vals();
    matrix.getLocalRowView(i, indsView, valsView);
    for(size_t j=0; ((size_t)j)<numRowEnts; ++j){
      mySum += valsView[j]*valsView[j];
    }
  }
  double totalSum = 0;
  Teuchos::reduceAll(*(matrix.getComm()), Teuchos::REDUCE_SUM, 1, &mySum, &totalSum);
  return sqrt(totalSum);

}



int main(int argc, char* argv[]){
	Teuchos::oblackholestream blackhole;
	Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
	typedef double Scalar;
	typedef Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
	typedef int Ordinal;
	using Tpetra::global_size_t;
	
	Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcp(&std::cout,false));
  //out->setOutputToRootOnly(comm->getRank());
	
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

	
	for (size_t i=0; i<numMyElements; i++) {
		if (myGlobalElements[i] == 0) {
			AT.insertGlobalValues( myGlobalElements[i],
			Teuchos::tuple<Ordinal>( myGlobalElements[i], myGlobalElements[i]+1 ),
			Teuchos::tuple<Scalar> ( two, three ) );
		}
		else if (static_cast<global_size_t>(myGlobalElements[i]) == numGlobalElements-1) {
			AT.insertGlobalValues( myGlobalElements[i],
			Teuchos::tuple<Ordinal>( myGlobalElements[i]-1, myGlobalElements[i] ),
			Teuchos::tuple<Scalar> ( negOne, two ) );
		}
		else if(static_cast<global_size_t>(myGlobalElements[i])==1){
			AT.insertGlobalValues( myGlobalElements[i],
			Teuchos::tuple<Ordinal>( myGlobalElements[i]-1, myGlobalElements[i], myGlobalElements[i]+1 ),
			Teuchos::tuple<Scalar> ( negOne, two, three ) );
		}
		else if(static_cast<global_size_t>(myGlobalElements[i])==2){
			AT.insertGlobalValues( myGlobalElements[i],
			Teuchos::tuple<Ordinal>( myGlobalElements[i]-1, myGlobalElements[i], myGlobalElements[i]+1 ),
			Teuchos::tuple<Scalar> ( negOne, two, negOne ) );
		}
	}
	
	// Finish up
	A.fillComplete(Tpetra::DoOptimizeStorage);
	AT.fillComplete(Tpetra::DoOptimizeStorage);

	
	Tpetra::RowMatrixTransposer<Scalar, Ordinal> transposer(A);
	TestMatrix = transposer.createTranspose(Tpetra::DoOptimizeStorage); //, TestMatrix/*, tMap*/);

  Teuchos::RCP<Tpetra::CrsMatrix<Scalar, Ordinal> > diffMatrix = Tpetra::createCrsMatrix<Scalar, Ordinal>(TestMatrix->getRowMap());
  //Apparently there is a problem with ADD because while these two matricies are 
  //identical when I add them together I don't get 0 like I should. Infact
  //I just get a matrix that has the exact same entries and sparsity structure.
  //I'll have to come back to this later. But RowMatrixTransposer is working right.
  //And all my other tests are telling me ADD works right too.
  //KLN 06/14/2011
  Tpetra::MatrixMatrix::Add(AT,false,-1.0,*TestMatrix,false, 1.0,diffMatrix);
  diffMatrix->fillComplete();
  //diffMatrix->describe(*out, Teuchos::VERB_EXTREME);
  double diffNorm = getNorm(*diffMatrix);
  double realNorm = getNorm(AT);
  double epsilon = diffNorm/realNorm;
  if(epsilon > 1e-10){
    *out << "The calculated A transpose and the real one don't match!" << std::endl;
    *out << "Diff Norm: " << diffNorm << std::endl;
    *out << "Real Norm: " << realNorm << std::endl;
    *out << "Epsilon: " << epsilon << std::endl;
    return 1;
  }
  



	return 0;
}


