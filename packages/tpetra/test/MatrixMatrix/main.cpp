#include "Tpetra_MMMultiply.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Tpetra_Map.hpp"

using namespace Teuchos;
using namespace Tpetra;
int main(int argc, char* argv[]){
	Teuchos::oblackholestream blackhole;
	Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
	/*Array<int> aNN (tuple<int>(2,2,2,1));
	Array<int> bNN (tuple<int>(1,2,2,1));
	const ArrayRCP<const int> aNumNnz = arcpFromArray(aNN);
	const ArrayRCP<const int> bNumNnz =arcpFromArray(bNN);*/

	RCP<const Comm<int> > comm = DefaultPlatform::getDefaultPlatform().getComm();
	const RCP<const Tpetra::Map<int> > rowMap = rcp(new Tpetra::Map<int>(4, 0, comm));

	RCP<CrsMatrix<double, int> > matrixA = rcp(new CrsMatrix<double, int>(rowMap, 2));
	RCP<CrsMatrix<double, int> > matrixB = rcp(new CrsMatrix<double, int>(rowMap, 2));


	matrixA->insertGlobalValues(0, tuple<int>(0,1), tuple<double>(1,4));
	matrixA->insertGlobalValues(1, tuple<int>(1,2), tuple<double>(2,6));
	matrixA->insertGlobalValues(2, tuple<int>(2,3), tuple<double>(3,7));
	matrixA->insertGlobalValues(3, tuple<int>(3), tuple<double>(4));


	matrixB->insertGlobalValues(0, tuple<int>(0), tuple<double>(8));
	matrixB->insertGlobalValues(1, tuple<int>(0,1), tuple<double>(9,1));
	matrixB->insertGlobalValues(2, tuple<int>(1,2), tuple<double>(2,7));
	matrixB->insertGlobalValues(3, tuple<int>(3), tuple<double>(3));

	std::cout << "just before fillCompletes \n";
	matrixA->fillComplete(rowMap, rowMap);
	matrixB->fillComplete(rowMap, rowMap);

	std::cout << "just after fillCompletes \n";
	RCP<CrsMatrix<double, int> > matrixC = rcp(new CrsMatrix<double,int>(rowMap, 4));

	MatrixMatrixMultiply<double, int> multiplier(matrixA, matrixB, matrixC);
	multiplier.multiply();
	Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
	matrixC->describe(*out, Teuchos::VERB_EXTREME);

	return 0;
}

