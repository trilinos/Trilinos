#include <string>

#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Tpetra_DefaultPlatform.hpp>
#include <MatrixMarket_Tpetra.hpp> // for loading matrices from file

#include "Amesos2_Factory.hpp"
#include "Amesos2_Util_is_same.hpp"
#include "Amesos2_Superlu.hpp"
#include "Amesos2_MatrixAdapter.hpp"
#include "Amesos2_MultiVecAdapter.hpp"

using std::cout;
using std::endl;
using std::string;

using Teuchos::as;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::Comm;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::tuple;
using Teuchos::ScalarTraits;
using Teuchos::OrdinalTraits;
using Teuchos::CONJ_TRANS;
using Teuchos::TRANS;
using Teuchos::NO_TRANS;


using Tpetra::global_size_t;
using Tpetra::CrsMatrix;
using Tpetra::MultiVector;
using Tpetra::Map;

typedef Tpetra::DefaultPlatform::DefaultPlatformType Platform;
typedef typename Platform::NodeType Node;

int main(int argc, char *argv[]){
  typedef float SCALAR;
  typedef int LO;
  typedef int GO;
  typedef std::complex<SCALAR> cmplx;
  typedef CrsMatrix<cmplx,LO,GO,Node> MAT;
  typedef ScalarTraits<cmplx> ST;
  typedef MultiVector<cmplx,LO,GO,Node> MV;
  typedef typename ST::magnitudeType Mag;
  typedef ScalarTraits<Mag> MT;
  const size_t numVecs = 1;

  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  RCP<const Comm<int> > comm = platform.getComm();
  RCP<Node>             node = platform.getNode();

  RCP<MAT> A =
    Tpetra::MatrixMarket::Reader<MAT>::readSparseFile("../test/matrices/amesos2_test_mat3.mtx",comm,node);

  RCP<const Map<LO,GO,Node> > rowmap = A->getRowMap();
  RCP<const Map<LO,GO,Node> > colmap = A->getColMap();

  RCP<MV> X = rcp(new MV(rowmap,
			 tuple<cmplx>(cmplx(-0.58657, 0.10646),
				      cmplx(0.86716, 0.75421),
				      cmplx(0.58970, 0.29876)),
			 3,
			 numVecs));
  RCP<MV> B = rcp(new MV(colmap,numVecs));
  RCP<MV> Xhat = rcp(new MV(rowmap,numVecs));
  X->setObjectLabel("X");
  B->setObjectLabel("B");
  Xhat->setObjectLabel("Xhat");

  A->apply(*X,*B,Teuchos::CONJ_TRANS); // use conjugate-transpose

  Xhat->randomize();

  // Solve A*Xhat = B for Xhat using the Superlu solver
  RCP<Amesos::SolverBase> solver
    = Amesos::Factory<MAT,MV>::create("Superlu",
				      A,
				      Xhat,
				      B );

  Teuchos::ParameterList params;
  //  params.set("Trans","TRANS","Solve with conjugate-transpose");

  solver->setParameters( rcpFromRef(params) );
  solver->symbolicFactorization().numericFactorization().solve();

  std::ostream &out = std::cout;
  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

#ifdef HAVE_AMESOS2_DEBUG
  B->describe(*fos, Teuchos::VERB_EXTREME);
  Xhat->describe(*fos, Teuchos::VERB_EXTREME);
  X->describe(*fos, Teuchos::VERB_EXTREME);
#endif

  // Check result of solve
  Array<Mag> xhatnorms(numVecs), xnorms(numVecs);
  Xhat->norm2(xhatnorms());
  X->norm2(xnorms());
}
