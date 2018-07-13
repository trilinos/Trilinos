
#include <string>

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <MatrixMarket_Tpetra.hpp> // for loading matrices from file

#include "Amesos2.hpp"
#include "Amesos2_Meta.hpp"

#include <Kokkos_Core.hpp>

using namespace std;

int main(int argc, char* argv[])
{
  typedef int                        LO;
  typedef int                        GO;
  typedef double                     VAL;

  typedef Tpetra::DefaultPlatform::DefaultPlatformType Platform;
  typedef Platform::NodeType          Node;
  //typedef Kokkos::OpenMP                         Node;
  typedef Tpetra::Map<LO,GO,Node>                 Map;

  typedef Tpetra::CrsMatrix<VAL,LO,GO,Node>           MAT;
  typedef Tpetra::MultiVector<VAL,LO,GO,Node>         VEC;
  typedef Teuchos::Comm<int>                         Comm;

  typedef Teuchos::RCP<Teuchos::ParameterList>       PList;


  Teuchos::GlobalMPISession mpisession(&argc,&argv,&std::cout);

  Kokkos::initialize(argc, argv);
  

  std::string solver_name= string(argv[1]);
  char *matrix_name = argv[2];

  
  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  Teuchos::RCP<const Comm> comm = platform.getComm();
  Teuchos::RCP<Node>            node = platform.getNode();

  Teuchos::RCP<MAT> A = Tpetra::MatrixMarket::Reader<MAT>::readSparseFile(matrix_name, comm, node);
  
  Teuchos::RCP<const Map > dmnmap = A->getDomainMap();
  Teuchos::RCP<const Map > rngmap = A->getRangeMap();
  
  Teuchos::RCP<VEC> X = Teuchos::rcp(new VEC(dmnmap,1));
  Teuchos::RCP<VEC> Y = Teuchos::rcp(new VEC(rngmap,1));

  X->setObjectLabel("X");
  Y->setObjectLabel("Y");

  X->randomize();
  A->apply(*X,*Y);

  Teuchos::RCP<Amesos2::Solver<MAT,VEC> > solver
    = Amesos2::create<MAT, VEC>("Basker", A, X, Y);

  //Paradiso options
  /*
  PList myplist = Teuchos::rcp(new Teuchos::ParameterList());
  if(solver_name.compare("pardiso_mkl")==0)
    {
      myplist->set("IPARM(2)", 3); //there ND
    }
  */

  cout << "Before Sfactor" << endl;
  //Teuchos::Time::Time ts("stime",true);
  solver->symbolicFactorization();
  //ts.stop();
  //cout << "After SFactor, Time : " 
  //     << ts.totalElapsedTime() << endl;
  cout << "Before Numerical Factorization" << endl;
  //Teuchos::Time::Time t("time", true);
  solver->numericFactorization();
  //t.stop();
  //cout << "After Numerical Factorization, Time: " 
  //     << t.totalElapsedTime() << endl;
  solver->solve();

  Kokkos::finalize();

  return 0;
}
