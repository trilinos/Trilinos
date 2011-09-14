// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_FancyOStream.hpp>

// Kokkos
#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>

// Xpetra
#include <Xpetra_Parameters.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsOperator.hpp>

// Gallery
#define XPETRA_ENABLED // == Gallery have to be build with the support of Xpetra matrices.
#include <MueLu_GalleryParameters.hpp>
#include <MueLu_MatrixFactory.hpp>

// MueLu
#include "MueLu_Hierarchy.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_IfpackSmoother.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

int main(int argc, char *argv[]) {
  using Teuchos::RCP;
  using std::cout;
  using std::endl;

  RCP<Teuchos::FancyOStream> fos = getFancyOStream(Teuchos::rcpFromRef(cout));

  //
  //
  //

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  Teuchos::CommandLineProcessor clp(false);
  
  MueLu::Gallery::Parameters<GO> matrixParameters(clp); // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);       // manage parameters of xpetra

  Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_DEFAULT;
  {
    Teuchos::EVerbosityLevel values[6] = {Teuchos::VERB_DEFAULT, Teuchos::VERB_NONE, Teuchos::VERB_LOW, Teuchos::VERB_MEDIUM, Teuchos::VERB_HIGH, Teuchos::VERB_EXTREME};
    const char* names[6] = {"VERB_DEFAULT", "VERB_NONE", "VERB_LOW", "VERB_MEDIUM", "VERB_HIGH", "VERB_EXTREME"};
    clp.setOption("verbLevel", &verbLevel, 6, values, names, "Verbose level");
    
  }
  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }
  
  matrixParameters.check();
  xpetraParameters.check();

  const RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
  RCP<Operator> A = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, CrsOperator>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());

  cout << endl << endl;
  //
  //
  //

  Level l;
  Hierarchy H;
  H.SetLevel(rcpFromRef(l));
  l.Request("A", NULL); //FIXME: remove line
  l.Set("A", A, NULL); //FIXME2: remove NULL

  //
  //
  //

  Teuchos::ParameterList smootherParamList;
  smootherParamList.set("relaxation: type", "symmetric Gauss-Seidel");
  smootherParamList.set("relaxation: sweeps", (LO) 1);
  smootherParamList.set("relaxation: damping factor", (SC) 1.0);
  
  IfpackSmoother ifpackSmoo("point relaxation stand-alone", smootherParamList);
  ifpackSmoo.setObjectLabel("My Ifpack Smoother");

  SmootherFactory smooFact(rcpFromRef(ifpackSmoo), rcpFromRef(ifpackSmoo));
  //  SmootherFactory smooFact(rcpFromRef(ifpackSmoo));

  smooFact.describe(*fos, verbLevel);

  //ifpackSmoo.describe(*fos, verbLevel);

  cout << endl << endl;  
  ifpackSmoo.Setup(l);
  cout << endl << endl;

  //ifpackSmoo.describe(*fos, verbLevel);

  //
  //
  //
  
  return EXIT_SUCCESS;

}

// setVerbLevel(Teuchos::VERB_HIGH);

//   cout << "---------- description() ---------- " << endl
//        << ifpackSmoo.description()
//        << endl;

//   cout << "---------- << operator ---------- " << endl
//        << ifpackSmoo;

//   cout << "---------- describe(NONE) ---------- " << endl;
//   ifpackSmoo.describe(*fos, Teuchos::VERB_NONE);

//   cout << "---------- describe(LOW) ---------- " << endl;
//   ifpackSmoo.describe(*fos, Teuchos::VERB_LOW);

//   cout << "---------- describe(MEDIUM) ---------- " << endl;
//   ifpackSmoo.describe(*fos, Teuchos::VERB_MEDIUM);

//   cout << "---------- describe(HIGH) ---------- " << endl;
//   ifpackSmoo.describe(*fos, Teuchos::VERB_HIGH);

//   cout << "---------- describe(EXTREME) ---------- " << endl;
//   ifpackSmoo.describe(*fos, Teuchos::VERB_EXTREME);

//   cout << "-----------" << endl;
