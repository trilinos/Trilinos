/*
 * testExtendedHashtable.cpp
 *
 *  Created on: 28.08.2011
 *      Author: tobias
 */



#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_Hashtable.hpp>
#include <Teuchos_HashUtils.hpp>

// Xpetra
#include <Xpetra_Parameters.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsOperator.hpp>

// Gallery
#define XPETRA_ENABLED // == Gallery have to be build with the support of Xpetra matrices.
#include <MueLu_GalleryParameters.hpp>
#include <MueLu_MatrixFactory.hpp>

#include "MueLu_Hierarchy.hpp"  // this is absolutely necessary???
#include "MueLu_SaPFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

#include "ExtendedHashtable.hpp"

#include "Teuchos_TabularOutputter.hpp"

int main(int argc, char *argv[]) {
  using Teuchos::RCP;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // Handles some I/O to the output screen
  RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

  //Teuchos::RCP<MueLu::UCAggregationFactory> UCAggFact = Teuchos::rcp(new UCAggregationFactory());
  Teuchos::RCP<MueLu::TentativePFactory< > > PtentFact = Teuchos::rcp(new MueLu::TentativePFactory<>());
  //Teuchos::RCP<MueLu::SaPFactory>       Pfact = Teuchos::rcp( new SaPFactory(TentPFact) );

  // build Operator
  Teuchos::ParameterList params;
  const RCP<const Map> map = MapFactory::Build(Xpetra::UseTpetra, 20, 0, comm);
  RCP<Operator> Op = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, CrsOperator>("Laplace1D", map, params);


  // an extended hashtable
  RCP<MueLu::UTILS::ExtendedHashtable> exh = Teuchos::rcp(new MueLu::UTILS::ExtendedHashtable());

  exh->Set<RCP<Operator> >("Hallo",Op,PtentFact.get());
  RCP<Operator> test = exh->Get<RCP<Operator> > ("Hallo",PtentFact.get());

  exh->Set("TEMPO",24,NULL);
  int test2 = exh->Get<int>("TEMPO",NULL);
  cout << test2 << endl;

  exh->Set("TEMPO",12,NULL);
  test2 = exh->Get<int>("TEMPO",NULL);
  cout << test2 << endl;

  exh->Print(*out); *out << endl;

  exh->Set<RCP<Operator> >("Hallo",Op,NULL);
  RCP<Operator> test3 = exh->Get<RCP<Operator> > ("Hallo",NULL);

  exh->Print(*out); *out << endl;

  exh->Remove("Hallo",NULL);

  exh->Print(*out); *out << endl;

  exh->Remove("Hallo",PtentFact.get());

  exh->Print(*out); *out << endl;

  exh->Remove("TEMPO",NULL);

  exh->Print(*out); *out << endl;

  return EXIT_SUCCESS;
}


