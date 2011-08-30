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

#include "MueLu_ExtendedHashtable.hpp"


int main(int argc, char *argv[]) {
  using Teuchos::RCP;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // Handles some I/O to the output screen
  RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

  Teuchos::RCP<MueLu::TentativePFactory< > > PtentFact  = Teuchos::rcp(new MueLu::TentativePFactory<>());
  Teuchos::RCP<MueLu::TentativePFactory< > > PtentFact2 = Teuchos::rcp(new MueLu::TentativePFactory<>());
  Teuchos::RCP<MueLu::TentativePFactory< > > PtentFact3 = Teuchos::rcp(new MueLu::TentativePFactory<>());

  // build Operator
  Teuchos::ParameterList params;
  const RCP<const Map> map = MapFactory::Build(Xpetra::UseTpetra, 20, 0, comm);
  RCP<Operator> Op = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, CrsOperator>("Laplace1D", map, params);

  // an extended hashtable
  RCP<MueLu::UTILS::ExtendedHashtable> exh = Teuchos::rcp(new MueLu::UTILS::ExtendedHashtable());

  exh->Set<RCP<Operator> >("op",Op,PtentFact);
  RCP<Operator> test = exh->Get<RCP<Operator> > ("op",PtentFact);
  if(test!=Op) cout << "error" << endl;

  exh->Set<RCP<Operator> >("op",Op,PtentFact2);
  test = exh->Get<RCP<Operator> > ("op",PtentFact2);
  if(test!=Op) cout << "error" << endl;

  exh->Set("op2",24,Teuchos::null);
  int test2 = exh->Get<int>("op2",Teuchos::null);
  if(test2 != 24) cout << "error" << endl;

  exh->Set("op2",12,Teuchos::null);
  test2 = exh->Get<int>("op2",Teuchos::null);
  if(test2 != 12) cout << "error" << endl;

  exh->Remove("op",PtentFact2);
  exh->Remove("op",PtentFact);

  exh->Set<std::string>("op","xxx",PtentFact3);
  std::string test3 = exh->Get<std::string> ("op",PtentFact3);
  if(test3!="xxx") cout << "error" << endl;

  return EXIT_SUCCESS;
}


