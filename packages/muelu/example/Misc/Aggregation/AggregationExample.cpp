// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>

// Xpetra
#include <Xpetra_Parameters.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsOperator.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraMatrixFactory.hpp>

// MueLu
#include "MueLu_Graph.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_UCAggregationFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

void dumpAggregates(Aggregates & aggregates);

int main(int argc, char *argv[]) {
  using Teuchos::RCP;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  /**********************************************************************************/
  /* SET TEST PARAMETERS                                                            */
  /**********************************************************************************/
  // Note: use --help to list available options.
  Teuchos::CommandLineProcessor clp(false);
  
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp); // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);       // manage parameters of xpetra
  
  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }
  
  matrixParameters.check();
  xpetraParameters.check();

  if (comm->getRank() == 0) {
    std::cout << xpetraParameters << matrixParameters;
  }

  /**********************************************************************************/
  /* CREATE INITAL MATRIX                                                           */
  /**********************************************************************************/
  const RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
  RCP<Operator> Op = Galeri::Xpetra::CreateCrsMatrix<SC, LO, GO, Map, CrsOperator>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList()); //TODO: Operator vs. CrsOperator
  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/

  int currentPrintLevel=10;
  int printFlag=6;
  
  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/
  
  if (comm->getRank() == 0 && printFlag < currentPrintLevel)
    printf("main() Aggregate_CoarsenUncoupled : \n");
 
  RCP<Graph> graph = rcp(new Graph(Op->getCrsGraph(), "Uncoupled"));
  
  RCP<UCAggregationFactory> AggFact = rcp(new UCAggregationFactory());
  AggFact->SetMinNodesPerAggregate(2);  
  AggFact->SetMaxNeighAlreadySelected(5);
  AggFact->SetOrdering(MueLu::AggOptions::GRAPH);
  //AggFact->SetPhase3AggCreation(0.5); TODO: error with explicit instantiation. I don't know why

#ifdef JG_TO_UPDATE
//   RCP<Aggregates> aggregates = rcp(new Aggregates(*graph, "UC")); 
//   AggFact->Build(*graph, *aggregates);
  
  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/
  
  RCP<LOVector> Final_ = LOVectorFactory::Build( aggregates->GetVertex2AggId()->getMap() );

  {
    Teuchos::ArrayRCP<LO> Final = Final_->getDataNonConst(0);
    Teuchos::ArrayRCP<const LO> vertex2AggId = aggregates->GetVertex2AggId()->getData(0);
    Teuchos::ArrayRCP<const LO> procWinner   = aggregates->GetProcWinner()->getData(0);

    for (size_t i = 0; i < aggregates->GetVertex2AggId()->getMap()->getNodeNumElements(); i++) 
      Final[i] = vertex2AggId[i] + procWinner[i]*1000;
  }

  if (comm->getRank() == 0)
      printf("finals\n");
  //cout << *Final_ << endl; sleep(2);


  RCP<Teuchos::FancyOStream> out = rcp(new Teuchos::FancyOStream(Teuchos::rcp(&std::cout,false)));

  Final_->describe(*out, Teuchos::VERB_EXTREME);

  // dumpAggregates(*aggregates);

#endif

  return EXIT_SUCCESS;
}


#include <iostream>
#include <fstream>

void dumpAggregates(Aggregates & aggregates) {
  using Teuchos::RCP;

  int myPid = aggregates.GetMap()->getComm()->getRank();

  Teuchos::ArrayRCP<const LO> vertex2AggId = aggregates.GetVertex2AggId()->getData(0);
  Teuchos::ArrayRCP<const LO> procWinner   = aggregates.GetProcWinner()->getData(0);
  size_t n = aggregates.GetVertex2AggId()->getMap()->getNodeNumElements();
  RCP<const Map> map = aggregates.GetVertex2AggId()->getMap();


  char filename[200];
  snprintf(filename, 200, "aggregates-%d.data", myPid);
  std::ofstream out(filename);
  if (!out.is_open())
    std::cout << "Unable to open file";
  else
    {

      for (size_t i = 0; i < n; i++) 
        out << map->getGlobalElement(i) << " " << vertex2AggId[i] << " " << procWinner[i] << std::endl;;
   
      out.close();
    }

}
