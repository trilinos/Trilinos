#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include "Teuchos_VerboseObject.hpp"
#include <Teuchos_FancyOStream.hpp>

// Cthulhu
#include <Cthulhu_Parameters.hpp>
#include <Cthulhu_Map.hpp>
#include <Cthulhu_MapFactory.hpp>
#include <Cthulhu_CrsOperator.hpp>
#include <Cthulhu_Example.hpp>

// Gallery
#define CTHULHU_ENABLED // == Gallery have to be build with the support of Cthulhu matrices.
#include <MueLu_GalleryParameters.hpp>
#include <MueLu_MatrixFactory.hpp>

// Aggregation
#include "MueLu_UCAggregationFactory.hpp"

#include "MueLu_UseShortNames.hpp"

void dumpAggregates(Aggregates & aggregates);

int main(int argc, char *argv[]) {
  
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  /**********************************************************************************/
  /* SET TEST PARAMETERS                                                            */
  /**********************************************************************************/
  // Note: use --help to list available options.
  Teuchos::CommandLineProcessor clp(false);
  
  MueLu::Gallery::Parameters matrixParameters(clp); // manage parameters of the test case
  Cthulhu::Parameters cthulhuParameters(clp);       // manage parameters of cthulhu
  
  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }
  
  matrixParameters.check();
  cthulhuParameters.check();

  if (comm->getRank() == 0) {
    matrixParameters.print();
    cthulhuParameters.print();
  }

  /**********************************************************************************/
  /* CREATE INITAL MATRIX                                                           */
  /**********************************************************************************/
  const RCP<const Map> map = MapFactory::Build(cthulhuParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
  RCP<CrsOperator> Op = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, CrsOperator>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList()); //TODO: Operator vs. CrsOperator
  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/

  int currentPrintLevel=10;
  int printFlag=6;
  
  MueLu::AggregationOptions aggOptions;
  
  aggOptions.SetPrintFlag(printFlag);      
  aggOptions.SetMinNodesPerAggregate(2);  
  aggOptions.SetMaxNeighAlreadySelected(5);
  // aggOptions.SetOrdering(1); //TODO: RandomReorder()
  aggOptions.SetOrdering(2);
  aggOptions.SetPhase3AggCreation(0.5);

  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/
  
  if (comm->getRank() == 0 && printFlag < currentPrintLevel)
    printf("main() Aggregate_CoarsenUncoupled : \n");
 
  RCP<Graph> graph = rcp(new Graph(Op->getCrsGraph(), "Uncoupled"));
  
  RCP<UCAggregationFactory> AggFact = rcp(new UCAggregationFactory(aggOptions));
  RCP<Aggregates> aggregates = AggFact->Build(*graph);
  
  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/
  
  RCP<Cthulhu::Vector<int> > Final_ = Cthulhu::VectorFactory<int>::Build( aggregates->GetVertex2AggId()->getMap() );

  {
    Teuchos::ArrayRCP<int> Final = Final_->getDataNonConst(0);
    Teuchos::ArrayRCP<const int> vertex2AggId = aggregates->GetVertex2AggId()->getData(0);
    Teuchos::ArrayRCP<const int> procWinner   = aggregates->GetProcWinner()->getData(0);

    for (size_t i = 0; i < aggregates->GetVertex2AggId()->getMap()->getNodeNumElements(); i++) 
      Final[i] = vertex2AggId[i] + procWinner[i]*1000;
  }

  if (comm->getRank() == 0)
      printf("finals\n");
  //cout << *Final_ << endl; sleep(2);


  RCP<Teuchos::FancyOStream> out = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcp(&std::cout,false)));

  Final_->describe(*out, Teuchos::VERB_EXTREME);

  // dumpAggregates(*aggregates);

  return EXIT_SUCCESS;
}


#include <iostream>
#include <fstream>

void dumpAggregates(Aggregates & aggregates) {
  int myPid = aggregates.GetMap()->getComm()->getRank();

  Teuchos::ArrayRCP<const int> vertex2AggId = aggregates.GetVertex2AggId()->getData(0);
  Teuchos::ArrayRCP<const int> procWinner   = aggregates.GetProcWinner()->getData(0);
  size_t n = aggregates.GetVertex2AggId()->getMap()->getNodeNumElements();
  RCP<const Map> map = aggregates.GetVertex2AggId()->getMap();


  char filename[200];
  snprintf(filename, 200, "aggregates-%d.data", myPid);
  std::ofstream out(filename);
  if (!out.is_open())
    cout << "Unable to open file";
  else
    {

      for (size_t i = 0; i < n; i++) 
        out << map->getGlobalElement(i) << " " << vertex2AggId[i] << " " << procWinner[i] << std::endl;;
   
      out.close();
    }

}
