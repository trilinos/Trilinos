// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_Array.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "MatrixMarket_Tpetra.hpp"

#include "Amesos2.hpp"
#include "Amesos2_Meta.hpp"

void process_command_line(int argc, char*argv[], std::string& xml_file);


int main(int argc, char *argv[])
{

  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
  bool success = true;

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try
    {
      Teuchos::Time timer("total");
      timer.start();

      typedef double Scalar;
      typedef Tpetra::Map<>::local_ordinal_type LO;
      typedef Tpetra::Map<>::global_ordinal_type GO;
      typedef Tpetra::Map<>::node_type Node;
      typedef Tpetra::MultiVector<Scalar,LO,GO,Node> TMV;
      typedef Tpetra::Operator<Scalar,LO,GO,Node>    TOP;

      typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node>  crs_matrix_type;
      typedef Tpetra::Map< LO, GO, Node>               map_type;



      Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

      std::string xml_file("shylubasker_test1_mm.xml");
      process_command_line(argc, argv, xml_file);

      *out << "Every proc reading parameters from xml_file: "
	   << xml_file << std::endl;
      Teuchos::ParameterList test_params =
	Teuchos::ParameterXMLFileReader(xml_file).getParameters();
      //Debug
      *out << "Testing inot parameter list" << std::endl;
      test_params.print();



      std::string mm_file("matrix1.mm"); //Default matrix
      mm_file = test_params.get<std::string>("mm_file");
      int num_solves = 1;
      num_solves = test_params.get<int>("num_solves");
      Teuchos::ParameterList amesos2_params =
	test_params.get<Teuchos::ParameterList>("Amesos2");
      //Test
      *out << "Test output of Amesos2 parameters" << std::endl;
      amesos2_params.print();


      Teuchos::RCP<crs_matrix_type> A =
	Tpetra::MatrixMarket::Reader<crs_matrix_type>::readSparseFile(mm_file,
								      comm);

      *out << "Done reading matrix market file" << std::endl;


      //Make some Vectors
      Teuchos::RCP<const map_type> dmnmap = A->getDomainMap();
      Teuchos::RCP<const map_type> rngmap = A->getRangeMap();


      Teuchos::RCP<TMV>  X = Teuchos::rcp(new TMV(rngmap,num_solves));
      X->setObjectLabel("X");
      X->randomize();
      Teuchos::RCP<TMV>  B = Teuchos::rcp(new TMV(dmnmap,num_solves));
      B->setObjectLabel("B");
      Teuchos::RCP<TMV> Xhat = Teuchos::rcp(new TMV(rngmap,num_solves));
      Xhat->setObjectLabel("Xhat");

      A->apply(*X,*B,Teuchos::TRANS);

      //Create Basker Solver
      Teuchos::RCP< Amesos2::Solver<crs_matrix_type,TMV> > solver =
	Amesos2::create<crs_matrix_type, TMV>("Basker", A, Xhat, B);

      //Add Parameters
      solver->setParameters(Teuchos::rcpFromRef(amesos2_params));

      //Call Symbolic
      solver->symbolicFactorization();
      //Call Numeric
      solver->numericFactorization();
      //Call Solve
      solver->solve();

      //Check if good (Not the best way!!!)
      Teuchos::Array<Scalar> xhatnorms(num_solves);
      Teuchos::Array<Scalar> xnorms(num_solves);

      Xhat->norm2(xhatnorms());
      X->norm2(xnorms());

    }//end -- try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success)

  if(success)
    {
      *out << "End Result: TEST PASSED\n";
    }
  else
    {
      *out << "End Result: TEST FAILED\n";
    }

  return ( success ? 0 : 1);
}//end-main


void process_command_line(int argc, char*argv[], std::string& xml_file)
{
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("xml_file", &xml_file, "XML Parameters file");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    throw std::runtime_error("Error parsing command-line.");
  }
}


