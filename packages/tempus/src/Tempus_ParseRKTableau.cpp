#include "Tempus_ParseRKTableau.hpp"

#include "Tempus_String_Utilities.hpp"
#include "Tempus_RKButcherTableau.hpp"

#include <string>
#include <vector>

namespace Tempus {

Teuchos::RCP<RKButcherTableau<double> >
parseRKTableau(const Teuchos::ParameterList & pl)
{
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;


  std::string error_msg = "The format of the Butcher Tableau parameter list is\n"
                          "  <Parameter name=\"A\" type=\"string\" value=\"# # # ;\n"
                          "                                                # # # ;\n"
                          "                                                # # #\"/>\n"
                          "  <Parameter name=\"b\" type=\"string\" value=\"# # #\"/>\n"
                          "  <Parameter name=\"c\" type=\"string\" value=\"# # #\"/>\n\n"
                          "Note the number of stages is implicit in the number of entries.\n" 
                          "The number of stages must be consistent.\n";

  std::string embedded_error_msg = "The format of the Butcher Tableau parameter list is\n"
                          "  <Parameter name=\"A\" type=\"string\" value=\"# # # ;\n"
                          "                                                # # # ;\n"
                          "                                                # # #\"/>\n"
                          "  <Parameter name=\"b\" type=\"string\" value=\"# # #\"/>\n"
                          "  <Parameter name=\"bhat\" type=\"string\" value=\"# # #\"/>\n"
                          "  <Parameter name=\"c\" type=\"string\" value=\"# # #\"/>\n\n"
                          "Note the number of stages is implicit in the number of entries.\n" 
                          "The number of stages must be consistent.\n";

  // do some error checking
    bool has_A = pl.isType<std::string>("A");
    bool has_b = pl.isType<std::string>("b");
    bool has_c = pl.isType<std::string>("c");

    TEUCHOS_TEST_FOR_EXCEPTION(!has_A,std::runtime_error,error_msg);
    TEUCHOS_TEST_FOR_EXCEPTION(!has_b,std::runtime_error,error_msg);
    TEUCHOS_TEST_FOR_EXCEPTION(!has_c,std::runtime_error,error_msg);

  std::size_t numStages = -1;
  Teuchos::SerialDenseMatrix<int,double> A_tbl;
  Teuchos::SerialDenseVector<int,double> b_tbl;
  Teuchos::SerialDenseVector<int,double> bhat_tbl;
  Teuchos::SerialDenseVector<int,double> c_tbl;


  // read in the A matrix
  {
    std::vector<std::string> A_row_tokens; 
    Tempus::StringTokenizer(A_row_tokens,pl.get<std::string>("A"),";",true);

    // this is the only place where numStages is set
    numStages = A_row_tokens.size();

    // allocate the matrix
    A_tbl.shape(as<int>(numStages),as<int>(numStages));

    // fill the rows
    for(std::size_t r=0;r<numStages;r++) {
      // parse the row (tokenize on space)
      std::vector<std::string> tokens;
      Tempus::StringTokenizer(tokens,A_row_tokens[r]," ",true);

      std::vector<double> values;
      Tempus::TokensToDoubles(values,tokens);

      TEUCHOS_TEST_FOR_EXCEPTION(values.size()!=numStages,std::runtime_error,
                                 "Error parsing A matrix, wrong number of stages in row " << r << "\n" + error_msg);

      for(std::size_t c=0;c<numStages;c++)
        A_tbl(r,c) = values[c]; 
    }
  }

  // size b and c vectors
  b_tbl.size(as<int>(numStages));
  c_tbl.size(as<int>(numStages));


  // read in the b vector
  {
    std::vector<std::string> tokens; 
    Tempus::StringTokenizer(tokens,pl.get<std::string>("b")," ",true);
    std::vector<double> values;
    Tempus::TokensToDoubles(values,tokens);

    TEUCHOS_TEST_FOR_EXCEPTION(values.size()!=numStages,std::runtime_error,
                               "Error parsing b vector, wrong number of stages.\n" + error_msg);

    for(std::size_t i=0;i<numStages;i++)
      b_tbl(i) = values[i]; 
  }

  // read in the c vector
  {
    std::vector<std::string> tokens; 
    Tempus::StringTokenizer(tokens,pl.get<std::string>("c")," ",true);
    std::vector<double> values;
    Tempus::TokensToDoubles(values,tokens);

    TEUCHOS_TEST_FOR_EXCEPTION(values.size()!=numStages,std::runtime_error,
                               "Error parsing c vector, wrong number of stages.\n" + error_msg);

    for(std::size_t i=0;i<numStages;i++)
      c_tbl(i) = values[i]; 
  }


  RCP<RKButcherTableau<double> > tableau = rcp(new RKButcherTableau<double>);

  tableau->initialize(A_tbl,b_tbl,c_tbl,1111111,"Built from drekar::parseRKTableau (Embedded) (order unspecified, not 1111111!)");

  return tableau;
}

}
