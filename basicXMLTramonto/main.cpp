#include "Optika_GUI.hpp"

int main(int argc, char* argv[]){
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using Teuchos::rcp;
  RCP<ParameterList> userInput = rcp(new ParameterList); 
  Optika::getInput("inputs.xml", userInput);
  return 0;
}
