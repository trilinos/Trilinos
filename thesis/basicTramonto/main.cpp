#include "Optika_GUI.hpp"

int main(int argc, char* argv[]){
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::sublist;
  RCP<ParameterList> tramontoList = 
    rcp(new ParameterList("Tramonto"));
  RCP<ParameterList> dimList = 
    sublist(tramontoList, "Dimension Control Parameters");

  dimList->set( "Length_ref", -1.0, 
    "Indicates how length parameters will be entered");

  dimList->set( "Density_ref", -1.0, 
    "Indicates how density parameters will be entered");

  dimList->set("Temp", -1.0, 
    "Temperature to be used. Set to -1.0 if using reduced units");

  dimList->set("Dielec_ref", -1.0, "The dielectric constant");

  dimList->set("VEXT_MAX", 0, 
    "The maximum external field where DFT equations will be solved.");

  Optika::getInput(dimList);
  return 0;
}
