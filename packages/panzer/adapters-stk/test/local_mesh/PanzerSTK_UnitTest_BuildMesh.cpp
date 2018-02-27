#include "PanzerSTK_UnitTest_BuildMesh.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_DefaultComm.hpp"

#include "Tpetra_DefaultPlatform.hpp"

#include "Panzer_STK_Interface.hpp"

#include "PanzerSTK_UnitTest_STKInterfaceGenerator.hpp"

namespace panzer_stk
{

Teuchos::RCP<panzer_stk::STK_Interface>
buildMesh(const std::vector<int> & N,
          const std::vector<int> & B,
          const std::vector<double> & L)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm();

  const int num_dims = N.size();

  TEUCHOS_ASSERT(L.size() == static_cast<std::size_t>(num_dims));

  Teuchos::ParameterList mesh_parameters("Mesh");

  mesh_parameters.set<double>("X0",-L[0]/2.);
  mesh_parameters.set<double>("Xf", L[0]/2.);
  mesh_parameters.set<int>("X Elements", N[0]);
  mesh_parameters.set<int>("X Blocks", B[0]);
  if(num_dims>1){
    mesh_parameters.set<double>("Y0",-L[1]/2.);
    mesh_parameters.set<double>("Yf", L[1]/2.);
    mesh_parameters.set<int>("Y Elements", N[1]);
    mesh_parameters.set<int>("Y Blocks", B[1]);
  }
  if(num_dims>2){
    mesh_parameters.set<double>("Z0",-L[2]/2.);
    mesh_parameters.set<double>("Zf", L[2]/2.);
    mesh_parameters.set<int>("Z Elements", N[2]);
    mesh_parameters.set<int>("Z Blocks", B[2]);
  }

  if(num_dims == 1){
    mesh_parameters.set<std::string>("Mesh Type","Line");
  } else if(num_dims == 2){
    mesh_parameters.set<std::string>("Mesh Type","Quad");
  } else if(num_dims == 3){
    mesh_parameters.set<std::string>("Mesh Type","Hex");
  } else {
    TEUCHOS_ASSERT(false);
  }

  return generateMesh(comm, mesh_parameters);
}

}
