// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"

int main(int argc, char * argv[]) 
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;
 
   Teuchos::GlobalMPISession mpiSession(&argc,&argv,&std::cout);

   std::string output_file_name = "square_mesh.gen";

   int xBlocks=1,yBlocks=1, zBlocks=1;
   int xElements=1,yElements=1, zElements=1;
   double x0=0.0, xf=1.0;
   double y0=0.0, yf=1.0;
   double z0=0.0, zf=1.0;
   bool threeD = false;
 
   // setup input arguments
   { 
     Teuchos::CommandLineProcessor clp;
     clp.throwExceptions(false);

     clp.setOption("o", &output_file_name, "Mesh output filename");

     clp.setOption("3d", "2d", &threeD, "Cube versus square mesh.");

     clp.setOption("x-blocks", &xBlocks, "Number of blocks in 'x' direction");
     clp.setOption("y-blocks", &yBlocks, "Number of blocks in 'y' direction");
     clp.setOption("z-blocks", &zBlocks, "Number of blocks in 'z' direction");

     clp.setOption("x-elmts", &xElements, "Number of elements in 'x' direction in each block");
     clp.setOption("y-elmts", &yElements, "Number of elements in 'y' direction in each block");
     clp.setOption("z-elmts", &zElements, "Number of elements in 'z' direction in each block");

     clp.setOption("x0", &x0, "Location of left edge");
     clp.setOption("xf", &xf, "Location of right edge");

     clp.setOption("y0", &y0, "Location of left edge");
     clp.setOption("yf", &yf, "Location of right edge");

     clp.setOption("z0", &z0, "Location of front(?) edge");
     clp.setOption("zf", &zf, "Location of back(?) edge");
 
     Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv,&std::cerr);

     if(parse_return==Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
       return -1;
      
     TEUCHOS_TEST_FOR_EXCEPTION(parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL, 
                                std::runtime_error, "Failed to parse command line!");
   }
  

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",xBlocks);
   pl->set("Y Blocks",yBlocks);
   pl->set("X Elements",xElements);
   pl->set("Y Elements",yElements);
   pl->set("X0",x0);
   pl->set("Y0",y0);
   pl->set("Xf",xf);
   pl->set("Yf",yf);

   if(threeD) {
     pl->set("Z Blocks",zBlocks);
     pl->set("Z Elements",zElements);
     pl->set("Z0",z0);
     pl->set("Zf",zf);
   }

   // int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
   // int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);

   RCP<panzer_stk::STK_MeshFactory> factory; 
   if(!threeD)
     factory = Teuchos::rcp(new panzer_stk::SquareQuadMeshFactory); 
   else
     factory = Teuchos::rcp(new panzer_stk::CubeHexMeshFactory); 

   factory->setParameterList(pl);

   RCP<panzer_stk::STK_Interface> mesh = factory->buildMesh(MPI_COMM_WORLD);
   mesh->writeToExodus(output_file_name);

   return 0;
}
