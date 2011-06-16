// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov),
//                    Denis Ridzal  (dridzal@sandia.gov),
//                    Kara Peterson (kjpeter@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   example_Poisson_stk.cpp
    \brief  Example solution of a Poisson equation on a hexahedral or 
            tetrahedral mesh using nodal (Hgrad) elements.

            This example requires a hexahedral or tetrahedral mesh in Exodus
            format with a nodeset containing boundary nodes. STK is used to
            read the mesh and populate a mesh database, Intrepid is used to
            build the stiffness matrix and right-hand side, and ML is used
            to solve the resulting linear system.

    \verbatim

     Poisson system:
 
            div A grad u = f in Omega
                       u = g on Gamma

       where
             A is a symmetric, positive definite material tensor
             f is a given source term


     Corresponding discrete linear system for nodal coefficients(x):

                 Kx = b

            K - HGrad stiffness matrix
            b - right hand side vector

    \endverbatim

    \author Created by P. Bochev, D. Ridzal, K. Peterson C. Siefert.

    \remark Usage:
    \code   ./example_Poisson_stk <meshfile>  \endcode

    \remark Example requires a hexahedral or tetrahedral mesh in Exodus format
*/
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

#define DUMP_DATA

/**************************************************************/
/*                          Includes                          */
/**************************************************************/

// Intrepid includes
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_Basis.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HGRAD_TET_C1_FEM.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"

// Epetra includes
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_Import.h"

// Teuchos includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

// AztecOO includes
#include "AztecOO.h"

// ML includes
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"

// Sacado includes
#include "Sacado.hpp"

// STK includes
#include "Ionit_Initializer.h"
#include "stk_io/IossBridge.hpp"
#include "stk_io/MeshReadWriteUtils.hpp"
#include "stk_util/parallel/Parallel.hpp"
#include "stk_mesh/base/FieldData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Comm.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/GetBuckets.hpp"
#include "stk_mesh/fem/CreateAdjacentEntities.hpp"

/*********************************************************/
/*                     Typedefs                          */
/*********************************************************/
typedef Sacado::Fad::SFad<double,3>      Fad3; //# ind. vars fixed at 3
typedef shards::CellTopology             ShardsCellTopology;
typedef Intrepid::FunctionSpaceTools     IntrepidFSTools;
typedef Intrepid::RealSpaceTools<double> IntrepidRSTools;
typedef Intrepid::CellTools<double>      IntrepidCTools;


/**********************************************************************************/
/******** FUNCTION DECLARATIONS FOR EXACT SOLUTION AND SOURCE TERMS ***************/
/**********************************************************************************/

/** \brief  User-defined exact solution.

    \param  x           [in]    x-coordinate of the evaluation point
    \param  y           [in]    y-coordinate of the evaluation point
    \param  z           [in]    z-coordinate of the evaluation point

    \return Value of the exact solution at (x,y,z)
 */
template<typename Scalar>
const Scalar exactSolution(const Scalar& x, const Scalar& y, const Scalar& z);

/** \brief  User-defined material tensor.

    \param  material    [out]   3 x 3 material tensor evaluated at (x,y,z)
    \param  x           [in]    x-coordinate of the evaluation point
    \param  y           [in]    y-coordinate of the evaluation point
    \param  z           [in]    z-coordinate of the evaluation point

    \warning Symmetric and positive definite tensor is required for every (x,y,z).
*/
template<typename Scalar>
void materialTensor(Scalar material[][3], const Scalar&  x, const Scalar&  y, const Scalar&  z);

/** \brief  Computes gradient of the exact solution. Requires user-defined exact solution.

    \param  gradExact  [out]   gradient of the exact solution evaluated at (x,y,z)
    \param  x          [in]    x-coordinate of the evaluation point
    \param  y          [in]    y-coordinate of the evaluation point
    \param  z          [in]    z-coordinate of the evaluation point
 */
template<typename Scalar>
void exactSolutionGrad(Scalar gradExact[3], const Scalar& x, const Scalar& y, const Scalar& z);


/** \brief Computes source term: f = -div(A.grad u).  Requires user-defined exact solution
           and material tensor. 

    \param  x          [in]    x-coordinate of the evaluation point
    \param  y          [in]    y-coordinate of the evaluation point
    \param  z          [in]    z-coordinate of the evaluation point

    \return Source term corresponding to the user-defined exact solution evaluated at (x,y,z)
 */
template<typename Scalar>
const Scalar sourceTerm(Scalar& x, Scalar& y, Scalar& z);


/** \brief Computation of the material tensor at array of points in physical space.

    \param worksetMaterialValues      [out]     Rank-2, 3 or 4 array with dimensions (C,P), (C,P,D) or (C,P,D,D)
                                                with the values of the material tensor
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void evaluateMaterialTensor(ArrayOut &        worksetMaterialValues,
                            const ArrayIn &   evaluationPoints);


/** \brief Computation of the source term at array of points in physical space.

    \param sourceTermValues           [out]     Rank-2 (C,P) array with the values of the source term
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void evaluateSourceTerm(ArrayOut &       sourceTermValues,
                        const ArrayIn &  evaluationPoints);

/** \brief Computation of the exact solution at array of points in physical space.

    \param exactSolutionValues        [out]     Rank-2 (C,P) array with the values of the exact solution
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void evaluateExactSolution(ArrayOut &       exactSolutionValues,
                           const ArrayIn &  evaluationPoints);


/** \brief Computation of the gradient of the exact solution at array of points in physical space.

    \param exactSolutionGradValues    [out]     Rank-3 (C,P,D) array with the values of the gradient of the exact solution
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void evaluateExactSolutionGrad(ArrayOut &       exactSolutionGradValues,
                               const ArrayIn &  evaluationPoints);


/**********************************************************************************/
/**************** FUNCTION DECLARATION FOR ML PRECONDITIONER *********************/
/**********************************************************************************/
int TestMultiLevelPreconditioner(char                        ProblemType[],
                                 Teuchos::ParameterList &    MLList,
                                 Epetra_CrsMatrix &          A,
                                 const Epetra_MultiVector &  xexact,
                                 Epetra_MultiVector &        b,
                                 Epetra_MultiVector &        uh,
                                 double &                    TotalErrorResidual,
                                 double &                    TotalErrorExactSol);


/**********************************************************************************/
/************* FUNCTION DECLARATIONS FOR SIMPLE BASIS FACTORY *********************/
/**********************************************************************************/

/** \brief  Simple factory that chooses basis function based on cell topology.

    \param  cellTopology  [in]    Shards cell topology
    \param  order         [in]    basis function order, currently unused
    \param  basis         [out]   pointer to Intrepid basis

    \return Intrepid basis 
 */

void getBasis(Teuchos::RCP<Intrepid::Basis<double,Intrepid::FieldContainer<double> > > &basis, 
               const shards::CellTopology & cellTopology,
               int order); 

/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/



/**********************************************************************************/
/******************************** MAIN ********************************************/
/**********************************************************************************/

int main(int argc, char *argv[]) {

  int numProcs=1;
  int rank=0;

#ifdef HAVE_MPI
  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);
  rank=mpiSession.getRank();
  numProcs=mpiSession.getNProc();
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();
  Epetra_Time Time(Comm);

   //Check number of arguments
  if ( (argc == 1) | (argc > 3)) {
     if (MyPID == 0) {
       std::cout <<"\n>>> ERROR: Invalid number of arguments.\n\n";
       std::cout <<"Usage:\n\n";
       std::cout <<"  ./Poisson <meshfile.exo> <solver.xml> \n\n";
       std::cout <<"   meshfile.exo            - Exodus hexhedral or tetrahedral mesh \n";
       std::cout <<"                             file with nodeset defined for boundary. \n\n";
       std::cout <<"   solver.xml (optional)   - xml file with ML solver options\n\n";
     }
     exit(1);
   }

 if (MyPID == 0){
  std::cout \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|              Example: Solve Poisson Equation                                |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
    << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" \
    << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  STK's website:      http://trilinos.sandia.gov/packages/stk                |\n" \
    << "|  ML's website:       http://trilinos.sandia.gov/packages/ml                 |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";
  }


#ifdef HAVE_MPI
  if (numProcs > 1) {
     if (MyPID == 0) {
       std::cout <<"\n>>> ERROR: Example will only run on a single processor. \n\n";
     }
     exit(1);
  }
#endif

/**********************************************************************************/
/********************************** GET XML INPUTS ********************************/
/**********************************************************************************/

   // get xml file from command line if provided, otherwise use default
    std::string  xmlSolverInFileName;
    if(argc>=3) xmlSolverInFileName=string(argv[2]);

   // Read xml file into parameter list
    Teuchos::ParameterList inputSolverList;

    if(xmlSolverInFileName.length()) {

       if (MyPID == 0)
        std::cout << "\nReading parameter list from the XML file \""<<xmlSolverInFileName<<"\" ...\n\n";

       Teuchos::updateParametersFromXmlFile(xmlSolverInFileName,&inputSolverList);

     } 
     else if (MyPID == 0) std::cout << "Using default solver values ..." << std::endl;


/**********************************************************************************/
/*********************************** READ MESH ************************************/
/**********************************************************************************/

   // 3-D meshes only
    int spaceDim = 3;

   // initialize io
    Ioss::Init::Initializer io;

   // define meta data
    stk::mesh::fem::FEMMetaData femMetaData(spaceDim);
    stk::mesh::MetaData &metaData = stk::mesh::fem::FEMMetaData::get_meta_data(femMetaData);

   // read in mesh from file
    stk::io::MeshData meshData;
    stk::io::create_input_mesh("exodusii",argv[1],MPI_COMM_WORLD,femMetaData,meshData);

   // commit meta data
    femMetaData.commit();

   // populate mesh entities (nodes, elements, etc.)
    stk::mesh::BulkData  bulkData(metaData,MPI_COMM_WORLD);
    stk::io::populate_bulk_data(bulkData, meshData);

  /*  Not necessary for Poisson problem
   // create adjacent entities
     stk::mesh::PartVector empty_add_parts;
     stk::mesh::create_adjacent_entities(bulkData, empty_add_parts);
  */

   // get entity ranks
    const stk::mesh::EntityRank elementRank = femMetaData.element_rank();
    const stk::mesh::EntityRank nodeRank    = femMetaData.node_rank();

   // get nodes
    std::vector<stk::mesh::Entity*> nodes;
    stk::mesh::get_entities(bulkData, nodeRank, nodes);
    int numNodes = nodes.size();

   // get elems
    std::vector<stk::mesh::Entity*> elems;
    stk::mesh::get_entities(bulkData, elementRank, elems);
    int numElems = elems.size();

    if (MyPID == 0) {
       std::cout << " Number of Elements: " << numElems << " \n";
       std::cout << "    Number of Nodes: " << numNodes << " \n\n";
    }

   // get coordinates field
    stk::mesh::Field<double, stk::mesh::Cartesian> *coords = 
       femMetaData.get_field<stk::mesh::Field<double, stk::mesh::Cartesian> >("coordinates");
    
   // get buckets containing entities of node rank
    const std::vector<stk::mesh::Bucket*> & nodeBuckets = bulkData.buckets( nodeRank );
    std::vector<stk::mesh::Entity*> bcNodes;

   // loop over all mesh parts 
    const stk::mesh::PartVector & all_parts = femMetaData.get_parts();
    for (stk::mesh::PartVector::const_iterator i  = all_parts.begin();
                                                i != all_parts.end(); ++i) {

         stk::mesh::Part & part = **i ;

        // if part only contains nodes, then it is a node set 
        //   ! this assumes that the only node set defined is the set 
        //   ! of boundary nodes
         if (part.primary_entity_rank() == nodeRank) {
            stk::mesh::Selector bcNodeSelector(part);
            stk::mesh::get_selected_entities(bcNodeSelector, nodeBuckets, bcNodes);
         }

    } // end loop over mesh parts

   // if no boundary node set was found give a warning
    if (bcNodes.size() == 0) {
       if (MyPID == 0) {
         std::cout << "\n     Warning! - No boundary node set found. \n";
         std::cout << "  Boundary conditions will not be applied correctly. \n\n";
       }
    }

    if(MyPID==0) {std::cout << "Read mesh                                   "
                   << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}

/**********************************************************************************/
/********************************SET CELL TOPOLOGY ********************************/
/**********************************************************************************/

   // Here the topology is defined from the mesh. Note that it is assumed
   // that there is a part labeled "block_1" and that the cell topology is
   // homogeneous over the entire mesh (i.e. there is not another block
   // containing elements with a different topology).

   // get the part labeled block_1
    stk::mesh::Part* const part = femMetaData.get_part("block_1");

   // get the topology of this part 
   //  (stk::mesh::fem::CellTopology is shards::CellTopology)
    stk::mesh::fem::CellTopology cellType = femMetaData.get_cell_topology( *part );

   // Get dimensions 
    int numNodesPerElem = cellType.getNodeCount();

    if(MyPID==0) {std::cout << "Get cell topology                           "
                   << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}

/**********************************************************************************/
/********************************* GET CUBATURE ***********************************/
/**********************************************************************************/

   // Define cubature of the specified degree for the cellType
    Intrepid::DefaultCubatureFactory<double>  cubFactory;
    int cubDegree = 2;
    Teuchos::RCP<Intrepid::Cubature<double> > cellCubature = cubFactory.create(cellType, cubDegree);

    int cubDim       = cellCubature -> getDimension();
    int numCubPoints = cellCubature -> getNumPoints();

   // Get numerical integration points and weights
    Intrepid::FieldContainer<double> cubPoints (numCubPoints, cubDim);
    Intrepid::FieldContainer<double> cubWeights(numCubPoints);

    cellCubature -> getCubature(cubPoints, cubWeights);

    if(MyPID==0) {std::cout << "Getting cubature                            " 
                 << Time.ElapsedTime() << " sec \n"  ; Time.ResetStartTime();}

/**********************************************************************************/
/*********************************** GET BASIS ************************************/
/**********************************************************************************/

   // Select basis from the cell topology
    int order = 1;
    Teuchos::RCP<Intrepid::Basis<double, Intrepid::FieldContainer<double> > >  HGradBasis;
    getBasis(HGradBasis, cellType, order);


    int numFieldsG = HGradBasis->getCardinality();
    Intrepid::FieldContainer<double> basisValues(numFieldsG, numCubPoints);
    Intrepid::FieldContainer<double> basisGrads(numFieldsG, numCubPoints, spaceDim);

  // Evaluate basis values and gradients at cubature points
    HGradBasis->getValues(basisValues, cubPoints, Intrepid::OPERATOR_VALUE);
    HGradBasis->getValues(basisGrads, cubPoints, Intrepid::OPERATOR_GRAD);

    if(MyPID==0) {std::cout << "Getting basis                               " 
                 << Time.ElapsedTime() << " sec \n"  ; Time.ResetStartTime();}


/**********************************************************************************/
/********************* BUILD MAPS FOR GLOBAL SOLUTION *****************************/
/**********************************************************************************/

  // For now assume local nodes = global nodes (and assume DOFs = nodes)

    // For parallel need a map like this:
    // Epetra_Map globalMapG(-1,ownedNodes,ownedGIDs,0,Comm);
 
    Epetra_Map   globalMapG(numNodes, 0, Comm);
    Epetra_FECrsMatrix StiffMatrix(Copy, globalMapG, numFieldsG);
    Epetra_FEVector rhsVector(globalMapG);
 
    if(MyPID==0) {std::cout <<  "Build global maps                           "
                 << Time.ElapsedTime() << " sec \n";  Time.ResetStartTime();}


#ifdef DUMP_DATA
/**********************************************************************************/
/**** PUT COORDINATES AND NODAL VALUES IN ARRAYS FOR OUTPUT (FOR PLOTTING ONLY) ***/
/**********************************************************************************/

  // Put coordinates in multivector for output
      Epetra_MultiVector nCoord(globalMapG,3);

  // Put element to node mapping in multivector for output
     Epetra_Map   globalMapElem(numElems, 0, Comm);
     Epetra_MultiVector elem2node(globalMapElem,numNodesPerElem);

    // Loop over elements
     for (size_t j = 0; j < elems.size(); j++) {

       // get nodes attached to this element 
        const stk::mesh::PairIterRelation elem_nodes = elems[j]->relations(nodeRank);
   
       // loop over nodes and fill element to node map
           // local ids
           //    element id : elems[j]->identifier()-1
           //       node id: elem_nodes[i].entity()->identifier()-1
        for (size_t i = 0; i < elem_nodes.size(); i++) {
          elem2node[i][elems[j]->identifier()-1] = elem_nodes[i].entity()->identifier() - 1;
          double * coord = stk::mesh::field_data(*coords, *elem_nodes[i].entity());
          nCoord[0][elem_nodes[i].entity()->identifier()-1] = coord[0];
          nCoord[1][elem_nodes[i].entity()->identifier()-1] = coord[1];
          nCoord[2][elem_nodes[i].entity()->identifier()-1] = coord[2];
        }

     } // end loop over elements

    // output multivectors
     EpetraExt::MultiVectorToMatrixMarketFile("elem2node.dat",elem2node,0,0,false);
     EpetraExt::MultiVectorToMatrixMarketFile("coords.dat",nCoord,0,0,false);

    if(MyPID==0) {Time.ResetStartTime();}

#endif

/**********************************************************************************/
/************************** DIRICHLET BC SETUP ************************************/
/**********************************************************************************/

 // Vector for use in applying BCs
    Epetra_MultiVector v(globalMapG,true);
    v.PutScalar(0.0);

    int * bcNodeVec = new int [bcNodes.size()];
    // Loop over boundary nodes
     for (unsigned i = 0; i < bcNodes.size(); i++) {
       
        int bcNodeId = bcNodes[i]->identifier() - 1;
        bcNodeVec[i] = bcNodeId;

        // get coordinates for this node
        stk::mesh::Entity * bcnode = bulkData.get_entity(nodeRank,bcNodes[i]->identifier());
        double * coord = stk::mesh::field_data(*coords, *bcnode);

        // look up exact value of function on boundary
        double x  = coord[0];
        double y  = coord[1];
        double z  = coord[2];
        v[0][bcNodeId]=exactSolution(x, y, z);

     } // end loop over boundary nodes

   if(MyPID==0) {std::cout << "Get Dirichlet boundary values               "
                 << Time.ElapsedTime() << " sec \n\n"; Time.ResetStartTime();}


/**********************************************************************************/
/******************** DEFINE WORKSETS AND LOOP OVER THEM **************************/
/**********************************************************************************/

   // Define desired workset size and count how many worksets there are on this processor's mesh block
    int desiredWorksetSize = numElems;                   // change to desired workset size!
    //int desiredWorksetSize = 100;                      // change to desired workset size!
    int numWorksets        = numElems/desiredWorksetSize;

    // When numElems is not divisible by desiredWorksetSize, increase workset count by 1
    if(numWorksets*desiredWorksetSize < numElems) numWorksets += 1;

    if (MyPID == 0) {
       std::cout << "\tDesired workset size:                 " << desiredWorksetSize <<"\n";
       std::cout << "\tNumber of worksets (per processor):   " << numWorksets <<"\n\n";
       Time.ResetStartTime();
     }

    for(int workset = 0; workset < numWorksets; workset++){

       // compute cell numbers where the workset starts and ends
       int worksetSize  = 0;
       int worksetBegin = (workset + 0)*desiredWorksetSize;
       int worksetEnd   = (workset + 1)*desiredWorksetSize;

       // when numElems is not divisible by desiredWorksetSize, the last workset ends at numElems
        worksetEnd   = (worksetEnd <= numElems) ? worksetEnd : numElems;
 
       // allocate the array for the cell nodes
        worksetSize  = worksetEnd - worksetBegin;
        Intrepid::FieldContainer<double> cellWorkset(worksetSize, numNodesPerElem, spaceDim);

       // copy coordinates into cell workset
        int cellCounter = 0;
        for(int cell = worksetBegin; cell < worksetEnd; cell++){

          // Get element entity from id of cell 
           stk::mesh::Entity * worksetElem = bulkData.get_entity(elementRank,cell+1);

          // get nodes attached to this element 
           const stk::mesh::PairIterRelation worksetNodes = worksetElem->relations(nodeRank);
   
          // loop over nodes and get coordinates to fill workset array
           for (size_t i = 0; i < worksetNodes.size(); i++) {
               double * coord = stk::mesh::field_data(*coords, *worksetNodes[i].entity());
               cellWorkset(cellCounter, i, 0) = coord[0];
               cellWorkset(cellCounter, i, 1) = coord[1];
               cellWorkset(cellCounter, i, 2) = coord[2];
           }
 
         cellCounter++;

      } // end cell loop

 /**********************************************************************************/
 /*                                Allocate arrays                                 */
 /**********************************************************************************/

     // Containers for Jacobians, integration measure & cubature points in workset cells
      Intrepid::FieldContainer<double> worksetJacobian  (worksetSize, numCubPoints, spaceDim, spaceDim);
      Intrepid::FieldContainer<double> worksetJacobInv  (worksetSize, numCubPoints, spaceDim, spaceDim);
      Intrepid::FieldContainer<double> worksetJacobDet  (worksetSize, numCubPoints);
      Intrepid::FieldContainer<double> worksetCubWeights(worksetSize, numCubPoints);
      Intrepid::FieldContainer<double> worksetCubPoints (worksetSize, numCubPoints, cubDim);

     // Containers for basis values transformed to workset cells and them multiplied by cubature weights
      Intrepid::FieldContainer<double> worksetBasisValues        (worksetSize, numFieldsG, numCubPoints);
      Intrepid::FieldContainer<double> worksetBasisValuesWeighted(worksetSize, numFieldsG, numCubPoints);
      Intrepid::FieldContainer<double> worksetBasisGrads         (worksetSize, numFieldsG, numCubPoints, spaceDim);
      Intrepid::FieldContainer<double> worksetBasisGradsWeighted (worksetSize, numFieldsG, numCubPoints, spaceDim);

     // Containers for diffusive & advective fluxes & non-conservative adv. term and reactive terms
      Intrepid::FieldContainer<double> worksetDiffusiveFlux(worksetSize, numFieldsG, numCubPoints, spaceDim);

     // Containers for material values and source term. Require user-defined functions
      Intrepid::FieldContainer<double> worksetMaterialVals (worksetSize, numCubPoints, spaceDim, spaceDim);
      Intrepid::FieldContainer<double> worksetSourceTerm   (worksetSize, numCubPoints);

     // Containers for workset contributions to the discretization matrix and the right hand side
      Intrepid::FieldContainer<double> worksetStiffMatrix (worksetSize, numFieldsG, numFieldsG);
      Intrepid::FieldContainer<double> worksetRHS         (worksetSize, numFieldsG);

     if(MyPID==0) {std::cout << "Allocate arrays                             "
                 << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}



 /**********************************************************************************/
 /*                                Calculate Jacobians                             */
 /**********************************************************************************/

      IntrepidCTools::setJacobian(worksetJacobian, cubPoints, cellWorkset, cellType);
      IntrepidCTools::setJacobianInv(worksetJacobInv, worksetJacobian );
      IntrepidCTools::setJacobianDet(worksetJacobDet, worksetJacobian );

      if(MyPID==0) {std::cout << "Calculate Jacobians                         "
                 << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}


 /**********************************************************************************/
 /*          Cubature Points to Physical Frame and Compute Data                    */
 /**********************************************************************************/

     // map cubature points to physical frame
      IntrepidCTools::mapToPhysicalFrame (worksetCubPoints, cubPoints, cellWorkset, cellType);

     // get A at cubature points
      evaluateMaterialTensor (worksetMaterialVals, worksetCubPoints);

     // get source term at cubature points
      evaluateSourceTerm (worksetSourceTerm, worksetCubPoints);

      if(MyPID==0) {std::cout << "Map to physical frame and get source term   "
                 << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}


 /**********************************************************************************/
 /*                         Compute Stiffness Matrix                               */
 /**********************************************************************************/

     // Transform basis gradients to physical frame:                        DF^{-T}(grad u)
      IntrepidFSTools::HGRADtransformGRAD<double>(worksetBasisGrads,                
                                                worksetJacobInv,   basisGrads);

     // Compute integration measure for workset cells:                      Det(DF)*w = J*w
      IntrepidFSTools::computeCellMeasure<double>(worksetCubWeights,              
                                                worksetJacobDet, cubWeights);


     // Multiply transformed (workset) gradients with weighted measure:     DF^{-T}(grad u)*J*w
      IntrepidFSTools::multiplyMeasure<double>(worksetBasisGradsWeighted,           
                                             worksetCubWeights, worksetBasisGrads);


     // Compute material tensor applied to basis grads:                     A*(DF^{-T}(grad u)
      IntrepidFSTools::tensorMultiplyDataField<double>(worksetDiffusiveFlux,      
                                                     worksetMaterialVals,
                                                     worksetBasisGrads);

     // Integrate to compute contribution to global stiffness matrix:      (DF^{-T}(grad u)*J*w)*(A*DF^{-T}(grad u))
      IntrepidFSTools::integrate<double>(worksetStiffMatrix,                      
                                       worksetBasisGradsWeighted,
                                       worksetDiffusiveFlux, Intrepid::COMP_BLAS);

      if(MyPID==0) {std::cout << "Compute stiffness matrix                    "
                 << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}

 /**********************************************************************************/
 /*                                   Compute RHS                                  */
 /**********************************************************************************/

     // Transform basis values to physical frame:                        clones basis values (u)
      IntrepidFSTools::HGRADtransformVALUE<double>(worksetBasisValues,              
                                                 basisValues);

     // Multiply transformed (workset) values with weighted measure:     (u)*J*w
      IntrepidFSTools::multiplyMeasure<double>(worksetBasisValuesWeighted,        
                                             worksetCubWeights, worksetBasisValues);

     // Integrate worksetSourceTerm against weighted basis function set:  f.(u)*J*w
      IntrepidFSTools::integrate<double>(worksetRHS,                             
                                       worksetSourceTerm,
                                       worksetBasisValuesWeighted, Intrepid::COMP_BLAS);

      if(MyPID==0) {std::cout << "Compute right-hand side                     "
                  << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}

 /**********************************************************************************/
 /*                         Assemble into Global Matrix                            */
 /**********************************************************************************/

     //"WORKSET CELL" loop: local cell ordinal is relative to numElems
      for(int cell = worksetBegin; cell < worksetEnd; cell++){

        // Compute cell ordinal relative to the current workset
         int worksetCellOrdinal = cell - worksetBegin;

        // Get element entity from id of cell 
         stk::mesh::Entity * worksetElem = bulkData.get_entity(elementRank,cell+1);

         // get nodes attached to this element 
          const stk::mesh::PairIterRelation worksetNodes = worksetElem->relations(nodeRank);
   
         // "CELL EQUATION" loop for the workset cell: cellRow is relative to the cell DoF numbering
          for (int cellRow = 0; cellRow < numFieldsG; cellRow++){

             int globalRow = worksetNodes[cellRow].entity()->identifier() - 1;
             double sourceTermContribution =  worksetRHS(worksetCellOrdinal, cellRow);
             rhsVector.SumIntoGlobalValues(1, &globalRow, &sourceTermContribution);

             // "CELL VARIABLE" loop for the workset cell: cellCol is relative to the cell DoF numbering
             for (int cellCol = 0; cellCol < numFieldsG; cellCol++){

               int globalCol = worksetNodes[cellCol].entity()->identifier() - 1;
               double operatorMatrixContribution = worksetStiffMatrix(worksetCellOrdinal, cellRow, cellCol);
               StiffMatrix.InsertGlobalValues(1, &globalRow, 1, &globalCol, &operatorMatrixContribution);
 
            }// end cell col loop 

         }// end cell row loop 

      }// end workset cell loop 

    } // end workset loop

    StiffMatrix.GlobalAssemble();
    StiffMatrix.FillComplete();
    rhsVector.GlobalAssemble();

    if(MyPID==0) {std::cout << "Global assembly                             "
                 << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}

/**********************************************************************************/
/************************ ADJUST MATRIX AND RHS FOR BCs ***************************/
/**********************************************************************************/

    // Apply stiffness matrix to v
     Epetra_MultiVector rhsDir(globalMapG,true);
     StiffMatrix.Apply(v,rhsDir);

    // Update right-hand side
     rhsVector.Update(-1.0,rhsDir,1.0);

    // Loop over boundary nodes and replace rhs values with boundary values
     for (size_t i = 0; i < bcNodes.size(); i++) {
       
        int bcNodeId = bcNodes[i]->identifier() - 1;
        rhsVector[0][bcNodeId]=v[0][bcNodeId];

     } // end loop over boundary nodes

     // Zero out rows and columns of stiffness matrix corresponding to Dirichlet edges
     //  and add one to diagonal.
      ML_Epetra::Apply_OAZToMatrix(bcNodeVec, bcNodes.size(), StiffMatrix);

     delete [] bcNodeVec;

     if(MyPID==0) {std::cout << "Adjust global matrix and rhs due to BCs     " << Time.ElapsedTime()
                  << " sec \n"; Time.ResetStartTime();}

#ifdef DUMP_DATA
    // Dump matrices to disk
     EpetraExt::RowMatrixToMatlabFile("stiff_matrix.dat",StiffMatrix);
     EpetraExt::MultiVectorToMatrixMarketFile("rhs_vector.dat",rhsVector,0,0,false);
#endif

/**********************************************************************************/
/*********************************** SOLVE ****************************************/
/**********************************************************************************/

    // Run the solver
     Teuchos::ParameterList MLList = inputSolverList;
     ML_Epetra::SetDefaults("SA", MLList, 0, 0, false);
     Epetra_FEVector exactNodalVals(globalMapG);
     Epetra_FEVector femCoefficients(globalMapG);
     double TotalErrorResidual = 0.0;
     double TotalErrorExactSol = 0.0;

    // Get exact solution at nodes
     for (unsigned inode = 0; inode < nodes.size(); inode++) {
       
        int nodeId = nodes[inode]->identifier() - 1;

        // get coordinates for this node
        stk::mesh::Entity * currentNode = bulkData.get_entity(nodeRank,nodes[inode]->identifier());
        double * coord = stk::mesh::field_data(*coords, *currentNode);

        // look up exact value of function 
        double x  = coord[0];
        double y  = coord[1];
        double z  = coord[2];
        exactNodalVals[0][nodeId]=exactSolution(x, y, z);

     } // end loop over nodes

     exactNodalVals.GlobalAssemble();

     char probType[10] = "laplace";

     TestMultiLevelPreconditioner(probType,             MLList,
                               StiffMatrix,          exactNodalVals,
                               rhsVector,            femCoefficients,
                               TotalErrorResidual,   TotalErrorExactSol);


/**********************************************************************************/
/**************************** CALCULATE ERROR *************************************/
/**********************************************************************************/

     if (MyPID == 0) {Time.ResetStartTime();}

     double L2err = 0.0;
     double L2errTot = 0.0;
     double H1err = 0.0;
     double H1errTot = 0.0;
     double Linferr = 0.0;
     double LinferrTot = 0.0;

#ifdef HAVE_MPI
    // Import solution onto current processor
     // FIXME
     int numNodesGlobal = globalMapG.NumGlobalElements();
     Epetra_Map     solnMap(numNodesGlobal, numNodesGlobal, 0, Comm);
     Epetra_Import  solnImporter(solnMap, globalMapG);
     Epetra_Vector  uCoeff(solnMap);
     uCoeff.Import(femCoefficients, solnImporter, Insert);
#endif

    // Define desired workset size
     desiredWorksetSize = numElems;
     int numWorksetsErr    = numElems/desiredWorksetSize;

    // When numElems is not divisible by desiredWorksetSize, increase workset count by 1
     if(numWorksetsErr*desiredWorksetSize < numElems) numWorksetsErr += 1;

    // Get cubature points and weights for error calc (may be different from previous)
     Intrepid::DefaultCubatureFactory<double>  cubFactoryErr;
     int cubDegErr = 3;
     Teuchos::RCP<Intrepid::Cubature<double> > cellCubatureErr = cubFactoryErr.create(cellType, cubDegErr);
     int cubDimErr       = cellCubatureErr->getDimension();
     int numCubPointsErr = cellCubatureErr->getNumPoints();
     Intrepid::FieldContainer<double> cubPointsErr(numCubPointsErr, cubDimErr);
     Intrepid::FieldContainer<double> cubWeightsErr(numCubPointsErr);
     cellCubatureErr->getCubature(cubPointsErr, cubWeightsErr);
 
    // Evaluate basis values and gradients at cubature points
     Intrepid::FieldContainer<double> uhGVals(numFieldsG, numCubPointsErr);
     Intrepid::FieldContainer<double> uhGrads(numFieldsG, numCubPointsErr, spaceDim);
     HGradBasis->getValues(uhGVals, cubPointsErr, Intrepid::OPERATOR_VALUE);
     HGradBasis->getValues(uhGrads, cubPointsErr, Intrepid::OPERATOR_GRAD);

    // Loop over worksets
     for(int workset = 0; workset < numWorksetsErr; workset++){

      // compute cell numbers where the workset starts and ends
       int worksetSize  = 0;
       int worksetBegin = (workset + 0)*desiredWorksetSize;
       int worksetEnd   = (workset + 1)*desiredWorksetSize;

      // when numElems is not divisible by desiredWorksetSize, the last workset ends at numElems
       worksetEnd   = (worksetEnd <= numElems) ? worksetEnd : numElems;

      // now we know the actual workset size and can allocate the array for the cell nodes
       worksetSize  = worksetEnd - worksetBegin;
       Intrepid::FieldContainer<double> cellWorksetEr(worksetSize, numNodesPerElem, spaceDim);
       Intrepid::FieldContainer<double> worksetApproxSolnCoef(worksetSize, numNodesPerElem);

      // loop over cells to fill arrays with coordinates and calculation solution coefficient
        int cellCounter = 0;
        for(int cell = worksetBegin; cell < worksetEnd; cell++){

          // Get element entity from id of cell 
           stk::mesh::Entity * worksetElem = bulkData.get_entity(elementRank,cell+1);

          // get nodes attached to this element 
           const stk::mesh::PairIterRelation worksetNodes = worksetElem->relations(nodeRank);
   
          // loop over nodes and get coordinates to fill workset array
           for (size_t i = 0; i < worksetNodes.size(); i++) {
               double * coord = stk::mesh::field_data(*coords, *worksetNodes[i].entity());
               cellWorksetEr(cellCounter, i, 0) = coord[0];
               cellWorksetEr(cellCounter, i, 1) = coord[1];
               cellWorksetEr(cellCounter, i, 2) = coord[2];
               int rowIndex = worksetNodes[i].entity()->identifier() - 1;
#ifdef HAVE_MPI
               worksetApproxSolnCoef(cellCounter, i) = uCoeff.Values()[rowIndex];
#else
               worksetApproxSolnCoef(cellCounter, i) = femCoefficients.Values()[rowIndex];
#endif
           }
 
         cellCounter++;

       } // end cell loop

      // Containers for Jacobian
       Intrepid::FieldContainer<double> worksetJacobianE(worksetSize, numCubPointsErr, spaceDim, spaceDim);
       Intrepid::FieldContainer<double> worksetJacobInvE(worksetSize, numCubPointsErr, spaceDim, spaceDim);
       Intrepid::FieldContainer<double> worksetJacobDetE(worksetSize, numCubPointsErr);
       Intrepid::FieldContainer<double> worksetCubWeightsE(worksetSize, numCubPointsErr);

      // Containers for basis values and gradients in physical space
       Intrepid::FieldContainer<double> uhGValsTrans(worksetSize,numFieldsG, numCubPointsErr);
       Intrepid::FieldContainer<double> uhGradsTrans(worksetSize, numFieldsG, numCubPointsErr, spaceDim);

      // compute cell Jacobians, their inverses and their determinants
       IntrepidCTools::setJacobian(worksetJacobianE, cubPointsErr, cellWorksetEr, cellType);
       IntrepidCTools::setJacobianInv(worksetJacobInvE, worksetJacobianE );
       IntrepidCTools::setJacobianDet(worksetJacobDetE, worksetJacobianE );

      // map cubature points to physical frame
       Intrepid::FieldContainer<double> worksetCubPoints(worksetSize, numCubPointsErr, cubDimErr);
       IntrepidCTools::mapToPhysicalFrame(worksetCubPoints, cubPointsErr, cellWorksetEr, cellType);

      // evaluate exact solution and gradient at cubature points
       Intrepid::FieldContainer<double> worksetExactSoln(worksetSize, numCubPointsErr);
       Intrepid::FieldContainer<double> worksetExactSolnGrad(worksetSize, numCubPointsErr, spaceDim);
       evaluateExactSolution(worksetExactSoln, worksetCubPoints);
       evaluateExactSolutionGrad(worksetExactSolnGrad, worksetCubPoints);

      // transform basis values to physical coordinates
       IntrepidFSTools::HGRADtransformVALUE<double>(uhGValsTrans, uhGVals);
       IntrepidFSTools::HGRADtransformGRAD<double>(uhGradsTrans, worksetJacobInvE, uhGrads);

      // compute weighted measure
       IntrepidFSTools::computeCellMeasure<double>(worksetCubWeightsE, worksetJacobDetE, cubWeightsErr);

      // evaluate the approximate solution and gradient at cubature points
       Intrepid::FieldContainer<double> worksetApproxSoln(worksetSize, numCubPointsErr);
       Intrepid::FieldContainer<double> worksetApproxSolnGrad(worksetSize, numCubPointsErr, spaceDim);
       IntrepidFSTools::evaluate<double>(worksetApproxSoln, worksetApproxSolnCoef, uhGValsTrans);
       IntrepidFSTools::evaluate<double>(worksetApproxSolnGrad, worksetApproxSolnCoef, uhGradsTrans);

      // get difference between approximate and exact solutions
       Intrepid::FieldContainer<double> worksetDeltaSoln(worksetSize, numCubPointsErr);
       Intrepid::FieldContainer<double> worksetDeltaSolnGrad(worksetSize, numCubPointsErr, spaceDim);
       IntrepidRSTools::subtract(worksetDeltaSoln, worksetApproxSoln, worksetExactSoln);
       IntrepidRSTools::subtract(worksetDeltaSolnGrad, worksetApproxSolnGrad, worksetExactSolnGrad);

      // take absolute values
       IntrepidRSTools::absval(worksetDeltaSoln);
       IntrepidRSTools::absval(worksetDeltaSolnGrad);

      // apply cubature weights to differences in values and grads for use in integration
       Intrepid::FieldContainer<double> worksetDeltaSolnWeighted(worksetSize, numCubPointsErr);
       Intrepid::FieldContainer<double> worksetDeltaSolnGradWeighted(worksetSize, numCubPointsErr, spaceDim);
       IntrepidFSTools::scalarMultiplyDataData<double>(worksetDeltaSolnWeighted, 
                                                worksetCubWeightsE, worksetDeltaSoln);
       IntrepidFSTools::scalarMultiplyDataData<double>(worksetDeltaSolnGradWeighted, 
                                                worksetCubWeightsE, worksetDeltaSolnGrad);

      // integrate to get errors on each element
       Intrepid::FieldContainer<double> worksetL2err(worksetSize);
       Intrepid::FieldContainer<double> worksetH1err(worksetSize);
       IntrepidFSTools::integrate<double>(worksetL2err, worksetDeltaSoln, 
                                          worksetDeltaSolnWeighted, Intrepid::COMP_BLAS);
       IntrepidFSTools::integrate<double>(worksetH1err, worksetDeltaSolnGrad, 
                                          worksetDeltaSolnGradWeighted, Intrepid::COMP_BLAS);

      // loop over cells to get errors for total workset
       cellCounter = 0;
       for(int cell = worksetBegin; cell < worksetEnd; cell++){

          // loop over cubature points
           for(int nPt = 0; nPt < numCubPointsErr; nPt++){

               Linferr = std::max(Linferr, worksetDeltaSoln(cellCounter,nPt));

            }

             L2err += worksetL2err(cellCounter); 
             H1err += worksetH1err(cellCounter);
 
         cellCounter++;

      } // end cell loop

    } // end loop over worksets


#ifdef HAVE_MPI
     // sum over all processors
     Comm.SumAll(&L2err,&L2errTot,1);
     Comm.SumAll(&H1err,&H1errTot,1);
     Comm.MaxAll(&Linferr,&LinferrTot,1);
#else
     L2errTot = L2err;
     H1errTot = H1err;
     LinferrTot = Linferr;
#endif

     if (MyPID == 0) {
       std::cout << "\n" << "L2 Error:  " << sqrt(L2errTot) <<"\n";
       std::cout << "H1 Error:  " << sqrt(H1errTot) <<"\n";
       std::cout << "LInf Error:  " << LinferrTot <<"\n\n";
     }


    if(MyPID==0) {std::cout << "Calculate error                             "
                  << Time.ElapsedTime() << " s \n"; Time.ResetStartTime();}


   return 0;

}
/**********************************************************************************/
/********************************* END MAIN ***************************************/
/**********************************************************************************/

/**********************************************************************************/
/************ USER DEFINED FUNCTIONS FOR EXACT SOLUTION ***************************/
/**********************************************************************************/

template<typename Scalar>
const Scalar exactSolution(const Scalar& x, const Scalar& y, const Scalar& z) {

  // Patch test - tet: function is in the FE space and should be recovered
     return 1. + x + y + z ;

  // Patch test - hex: tri-linear function is in the FE space and should be recovered
  // return 1. + x + y + z + x*y + x*z + y*z + x*y*z;

  // Analytic solution with homogeneous Dirichlet boundary data for [0 1]x[0 1]x[0 1]
  // return sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z)*exp(x+y+z);

  // Analytic solution with inhomogeneous Dirichlet boundary data
  // return exp(x + y + z)/(1. + x*y + y*z + x*y*z);
}


template<typename Scalar>
void materialTensor(Scalar material[][3], const Scalar& x, const Scalar& y, const Scalar& z) {

  material[0][0] = 1.;
  material[0][1] = 0.;
  material[0][2] = 0.;
  //
  material[1][0] = 0.;
  material[1][1] = 1.;
  material[1][2] = 0.;
  //
  material[2][0] = 0.;
  material[2][1] = 0.;
  material[2][2] = 1.;
}

/**********************************************************************************/
/************** AUXILIARY FUNCTIONS FROM EXACT SOLUTION ***************************/
/**********************************************************************************/

/************ Grad of Exact Solution ****************/
template<typename Scalar>
void exactSolutionGrad(Scalar gradExact[3], const Scalar& x, const Scalar& y, const Scalar& z) {

  // To enable derivatives of the gradient (i.e., 2nd derivatives of the exact solution) need 2 levels of fad types
  Sacado::Fad::SFad<Scalar,3> fad_x = x;
  Sacado::Fad::SFad<Scalar,3> fad_y = y;
  Sacado::Fad::SFad<Scalar,3> fad_z = z;
  Sacado::Fad::SFad<Scalar,3> u;

  // Indicate the independent variables
  fad_x.diff(0,3);
  fad_y.diff(1,3);
  fad_z.diff(2,3);

  u = exactSolution(fad_x, fad_y, fad_z);

  gradExact[0] = u.dx(0);
  gradExact[1] = u.dx(1);
  gradExact[2] = u.dx(2);
}

/************ Source Term (RHS) ****************/
template<typename Scalar>
const Scalar sourceTerm(Scalar& x, Scalar& y, Scalar& z){

  Scalar u;
  Scalar grad_u[3];
  Scalar flux[3];
  Scalar material[3][3];
  Scalar f = 0.;

  // Indicate the independent variables
  x.diff(0,3);
  y.diff(1,3);
  z.diff(2,3);

  // Get exact solution and its gradient
  u = exactSolution(x, y, z);
  exactSolutionGrad(grad_u, x, y, z);

  // Get material tensor
  materialTensor<Scalar>(material, x, y, z);

  // Compute total flux = (A.grad u)
  for(int i = 0; i < 3; i++){

    // Add diffusive flux
    for(int j = 0; j < 3; j++){
      flux[i] += material[i][j]*grad_u[j];
    }
  }

  // Compute source term (right hand side): f = -div(A.grad u) 
  f = -(flux[0].dx(0) + flux[1].dx(1) + flux[2].dx(2));

  return f;
}

/**********************************************************************************/
/*************************** EVALUATION METHODS ***********************************/
/**********************************************************************************/

/************ Material Tensor ****************/
template<class ArrayOut, class ArrayIn>
void evaluateMaterialTensor(ArrayOut &        matTensorValues,
                             const ArrayIn &   evaluationPoints){

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints        = evaluationPoints.dimension(1);
  int spaceDim         = evaluationPoints.dimension(2);

  double material[3][3];

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);
      double z = evaluationPoints(cell, pt, 2);

      materialTensor<double>(material, x, y, z);

      for(int row = 0; row < spaceDim; row++){
        for(int col = 0; col < spaceDim; col++){
          matTensorValues(cell, pt, row, col) = material[row][col];
        }
      }
    }
  }
}

/************ Source Term (RHS) ****************/
template<class ArrayOut, class ArrayIn>
void evaluateSourceTerm(ArrayOut &       sourceTermValues,
                        const ArrayIn &  evaluationPoints){

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints = evaluationPoints.dimension(1);

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      Sacado::Fad::SFad<double,3> x = evaluationPoints(cell, pt, 0);
      Sacado::Fad::SFad<double,3> y = evaluationPoints(cell, pt, 1);
      Sacado::Fad::SFad<double,3> z = evaluationPoints(cell, pt, 2);

      sourceTermValues(cell, pt) = sourceTerm<Sacado::Fad::SFad<double,3> >(x, y, z).val();
    }
  }
}

/************ Exact Solution ****************/
template<class ArrayOut, class ArrayIn>
void evaluateExactSolution(ArrayOut &       exactSolutionValues,
                           const ArrayIn &  evaluationPoints){

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints = evaluationPoints.dimension(1);

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);
      double z = evaluationPoints(cell, pt, 2);

      exactSolutionValues(cell, pt) = exactSolution<double>(x, y, z);
    }
  }
}


/************ Grad of Exact Solution ****************/
template<class ArrayOut, class ArrayIn>
void evaluateExactSolutionGrad(ArrayOut &       exactSolutionGradValues,
                               const ArrayIn &  evaluationPoints){

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints = evaluationPoints.dimension(1);
  int spaceDim  = evaluationPoints.dimension(2);

  double gradient[3];

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);
      double z = evaluationPoints(cell, pt, 2);

      exactSolutionGrad<double>(gradient, x, y, z);

      for(int row = 0; row < spaceDim; row++){
        exactSolutionGradValues(cell, pt, row) = gradient[row];
      }
    }
  }
}

/**********************************************************************************/
/******************************* TEST ML ******************************************/
/**********************************************************************************/

// Test ML
int TestMultiLevelPreconditioner(char ProblemType[],
                                 Teuchos::ParameterList   & MLList,
                                 Epetra_CrsMatrix   & A,
                                 const Epetra_MultiVector & xexact,
                                 Epetra_MultiVector & b,
                                 Epetra_MultiVector & uh,
                                 double & TotalErrorResidual,
                                 double & TotalErrorExactSol)
{
  Epetra_MultiVector x(xexact);
  x.PutScalar(0.0);

  Epetra_LinearProblem Problem(&A,&x,&b);
  Epetra_MultiVector* lhs = Problem.GetLHS();
  Epetra_MultiVector* rhs = Problem.GetRHS();

  Epetra_Time Time(A.Comm());

  // =================== //
  // call ML and AztecOO //
  // =================== //

  AztecOO solver(Problem);
  ML_Epetra::MultiLevelPreconditioner *MLPrec = new ML_Epetra::MultiLevelPreconditioner(A, MLList, true);

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, 1);

  solver.Iterate(200, 1e-10);

  delete MLPrec;

  uh = *lhs;

  // ==================================================== //
  // compute difference between exact solution and ML one //
  // ==================================================== //
  double d = 0.0, d_tot = 0.0;
  for( int i=0 ; i<lhs->Map().NumMyElements() ; ++i )
    d += ((*lhs)[0][i] - xexact[0][i]) * ((*lhs)[0][i] - xexact[0][i]);

  A.Comm().SumAll(&d,&d_tot,1);

  // ================== //
  // compute ||Ax - b|| //
  // ================== //
  double Norm;
  Epetra_Vector Ax(rhs->Map());
  A.Multiply(false, *lhs, Ax);
  Ax.Update(1.0, *rhs, -1.0);
  Ax.Norm2(&Norm);

  string msg = ProblemType;

  if (A.Comm().MyPID() == 0) {
    cout << msg << std::endl << "......Using " << A.Comm().NumProc() << " processes" << std::endl;
    cout << msg << "......||A x - b||_2 = " << Norm << std::endl;
    cout << msg << "......||x_exact - x||_2 = " << sqrt(d_tot) << std::endl;
    cout << msg << "......Total Time = " << Time.ElapsedTime() << std::endl << std::endl;
  }

  TotalErrorExactSol += sqrt(d_tot);
  TotalErrorResidual += Norm;

  return( solver.NumIters() );

}

/**********************************************************************************/
/**************************** SIMPLE BASIS FACTORY ********************************/
/**********************************************************************************/


void getBasis(Teuchos::RCP<Intrepid::Basis<double,Intrepid::FieldContainer<double> > > &basis, 
               const shards::CellTopology & cellTopology,
               int order)  {


 // select basis based on cell topology only for now, and assume first order basis
    switch (cellTopology.getKey()) {
       
       case shards::Tetrahedron<4>::key:
         basis = Teuchos::rcp(new Intrepid::Basis_HGRAD_TET_C1_FEM<double, Intrepid::FieldContainer<double> > );
         break;

       case shards::Hexahedron<8>::key:
         basis = Teuchos::rcp(new Intrepid::Basis_HGRAD_HEX_C1_FEM<double, Intrepid::FieldContainer<double> > );
         break;

       default:
         TEST_FOR_EXCEPTION( ( (cellTopology.getKey() != shards::Tetrahedron<4>::key)             &&
                               (cellTopology.getKey() != shards::Hexahedron<8>::key) ),
                                std::invalid_argument,
                               "Unknown cell topology for basis selction. Please use Hexahedron_8 or Tetrahedron_4.");

     }

}

  
     


