/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

void Original(int beg , int end)
{
	std::ifstream readfile;
	timeval start , stop , result;

	int numCells = 0;
	for(int i = beg ; i < end; i++) {	
		// Get cell topology for base hexahedron
		shards::CellTopology cellType(shards::getCellTopologyData<shards::Hexahedron<8> >() );
		typedef shards::CellTopology CellType;
		// Get dimensions
		int numNodesPerCell = cellType.getNodeCount();
		int spaceDim = cellType.getDimension();


		/*********************************************************/
		/*    Create hexahedral cell and copy into batch         */
		/*********************************************************/

		// define generic hexahedron
		Intrepid::FieldContainer<double> nodeCoords( numNodesPerCell, spaceDim);
		nodeCoords(0,0) = 0.0; nodeCoords(0,1) = 0.0;  nodeCoords(0,2) = 0;
		nodeCoords(1,0) = 1.2; nodeCoords(1,1) = -0.1; nodeCoords(1,2) = 0.2;
		nodeCoords(2,0) = 1.1; nodeCoords(2,1) = 1.1;  nodeCoords(2,2) = -0.1;
		nodeCoords(3,0) = 0.1; nodeCoords(3,1) = 0.9;  nodeCoords(3,2) = 0.1;
		nodeCoords(4,0) = 0.0; nodeCoords(4,1) = -0.1; nodeCoords(4,2) = 1.1;
		nodeCoords(5,0) = 1.1; nodeCoords(5,1) = 0.1;  nodeCoords(5,2) = 1.2;
		nodeCoords(6,0) = 0.9; nodeCoords(6,1) = 0.9;  nodeCoords(6,2) = 0.9;
		nodeCoords(7,0) = -0.1; nodeCoords(7,1) = 1.0;  nodeCoords(7,2) = 1.0;


		// set number of cells in cell workset

		numCells = (1<<i)+7;
		// copy generic hex coords into array 
		Intrepid::FieldContainer<double> cellWorkset(numCells, numNodesPerCell, spaceDim);
		for (int icell = 0; icell < numCells; icell++) {
			for (int inode = 0; inode < numNodesPerCell; inode++) { 
				for (int idim = 0; idim < spaceDim; idim++) { 
					cellWorkset(icell, inode, idim) = nodeCoords(inode,idim);
				}
			}
		}   


		/*********************************************************/
		/*               Define cell cubature                    */
		/*********************************************************/

		// Define cubature 
		int cubDegree = 2;
		std::vector< Teuchos::RCP< Intrepid::Cubature<double> > > lineCubs(3);
		lineCubs[0]  = Teuchos::rcp(new Intrepid::CubatureDirectLineGauss<double>(cubDegree));
		lineCubs[1]  = Teuchos::rcp(new Intrepid::CubatureDirectLineGauss<double>(cubDegree));
		lineCubs[2]  = Teuchos::rcp(new Intrepid::CubatureDirectLineGauss<double>(cubDegree));
		Intrepid::CubatureTensor<double> cellCubature(lineCubs);

		int cubDim       = cellCubature.getDimension(); 
		int numCubPoints = cellCubature.getNumPoints();
		// Get numerical integration points and weights
		Intrepid::FieldContainer<double> cubPoints (numCubPoints, cubDim);
		Intrepid::FieldContainer<double> cubWeights(numCubPoints);

		cellCubature.getCubature(cubPoints, cubWeights);

		/*********************************************************/
		/*                 Define cell basis                     */
		/*********************************************************/

		// Select basis
		Intrepid::Basis_HGRAD_HEX_C1_FEM<double, Intrepid::FieldContainer<double> > hexHGradBasis;
		int numFields = hexHGradBasis.getCardinality();

		Intrepid::FieldContainer<double> basisValues(numFields, numCubPoints);
		Intrepid::FieldContainer<double> basisGrads (numFields, numCubPoints, spaceDim);

		// Evaluate basis values and gradients at cubature points
		hexHGradBasis.getValues(basisValues, cubPoints, Intrepid::OPERATOR_VALUE);
		hexHGradBasis.getValues(basisGrads, cubPoints, Intrepid::OPERATOR_GRAD);

		/*********************************************************/
		/*            Calculate cell Jacobians                   */
		/*********************************************************/

		// Containers for Jacobians, inverses and determinants
		Intrepid::FieldContainer<double> worksetJacobian  (numCells, numCubPoints, spaceDim, spaceDim);
		Intrepid::FieldContainer<double> worksetJacobInv  (numCells, numCubPoints, spaceDim, spaceDim);
		Intrepid::FieldContainer<double> worksetJacobDet  (numCells, numCubPoints);

		gettimeofday(&start,NULL);

		Intrepid::CellTools<double>::setJacobian(worksetJacobian, cubPoints, cellWorkset, cellType);  

		Intrepid::CellTools<double>::setJacobianInv(worksetJacobInv, worksetJacobian );
		Intrepid::CellTools<double>::setJacobianDet(worksetJacobDet, worksetJacobian );


		/*********************************************************/
		/*             Compute Stiffness Matrix                  */
		/*********************************************************/



		// Containers for stiffness matrix
		Intrepid::FieldContainer<double> worksetStiffMatrix        (numCells, numFields, numFields);
		Intrepid::FieldContainer<double> worksetCubWeights         (numCells, numCubPoints);
		Intrepid::FieldContainer<double> worksetBasisGrads         (numCells, numFields, numCubPoints, spaceDim);
		Intrepid::FieldContainer<double> worksetBasisGradsWeighted (numCells, numFields, numCubPoints, spaceDim);

		// Transform basis gradients to physical frame:                        DF^{-T}(grad u)
		Intrepid::FunctionSpaceTools::HGRADtransformGRAD<double>(worksetBasisGrads,                
				                            worksetJacobInv,   basisGrads);
				                            
				                            	



		// Compute integration measure for workset cells:                      Det(DF)*w = J*w
		Intrepid::FunctionSpaceTools::computeCellMeasure<double>(worksetCubWeights,              
				                            worksetJacobDet, cubWeights);


		// Multiply transformed (workset) gradients with weighted measure:     DF^{-T}(grad u)*J*w
		Intrepid::FunctionSpaceTools::multiplyMeasure<double>(worksetBasisGradsWeighted,           
				                         worksetCubWeights, worksetBasisGrads);

		// Integrate to compute contribution to global stiffness matrix:      (DF^{-T}(grad u)*J*w)*(DF^{-T}(grad u))
		Intrepid::FunctionSpaceTools::integrate<double>(worksetStiffMatrix,                      
				                   worksetBasisGradsWeighted,
				                   worksetBasisGrads,Intrepid::COMP_CPP);
				    

		/**********************************************************************************/
		/*                                   Compute RHS                                  */
		/**********************************************************************************/

		// Containers for rhs
		Intrepid::FieldContainer<double> worksetRHS                (numCells, numFields);
		Intrepid::FieldContainer<double> worksetBasisValues        (numCells, numFields, numCubPoints);
		Intrepid::FieldContainer<double> worksetBasisValuesWeighted(numCells, numFields, numCubPoints);

		// Fill array with source term data values - for now fill with all ones
		Intrepid::FieldContainer<double> worksetSourceTerm(numCells, numCubPoints);
		for (int icell = 0; icell < numCells; icell++){
			for (int ipts = 0; ipts < numCubPoints; ipts++){
			 	worksetSourceTerm(icell,ipts) = 1.0;
			}
		}

		// Transform basis values to physical frame: clones basis values (u)
		Intrepid::FunctionSpaceTools::HGRADtransformVALUE<double>(worksetBasisValues,              
				                             basisValues);

		// Multiply transformed (workset) values with weighted measure:     (u)*J*w
		Intrepid::FunctionSpaceTools::multiplyMeasure<double>(worksetBasisValuesWeighted,        
				                         worksetCubWeights, worksetBasisValues);

		// Integrate worksetSourceTerm against weighted basis function set:  f.(u)*J*w
		Intrepid::FunctionSpaceTools::integrate<double>(worksetRHS,                             
				                   worksetSourceTerm,
				                   worksetBasisValuesWeighted,Intrepid::COMP_CPP);
				             
		gettimeofday(&stop, NULL);
		timersub(&stop , &start, &result);
		double time = (result.tv_sec + (result.tv_usec/1000000.0));

		std::cout	<<"Original kernel , "<<numCells<<" , "
			<< time
			<< " , " << (time/numCells)
			<<std::endl;                                

	}
}
