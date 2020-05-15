#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

#include <iostream>

#ifdef KOKKOS_ENABLE_CUDA
  typedef Kokkos::Cuda   ExecSpace ;
  typedef Kokkos::CudaSpace   MemSpace ;
#elif defined(KOKKOS_ENABLE_OPENMP)
  typedef Kokkos::OpenMP     ExecSpace ;
  typedef Kokkos::OpenMP     MemSpace ;
#else
  typedef Kokkos::Serial   ExecSpace ;
  typedef Kokkos::HostSpace   MemSpace ;
#endif

struct node {
	node() {
		x = 0.0;
		y = 0.0;
	}

	node(double _x, double _y) {
		x = _x;
		y = _y;
	}
	~node() {
	}

	double x;
	double y;

};

struct quad {


	int nodes[4]; //AN ARRAY OF INDICES TO OFFSET INTO THE NODE VIEW

	unsigned ID;

	quad * leftAdjElem; //used in fixup (removing redundant nodes)
	quad * rightAdjElem;
	quad * frontAdjElem;
	quad * backAdjElem;

	quad * topLeftChild; //used in refinement
	quad * topRightChild;
	quad * bottomLeftChild;
	quad * bottomRightChild;

	quad() {




		ID = 0;

		leftAdjElem = NULL; //a quad doesn't strictly own these so we won't delete the data they point to.
		rightAdjElem = NULL; //other classes should worry about that
		frontAdjElem = NULL;
		backAdjElem = NULL;

		topLeftChild = NULL;
		topRightChild = NULL;
		bottomLeftChild = NULL;
		bottomRightChild = NULL;
	}

	~quad() {


		leftAdjElem = NULL; //a quad doesn't strictly own these so we won't delete the data they point to.
		rightAdjElem = NULL; //other classes should worry about that
		frontAdjElem = NULL;
		backAdjElem = NULL;

		topLeftChild = NULL; //same as above
		topRightChild = NULL;
		bottomLeftChild = NULL;
		bottomRightChild = NULL;

	}
	//	  B
	//	0 - - 1
	//	|     |
	//	|     |
	//	3 - - 2
	//	  F
	//
	//
	//back = nodes(0 to 1)
	//front = nodes(2 to 3)

	//left = (3 to 0)
	//right = (1 to 2)
};







//BEGIN UnstructuredGrid///////////////////////////////////////////////////////////////////////////
//
//
//
class UnstructuredGrid {

//	n[0]	n[1]
//	     Q
//	n[3]	n[2]
//
//	GRID BASED NUMBERING FOR NODES ON HOST
//
//	0	1	2	3
//	 q1	 q2	 q3		....
//	4	5	6	7
//	 q4	 q5	 q6
//	8	9	10	11
//	 q7	 q8	  q9
//	12	13	14	15
//
//		.	  	    .
//		.	    	   	.
//		.	      	       .

public:

	int mRows;
	int mColumns;
	int numQuads;
	int numNodes;
	bool isInitialized;
	Kokkos::View<quad*, Kokkos::LayoutRight, Kokkos::HostSpace> quadView;
	Kokkos::View<node*, Kokkos::LayoutRight, Kokkos::HostSpace> hostNodes;

	Kokkos::View<quad*, Kokkos::LayoutRight, MemSpace> refinedGrid;
	Kokkos::View<node*, Kokkos::LayoutRight, MemSpace> refinedNodes;

	Kokkos::View<quad*, Kokkos::LayoutRight, MemSpace> quadViewDevice;
	Kokkos::View<node*, Kokkos::LayoutRight, MemSpace> deviceNodes;

	UnstructuredGrid(int rows, int columns, Kokkos::View<quad*, Kokkos::LayoutRight, Kokkos::HostSpace> & interimQuadViewHost, Kokkos::View<node*, Kokkos::LayoutRight, Kokkos::HostSpace> & interimNodeViewHost,
				Kokkos::View<quad*, Kokkos::LayoutRight, MemSpace> & interimQuadViewDevice, Kokkos::View<node*, Kokkos::LayoutRight, MemSpace> & interimNodeViewDevice,
				Kokkos::View<quad*, Kokkos::LayoutRight, MemSpace> & interimRefinedViewQuads, Kokkos::View<node*, Kokkos::LayoutRight, MemSpace> & interimRefinedViewNodes);
	void initialize(double height, double width);
	void constructAdjacency();
	void uniformRefine();
	void initializeQuads();


private:





};


UnstructuredGrid::UnstructuredGrid(int rows, int columns, Kokkos::View<quad*, Kokkos::LayoutRight, Kokkos::HostSpace> & interimQuadViewHost, Kokkos::View<node*, Kokkos::LayoutRight, Kokkos::HostSpace> & interimNodeViewHost,
				Kokkos::View<quad*, Kokkos::LayoutRight, MemSpace> & interimQuadViewDevice, Kokkos::View<node*, Kokkos::LayoutRight, MemSpace> & interimNodeViewDevice,
				Kokkos::View<quad*, Kokkos::LayoutRight, MemSpace> & interimRefinedViewQuads, Kokkos::View<node*, Kokkos::LayoutRight, MemSpace> & interimRefinedViewNodes) {
	//Kokkos::View<quad*, Kokkos::LayoutRight, Kokkos::HostSpace> interimQuadViewHost(
	//		"quadView", rows * columns); //perhaps there's something about the way I'm assigning them here that's off... But in my little test I assign the same fashion. Perhaps something goes out of scope
	quadView = interimQuadViewHost;																			//in this constuctor!?
	//Kokkos::View<node*, Kokkos::LayoutRight, Kokkos::HostSpace> interimNodeViewHost(
	//		"nodeView", (1 + rows) * (1 + columns));
	hostNodes = interimNodeViewHost;

	//Kokkos::View<quad*, Kokkos::LayoutRight, MemSpace> interimQuadViewDevice(
	//		"quadViewDevice", rows * columns);
	quadViewDevice = interimQuadViewDevice;
	//Kokkos::View<node*, Kokkos::LayoutRight, MemSpace> interimNodeViewDevice(
	//		"nodeViewDevice", 4 * rows * columns);
	deviceNodes = interimNodeViewDevice;					//4 distinct nodes per quad

	//Kokkos::View<quad*, Kokkos::LayoutRight, MemSpace> interimRefinedViewQuads(
	//		"refinedGrid", 4 * rows * columns); //each quad refines into four children
	refinedGrid = interimRefinedViewQuads;
	//Kokkos::View<node*, Kokkos::LayoutRight, MemSpace> interimRefinedViewNodes(
	//		"refinedNodes", 9 * rows * columns);
	refinedNodes = interimRefinedViewNodes;	//every nine of these represents the 4 children from
											//refinement of each quad
	mRows = rows;
	mColumns = columns;
	isInitialized = false;


}

//Some non-class helper functions while I try and figure out why PF's don't like being called within my class's function

//helper function. Should be eliminated once I get rid of lambdas
void PF_1_Connectivity(int rows, int columns, Kokkos::View<quad*, Kokkos::LayoutRight, MemSpace> & quadViewDevice)
{
	Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,rows*columns), KOKKOS_LAMBDA (int i) {
			quadViewDevice(i).nodes[0] = 4*i;
			quadViewDevice(i).nodes[1] = 4*i + 1;
			quadViewDevice(i).nodes[2] = 4*i + 2;
			quadViewDevice(i).nodes[3] = 4*i + 3;
			printf("finished connectivity for disjoint quads view on thread %i\n",i);
		});
}

//helper function. Should be eliminated once I get rid of lambdas
void PF_2_Connectivity(int rows, int columns, Kokkos::View<quad*, Kokkos::LayoutRight, MemSpace> & refinedGrid)
{
	Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,rows*columns), KOKKOS_LAMBDA (int i)//zero to eight, left to right, top to bottom
			{

			int iLocal = 4*i;
			int iNode = 9*i;
			refinedGrid(iLocal).nodes[0] = iNode;
			refinedGrid(iLocal).nodes[1] = iNode + 1;
			refinedGrid(iLocal).nodes[2] = iNode + 4;
			refinedGrid(iLocal).nodes[3] = iNode + 3;


			refinedGrid(iLocal+1).nodes[0] = iNode + 1;
			refinedGrid(iLocal+1).nodes[1] = iNode + 2;
			refinedGrid(iLocal+1).nodes[2] = iNode + 5;
			refinedGrid(iLocal+1).nodes[3] = iNode + 4;


			refinedGrid(iLocal+2).nodes[0] = iNode + 3;
			refinedGrid(iLocal+2).nodes[1] = iNode + 4;
			refinedGrid(iLocal+2).nodes[2] = iNode + 7;
			refinedGrid(iLocal+2).nodes[3] = iNode + 6;



			refinedGrid(iLocal+3).nodes[0] = iNode + 4;
			refinedGrid(iLocal+3).nodes[1] = iNode + 5;
			refinedGrid(iLocal+3).nodes[2] = iNode + 8;
			refinedGrid(iLocal+3).nodes[3] = iNode + 7;



	//		0	1	2
	//		 c0	 c1
	//		3	4	5
	//		 c2	 c3
	//		6	7	8
		});



}



void UnstructuredGrid::initializeQuads() {



	PF_1_Connectivity(mRows, mColumns, quadViewDevice);

	PF_2_Connectivity(mRows, mColumns, refinedGrid);





	//quad * currentEntityView;
	quad * prevHorizontal;
	quad * prevVertical;
	int nodesHorizontal = mColumns + 1;
	//int nodesVertical = mRows +1;
	for(int i=0;i<mRows;i++)
	{

		for(int j=0;j<mColumns;j++)
		{
			quad * currentEntityView = &(quadView(mColumns * i + j));

			if(j==0)
				prevHorizontal = NULL;
			else
				prevHorizontal = &(quadView(mColumns * i + j - 1));

			if(i==0)
				prevVertical=NULL;
			else
				prevVertical =  &(quadView(mColumns * (i - 1) + j));



			if (!prevVertical && !prevHorizontal)
			{
				(*currentEntityView).nodes[0] = 0;
				(*currentEntityView).nodes[1] = 1;
				(*currentEntityView).nodes[2] = nodesHorizontal+1;
				(*currentEntityView).nodes[3] = nodesHorizontal;
			}
			else if(!prevVertical && prevHorizontal)
			{
				(*currentEntityView).nodes[0] = (*prevHorizontal).nodes[1]; //these inherit
				(*currentEntityView).nodes[3] = (*prevHorizontal).nodes[2];

				(*currentEntityView).nodes[1] = (*prevHorizontal).nodes[1] + 1; //get their own offset
				(*currentEntityView).nodes[2] = (*prevHorizontal).nodes[2] + 1;
			}
			else if(prevVertical && !prevHorizontal)
			{
				(*currentEntityView).nodes[0] = (*prevVertical).nodes[3]; //these inherit
				(*currentEntityView).nodes[1] = (*prevVertical).nodes[2];

				(*currentEntityView).nodes[2] = (i+1)*nodesHorizontal+1; //get their own offset
				(*currentEntityView).nodes[3] = (i+1)*nodesHorizontal;
			}
			else // prevVertical && prevHorizontal
			{
				(*currentEntityView).nodes[0] = (*prevVertical).nodes[3]; //these inherit
				(*currentEntityView).nodes[1] = (*prevVertical).nodes[2];
				(*currentEntityView).nodes[3] = (*prevHorizontal).nodes[2];


				(*currentEntityView).nodes[2] = (i+1)*nodesHorizontal + (j+1); //gets its own offset
			}

		}

	}



}


//helper function. Should be eliminated once I get rid of lambdas
void PF_Coordinates(int rows, int columns,Kokkos::View<quad*, Kokkos::LayoutRight, MemSpace> quadViewDevice,
			Kokkos::View<node*, Kokkos::LayoutRight, MemSpace> deviceNodes, Kokkos::View<double*, Kokkos::LayoutRight, MemSpace> intervals )
{

//	Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,rows*columns), KOKKOS_LAMBDA (int i)
//	interval(0) = xInterval
//	interval(1) = yInterval
	Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,rows), KOKKOS_LAMBDA (int i)
		{
			for(int j=0; j<(columns);j++)
			{

				deviceNodes(quadViewDevice(i*columns + j).nodes[0]).x=j*intervals(0);
				deviceNodes(quadViewDevice(i*columns + j).nodes[0]).y=i*intervals(1);

				deviceNodes(quadViewDevice(i*columns + j).nodes[1]).x=(j+1)*intervals(0);
				deviceNodes(quadViewDevice(i*columns + j).nodes[1]).y=i*intervals(1);

				deviceNodes(quadViewDevice(i*columns + j).nodes[2]).x=(j+1)*intervals(0);
				deviceNodes(quadViewDevice(i*columns + j).nodes[2]).y=(i+1)*intervals(1);

				deviceNodes(quadViewDevice(i*columns + j).nodes[3]).x=j*intervals(0);
				deviceNodes(quadViewDevice(i*columns + j).nodes[3]).y=(i+1)*intervals(1);


			}
		});
}

void UnstructuredGrid::initialize(double height, double width) {
	numQuads = mRows * mColumns;
	numNodes = (1 + mRows) * (1 + mColumns);

	std::cout << "About to call function to setup connectivity from within initialization function" << std::endl;
	initializeQuads();
	std::cout << "Connectivity Set" << std::endl;

	double yInterval = height / (double) mRows; //how far apart each node will be in their given rows
	double xInterval = width / (double) mColumns;

	Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> intervalsHost("intervalsHost", 2);
	Kokkos::View<double*, Kokkos::LayoutRight, MemSpace> intervalsDevice("intervalsDevice", 2);
	intervalsHost(0) = xInterval;
	intervalsHost(1) = yInterval;

	for(int i=0; i<(1 + mRows);i++)
	{
		for(int j=0; j<(1 + mColumns);j++)
		{
			hostNodes(i*(1 + mColumns) + j).x = j*intervalsHost(0);
			hostNodes(i*(1 + mColumns) + j).y = i*intervalsHost(1);
		}
	}

	Kokkos::deep_copy(intervalsDevice, intervalsHost);

	PF_Coordinates(mRows, mColumns, quadViewDevice, deviceNodes, intervalsDevice);
	std::cout << "Finished assigning coordinates to initial grids" <<std::endl;


	isInitialized = true;

}

void UnstructuredGrid::constructAdjacency() {
	//adj matrix will be used to expedite fixing up of shared nodes (making unshared nodes shared again)

	//for each elem
	//for each adj elem
	//if nodes' field data match, make them shared
}

void PF_UMR //helper function. Should be eliminated once I get rid of lambdas
(
		int numQuads,
		Kokkos::View<node*, Kokkos::LayoutRight, MemSpace> refinedNodes,
		Kokkos::View<quad*, Kokkos::LayoutRight, MemSpace> quadViewDevice,
		Kokkos::View<node*, Kokkos::LayoutRight, MemSpace> deviceNodes
)
{
	std::cout << "About to enter PF ... " ;
	Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,numQuads), KOKKOS_LAMBDA (const int iQuad/*, int iNode, double & xCoord, double & yCoord*/ /*tried adding some variables up here for assignment*/){
														// doesn't like this...

			printf("Starting interation %i\n", iQuad ); //errors after reaching here, each thread is reaching
			int iNode = 9*iQuad; //how costly is stuff like this? should I try and do this earlier?
			//gather then scatter to make more readable

			printf("Set local node index to %i, about to start calculating coordinates on thread %i\n", iNode, iQuad );

			refinedNodes(iNode).x = deviceNodes(quadViewDevice(iQuad).nodes[0]).x;
			refinedNodes(iNode).y = deviceNodes(quadViewDevice(iQuad).nodes[0]).y;

			printf("Calculated node %i on thread  %i\n", iNode, iQuad );



			refinedNodes(iNode+1).x = ( deviceNodes(quadViewDevice(iQuad).nodes[0]).x + deviceNodes(quadViewDevice(iQuad).nodes[1]).x )/2.0;
			refinedNodes(iNode+1).y = ( deviceNodes(quadViewDevice(iQuad).nodes[0]).y + deviceNodes(quadViewDevice(iQuad).nodes[1]).y )/2.0;
			printf("Calculated node %i on thread  %i\n", iNode+1, iQuad );

			refinedNodes(iNode+2).x = deviceNodes(quadViewDevice(iQuad).nodes[1]).x;
			refinedNodes(iNode+2).y = deviceNodes(quadViewDevice(iQuad).nodes[1]).y;
			printf("Calculated node %i on thread  %i\n", iNode+2, iQuad );


			refinedNodes(iNode+3).x = ( deviceNodes(quadViewDevice(iQuad).nodes[0]).x + deviceNodes(quadViewDevice(iQuad).nodes[3]).x )/2.0;
			refinedNodes(iNode+3).y = ( deviceNodes(quadViewDevice(iQuad).nodes[0]).y + deviceNodes(quadViewDevice(iQuad).nodes[3]).y )/2.0;
			printf("Calculated node %i on thread  %i\n", iNode+3, iQuad );

			refinedNodes(iNode+4).x = ( deviceNodes(quadViewDevice(iQuad).nodes[0]).x + deviceNodes(quadViewDevice(iQuad).nodes[1]).x
					+ deviceNodes(quadViewDevice(iQuad).nodes[2]).x + deviceNodes(quadViewDevice(iQuad).nodes[3]).x)/4.0;
			refinedNodes(iNode+4).y = ( deviceNodes(quadViewDevice(iQuad).nodes[0]).y + deviceNodes(quadViewDevice(iQuad).nodes[1]).y
							+ deviceNodes(quadViewDevice(iQuad).nodes[2]).y + deviceNodes(quadViewDevice(iQuad).nodes[3]).y)/4.0;
			printf("Calculated node %i on thread  %i\n", iNode+4, iQuad );

			refinedNodes(iNode+5).x = ( deviceNodes(quadViewDevice(iQuad).nodes[1]).x + deviceNodes(quadViewDevice(iQuad).nodes[2]).x )/2.0;
			refinedNodes(iNode+5).y = ( deviceNodes(quadViewDevice(iQuad).nodes[1]).y + deviceNodes(quadViewDevice(iQuad).nodes[2]).y )/2.0;
			printf("Calculated node %i on thread  %i\n", iNode+5, iQuad );


			refinedNodes(iNode+6).x = deviceNodes(quadViewDevice(iQuad).nodes[3]).x;
			refinedNodes(iNode+6).y = deviceNodes(quadViewDevice(iQuad).nodes[3]).y;
			printf("Calculated node %i on thread  %i\n", iNode+6, iQuad );

			refinedNodes(iNode+7).x = ( deviceNodes(quadViewDevice(iQuad).nodes[3]).x + deviceNodes(quadViewDevice(iQuad).nodes[2]).x )/2.0;
			refinedNodes(iNode+7).y = ( deviceNodes(quadViewDevice(iQuad).nodes[3]).y + deviceNodes(quadViewDevice(iQuad).nodes[2]).y )/2.0;
			printf("Calculated node %i on thread  %i\n", iNode+7, iQuad );

			refinedNodes(iNode+8).x = deviceNodes(quadViewDevice(iQuad).nodes[2]).x;
			refinedNodes(iNode+8).y = deviceNodes(quadViewDevice(iQuad).nodes[2]).y;
			printf("Calculated node %i on thread  %i\n", iNode+8, iQuad );

			printf("All shmandy dandy and done on iteration %i\n", iQuad );//errors before reaching here. No thread seems to reach here.

		});
		std::cout << " Done with PF. Refinement Function ending, returning to main. " << std::endl;



}

void UnstructuredGrid::uniformRefine() {

	if (!isInitialized)
		return; //can't refine if you have no mesh

	//refinedNodes(0).y = 5.0; //this works just fine
	//refinedNodes(35).y = 5.0;
	//refinedNodes(100000).y = 5.0; //as does this, which is weird because my input was a 2 by 2 meaning the largest this view should be is 36 (0 to 35)
	PF_UMR(numQuads, refinedNodes, quadViewDevice, deviceNodes);



}
//
//
//
//End UnstructuredGrid///////////////////////////////////////////////////////////////////////////




void printBaseGridFromDevice(int rows, int columns, Kokkos::View<quad*, Kokkos::LayoutRight, MemSpace> & quadViewDevice, Kokkos::View<node*,
		Kokkos::LayoutRight, MemSpace> & deviceNodes){

		std::cout << "Output from unshared node view" <<std::endl;

		Kokkos::View<quad*, Kokkos::LayoutRight, Kokkos::HostSpace> unsharedHostCopy("unsharedCopy",rows*columns);
		Kokkos::View<node*, Kokkos::LayoutRight, Kokkos::HostSpace> theNodes("nodes",4*rows*columns);
		Kokkos::deep_copy(unsharedHostCopy,quadViewDevice);
		Kokkos::deep_copy(theNodes,deviceNodes);


		for(int i = 0; i<rows;i++ )
		{

			for(int j = 0; j< columns;j++)
			{
				std::cout << "X: " <<theNodes(unsharedHostCopy(i*columns + j).nodes[0]).x << " "
						<< "Y: "<<theNodes(unsharedHostCopy(i*columns + j).nodes[0]).y << std::endl;

				std::cout << "X: " <<theNodes(unsharedHostCopy(i*columns + j).nodes[1]).x << " "
									<< "Y: "<<theNodes(unsharedHostCopy(i*columns + j).nodes[1]).y << std::endl;

				std::cout << "X: " <<theNodes(unsharedHostCopy(i*columns + j).nodes[2]).x << " "
									<< "Y: "<<theNodes(unsharedHostCopy(i*columns + j).nodes[2]).y << std::endl;

				std::cout << "X: " <<theNodes(unsharedHostCopy(i*columns + j).nodes[3]).x << " "
									<< "Y: "<<theNodes(unsharedHostCopy(i*columns + j).nodes[3]).y << std::endl;
				std::cout<<std::endl;
			}

		}

}

void printRefinedGridFromDevice(int rows, int columns, Kokkos::View<node*, Kokkos::LayoutRight, MemSpace> refinedViewNodes)
{
	Kokkos::View<node*, Kokkos::LayoutRight, Kokkos::HostSpace> refNodesHost("refNodesHost", 9*rows*columns);
	Kokkos::deep_copy(refNodesHost, refinedViewNodes);

	for(int i=0;i<rows*columns;i++)
	{
		int j=9*i;
		std::cout << "from element " << i << " in oldQuads" <<std::endl;
		std::cout << "(" <<refNodesHost(j).x <<", " <<refNodesHost(j).y << ")";
		std::cout << " (" <<refNodesHost(j+1).x <<", " <<refNodesHost(j+1).y << ")";
		std::cout << " (" <<refNodesHost(j+2).x <<", " <<refNodesHost(j+2).y << ")";
		std:: cout << std::endl;
		std:: cout << std::endl;

		std::cout << "from element " << i << " in oldQuads" <<std::endl;
		std::cout << "("  <<refNodesHost(j+3).x <<", " <<refNodesHost(j+3).y << ")";
		std::cout << " (" <<refNodesHost(j+4).x <<", " <<refNodesHost(j+4).y << ")";
		std::cout << " (" <<refNodesHost(j+5).x <<", " <<refNodesHost(j+5).y << ")";
		std:: cout << std::endl;
		std:: cout << std::endl;

		std::cout << "from element " << i << " in oldQuads" <<std::endl;
		std::cout << "("  <<refNodesHost(j+6).x <<", " <<refNodesHost(j+6).y << ")";
		std::cout << " (" <<refNodesHost(j+7).x <<", " <<refNodesHost(j+7).y << ")";
		std::cout << " (" <<refNodesHost(j+8).x <<", " <<refNodesHost(j+8).y << ")";
		std:: cout << std::endl;
		std:: cout << std::endl;
		std::cout <<"------------------------------------------------------------" << std::endl;
	}

}
/*
TEST(refine, unstructured)
{

#ifdef KOKKOS_ENABLE_CUDA
	std::cout << "Running Test with Kokkos::Cuda execution space and Kokkos::CudaSpace memory space" <<std::endl;

#elif defined(KOKKOS_ENABLE_OPENMP)
	std::cout << "Running Test with Kokkos::OpenMP execution space and Kokkos::OpenMP memory space" <<std::endl;
#else
	std::cout << "Running Test with Kokkos::Serial execution space and Kokkos::HostSpace memory space" <<std::endl;
#endif

	//std::cout << "Enter desired grid dimensions: " << std::endl;
	int x = 4;
	int y = 4;

	//std::cin>>x;
	//std::cin>>y;

	Kokkos::View<quad*, Kokkos::LayoutRight, Kokkos::HostSpace> interimQuadViewHost("quadViewHost",x*y);
	Kokkos::View<node*, Kokkos::LayoutRight, Kokkos::HostSpace> interimNodeViewHost("nodeViewHost",(x+1)*(y+1));
	Kokkos::View<quad*, Kokkos::LayoutRight, MemSpace> interimQuadViewDevice("quadViewDevice",x*y);
	Kokkos::View<node*, Kokkos::LayoutRight, MemSpace> interimNodeViewDevice("nodeViewDevice",4*x*y);
	Kokkos::View<quad*, Kokkos::LayoutRight, MemSpace> interimRefinedViewQuads("refinedGrid",4*y*x);
	Kokkos::View<node*, Kokkos::LayoutRight, MemSpace> interimRefinedViewNodes("refinedNodes",9*y*x);

	UnstructuredGrid testGrid(x,y, interimQuadViewHost, interimNodeViewHost,
			interimQuadViewDevice, interimNodeViewDevice, interimRefinedViewQuads, interimRefinedViewNodes);
	std::cout << "Passed constructor" << std::endl;
	testGrid.initialize(1.0,1.0);




	std::cout << "Output from shared node view" << std::endl;
	for(int i = 0; i<testGrid.mRows;i++ )
	{
		for(int j = 0; j< testGrid.mColumns;j++)
		{
			std::cout << "X: " <<testGrid.hostNodes(testGrid.quadView(i*testGrid.mColumns + j).nodes[0]).x << " "
					<< "Y: "<<testGrid.hostNodes(testGrid.quadView(i*testGrid.mColumns + j).nodes[0]).y << std::endl;
			std::cout << "X: " <<testGrid.hostNodes(testGrid.quadView(i*testGrid.mColumns + j).nodes[1]).x << " "
								<< "Y: "<<testGrid.hostNodes(testGrid.quadView(i*testGrid.mColumns + j).nodes[1]).y << std::endl;
			std::cout << "X: " <<testGrid.hostNodes(testGrid.quadView(i*testGrid.mColumns + j).nodes[2]).x << " "
								<< "Y: "<<testGrid.hostNodes(testGrid.quadView(i*testGrid.mColumns + j).nodes[2]).y << std::endl;
			std::cout << "X: " <<testGrid.hostNodes(testGrid.quadView(i*testGrid.mColumns + j).nodes[3]).x << " "
								<< "Y: "<<testGrid.hostNodes(testGrid.quadView(i*testGrid.mColumns + j).nodes[3]).y << std::endl;
			std::cout<<std::endl;
		}

	}

	printBaseGridFromDevice(testGrid.mRows, testGrid.mColumns, interimQuadViewDevice, interimNodeViewDevice );

	testGrid.uniformRefine();

	printRefinedGridFromDevice(testGrid.mRows, testGrid.mColumns, interimRefinedViewNodes);

}
*/

