#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

#include <iostream>



#ifdef KOKKOS_ENABLE_CUDA
  typedef Kokkos::Cuda        ExecSpace ;
  typedef Kokkos::CudaSpace   MemSpace ;
#elif defined(KOKKOS_ENABLE_OPENMP)
  typedef Kokkos::OpenMP     ExecSpace ;
  typedef Kokkos::OpenMP     MemSpace ;
#else
  typedef Kokkos::Serial   ExecSpace ;
  typedef Kokkos::HostSpace   MemSpace ;
#endif


typedef Kokkos::Serial	sExec;
typedef Kokkos::HostSpace sMem;


//Kokkos::View<quad*, Kokkos::LayoutRight, MSpace>


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



//////////////////////////////////// Using the memory space set by compile flag

template <typename MSpace = MemSpace, typename ExSpace = ExecSpace>
struct baseGridConnector
{
	Kokkos::View<quad*, Kokkos::LayoutRight, MSpace> quads;
	Kokkos::View<node*, Kokkos::LayoutRight, MSpace> nodes;


	baseGridConnector(int rows, int columns,Kokkos::View<quad*, Kokkos::LayoutRight, MSpace> & _quads,Kokkos::View<node*, Kokkos::LayoutRight, MSpace> & _nodes )
	{	
		mRows=rows;
		mColumns=columns;
		nodes=_nodes;
		quads=_quads;
	}

	int mRows;
	int mColumns;

	KOKKOS_INLINE_FUNCTION void operator() (const int& i) const {
		quads(i).nodes[0]=4*i;
		quads(i).nodes[1]=4*i+1;
		quads(i).nodes[2]=4*i+2;
		quads(i).nodes[3]=4*i+3;}

	void connect()
	{
		Kokkos::parallel_for(Kokkos::RangePolicy<ExSpace>(0,mRows*mColumns),*this);
	}

}; //seems to run

template <typename MSpace = MemSpace, typename ExSpace = ExecSpace>
struct baseGridInitializer
{
	Kokkos::View<quad*, Kokkos::LayoutRight, MSpace> quads;
	Kokkos::View<node*, Kokkos::LayoutRight, MSpace> nodes;


	int mRows;
	int mColumns;

	double mWidth;
	double mHeight;
	
	Kokkos::View<double*, Kokkos::LayoutRight, MSpace> intervalsDevice;

	baseGridInitializer(int rows, int columns, Kokkos::View<quad*, Kokkos::LayoutRight, MSpace> & _quads, Kokkos::View<node*, Kokkos::LayoutRight, MSpace> & _nodes, 
		double width, double height){
		quads= _quads;
		nodes= _nodes;
		mRows=rows;
		mColumns=columns;
		mWidth=width;
		mHeight=height;

		double yInterval = height / (double) mRows; 
		double xInterval = width / (double) mColumns;
		
		Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> intervalsHost("intervalsHost",2);
		
		intervalsHost(0)=xInterval;
		intervalsHost(1)=yInterval;

		Kokkos::View<double*, Kokkos::LayoutRight, MSpace> interimIntDev("intervalsDevice",2);
		intervalsDevice = interimIntDev;

		Kokkos::deep_copy(intervalsDevice, intervalsHost);
	}

	

	KOKKOS_INLINE_FUNCTION void operator() (const int& i) const {
		for(int j=0; j<(mColumns);j++)
		{
			nodes(quads(i*mColumns + j).nodes[0]).x=j*intervalsDevice(0);
			nodes(quads(i*mColumns + j).nodes[0]).y=i*intervalsDevice(1);

			nodes(quads(i*mColumns + j).nodes[1]).x=(j+1)*intervalsDevice(0);
			nodes(quads(i*mColumns + j).nodes[1]).y=i*intervalsDevice(1);

			nodes(quads(i*mColumns + j).nodes[2]).x=(j+1)*intervalsDevice(0);
			nodes(quads(i*mColumns + j).nodes[2]).y=(i+1)*intervalsDevice(1);

			nodes(quads(i*mColumns + j).nodes[3]).x=j*intervalsDevice(0);
			nodes(quads(i*mColumns + j).nodes[3]).y=(i+1)*intervalsDevice(1);
		}
	}


	void initializeCoords()
	{
		Kokkos::parallel_for(Kokkos::RangePolicy<ExSpace>(0,mRows), *this);

	}

};


template <typename MSpace = MemSpace, typename ExSpace = ExecSpace>
struct refinedGridConnector
{
	Kokkos::View<quad*, Kokkos::LayoutRight, MSpace> quads;
	Kokkos::View<node*, Kokkos::LayoutRight, MSpace> nodes;

	int mRows;
	int mColumns;

	

	refinedGridConnector(int rows, int columns, Kokkos::View<quad*, Kokkos::LayoutRight, MSpace> & _quads, Kokkos::View<node*, Kokkos::LayoutRight, MSpace> & _nodes)
	{
		mRows = rows;
		mColumns = columns;
		quads=_quads;
		nodes=_nodes;
	}

	

	KOKKOS_INLINE_FUNCTION void operator() (const int& i) const {
		int iLocal = 4*i;
		int iNode = 9*i;
		quads(iLocal).nodes[0] = iNode;
		quads(iLocal).nodes[1] = iNode+1;
		quads(iLocal).nodes[2] = iNode+4;
		quads(iLocal).nodes[3] = iNode+3;

		quads(iLocal+1).nodes[0] = iNode+1;
		quads(iLocal+1).nodes[1] = iNode+2;
		quads(iLocal+1).nodes[2] = iNode+5;
		quads(iLocal+1).nodes[3] = iNode+4;

		quads(iLocal+2).nodes[0] = iNode+3;
		quads(iLocal+2).nodes[1] = iNode+4;
		quads(iLocal+2).nodes[2] = iNode+7;
		quads(iLocal+2).nodes[3] = iNode+6;

		quads(iLocal+3).nodes[0] = iNode+4;
		quads(iLocal+3).nodes[1] = iNode+5;
		quads(iLocal+3).nodes[2] = iNode+8;
		quads(iLocal+3).nodes[3] = iNode+7;
	}
	void connect()
	{
		Kokkos::parallel_for(Kokkos::RangePolicy<ExSpace>(0,mRows*mColumns), *this);
	}
};



template <typename MSpace = MemSpace, typename ExSpace = ExecSpace>
struct UniformRefiner
{
	Kokkos::View<quad*, Kokkos::LayoutRight, MSpace> quadViewDevice;
	Kokkos::View<node*, Kokkos::LayoutRight, MSpace> deviceNodes;

	//Kokkos::View<quad*, Kokkos::LayoutRight, MemSpace> refQuads;
	Kokkos::View<node*, Kokkos::LayoutRight, MSpace> refinedNodes;

	int mRows;
	int mColumns;

	double mWidth;
	double mHeight;
	
	

	UniformRefiner(int rows, int columns, double width, double height, 
	Kokkos::View<quad*, Kokkos::LayoutRight, MSpace> & _orgQuads,
	Kokkos::View<node*, Kokkos::LayoutRight, MSpace> & _orgNodes,
	/*Kokkos::View<quad*, Kokkos::LayoutRight, MemSpace> & _refQuads,*/
	Kokkos::View<node*, Kokkos::LayoutRight, MSpace> & _refNodes
	)
	{
		quadViewDevice = _orgQuads;
		deviceNodes = _orgNodes;

		//refQuads = _refQuads;
		refinedNodes = _refNodes;

		mRows=rows;
		mColumns=columns;

		mWidth=width;
		mHeight=height;
	}

	KOKKOS_INLINE_FUNCTION void operator() (const int& iQuad) const {
			int iNode = 9*iQuad; 

			refinedNodes(iNode).x = deviceNodes(quadViewDevice(iQuad).nodes[0]).x;
			refinedNodes(iNode).y = deviceNodes(quadViewDevice(iQuad).nodes[0]).y;



			refinedNodes(iNode+1).x = ( deviceNodes(quadViewDevice(iQuad).nodes[0]).x + deviceNodes(quadViewDevice(iQuad).nodes[1]).x )/2.0;
			refinedNodes(iNode+1).y = ( deviceNodes(quadViewDevice(iQuad).nodes[0]).y + deviceNodes(quadViewDevice(iQuad).nodes[1]).y )/2.0;

			refinedNodes(iNode+2).x = deviceNodes(quadViewDevice(iQuad).nodes[1]).x;
			refinedNodes(iNode+2).y = deviceNodes(quadViewDevice(iQuad).nodes[1]).y;


			refinedNodes(iNode+3).x = ( deviceNodes(quadViewDevice(iQuad).nodes[0]).x + deviceNodes(quadViewDevice(iQuad).nodes[3]).x )/2.0;
			refinedNodes(iNode+3).y = ( deviceNodes(quadViewDevice(iQuad).nodes[0]).y + deviceNodes(quadViewDevice(iQuad).nodes[3]).y )/2.0;

			refinedNodes(iNode+4).x = ( deviceNodes(quadViewDevice(iQuad).nodes[0]).x + deviceNodes(quadViewDevice(iQuad).nodes[1]).x
					+ deviceNodes(quadViewDevice(iQuad).nodes[2]).x + deviceNodes(quadViewDevice(iQuad).nodes[3]).x)/4.0;
			refinedNodes(iNode+4).y = ( deviceNodes(quadViewDevice(iQuad).nodes[0]).y + deviceNodes(quadViewDevice(iQuad).nodes[1]).y
							+ deviceNodes(quadViewDevice(iQuad).nodes[2]).y + deviceNodes(quadViewDevice(iQuad).nodes[3]).y)/4.0;

			refinedNodes(iNode+5).x = ( deviceNodes(quadViewDevice(iQuad).nodes[1]).x + deviceNodes(quadViewDevice(iQuad).nodes[2]).x )/2.0;
			refinedNodes(iNode+5).y = ( deviceNodes(quadViewDevice(iQuad).nodes[1]).y + deviceNodes(quadViewDevice(iQuad).nodes[2]).y )/2.0;


			refinedNodes(iNode+6).x = deviceNodes(quadViewDevice(iQuad).nodes[3]).x;
			refinedNodes(iNode+6).y = deviceNodes(quadViewDevice(iQuad).nodes[3]).y;

			refinedNodes(iNode+7).x = ( deviceNodes(quadViewDevice(iQuad).nodes[3]).x + deviceNodes(quadViewDevice(iQuad).nodes[2]).x )/2.0;
			refinedNodes(iNode+7).y = ( deviceNodes(quadViewDevice(iQuad).nodes[3]).y + deviceNodes(quadViewDevice(iQuad).nodes[2]).y )/2.0;

			refinedNodes(iNode+8).x = deviceNodes(quadViewDevice(iQuad).nodes[2]).x;
			refinedNodes(iNode+8).y = deviceNodes(quadViewDevice(iQuad).nodes[2]).y;
	
	}


	void refine()
	{
		Kokkos::parallel_for(Kokkos::RangePolicy<ExSpace>(0,mRows*mColumns), *this);
	}
};
////////////////////////////////////




//helper functions for printing output to verify correctness
void printBaseGridFromGPU(int rows, int columns, Kokkos::View<quad*, Kokkos::LayoutRight, MemSpace> & quadViewDevice, Kokkos::View<node*,
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

void printRefinedGridFromGPU(int rows, int columns, Kokkos::View<node*, Kokkos::LayoutRight, MemSpace> refinedViewNodes)
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
/////////////////////////////////////////////////////


TEST(refine, unstructuredFunctors)
{	
	#ifdef KOKKOS_ENABLE_CUDA
	std::cout << "Running Test with Kokkos::Cuda execution space and Kokkos::CudaSpace memory space" <<std::endl;

	#elif defined(KOKKOS_ENABLE_OPENMP)
	std::cout << "Running Test with Kokkos::OpenMP execution space and Kokkos::OpenMP memory space" <<std::endl;
	#else
	std::cout << "Running Test with Kokkos::Serial execution space and Kokkos::HostSpace memory space" <<std::endl;
	#endif

	std::cout << "Enter desired grid dimensions: " << std::endl;

	int x=1;
	int y=1;

	Kokkos::View<quad*, Kokkos::LayoutRight, MemSpace> baseQuadViewDevice("quadViewDevice",x*y);
	Kokkos::View<node*, Kokkos::LayoutRight, MemSpace> baseNodeViewDevice("nodeViewDevice",4*x*y);
	Kokkos::View<quad*, Kokkos::LayoutRight, MemSpace> refinedViewQuads("refinedGrid",4*y*x);
	Kokkos::View<node*, Kokkos::LayoutRight, MemSpace> refinedViewNodes("refinedNodes",9*y*x);

	baseGridConnector<> connect(x,y,baseQuadViewDevice,baseNodeViewDevice);

	connect.connect();

	baseGridInitializer<> baseInit(x,y,baseQuadViewDevice,baseNodeViewDevice, 1.0, 1.0);
	baseInit.initializeCoords();
	
	printBaseGridFromGPU(x,y,baseQuadViewDevice,baseNodeViewDevice);

	refinedGridConnector<> refInit(x,y,refinedViewQuads, refinedViewNodes);
	refInit.connect();
	
	UniformRefiner<> ref(x,y,1.0,1.0,baseQuadViewDevice, baseNodeViewDevice, refinedViewNodes);

	ref.refine();

	printRefinedGridFromGPU(x,y,refinedViewNodes);


	
	int multiplier = 128;
	int numRepeats = 5;

	

	
	
		std::cout << "Data for refinement on device" << std::endl;
		for(int i=1;i<=3;i++) //growing size of problem
		{

		double avgTime=0.0;
		double avgNPS=0.0;
		
		int dimOne = 2*i*multiplier;
		int dimTwo = i*multiplier;


		
		for(int j=0;j<numRepeats;j++){
			
			Kokkos::View<quad*, Kokkos::LayoutRight, MemSpace> baseQuads("baseQuads",dimOne*dimTwo);
			Kokkos::View<node*, Kokkos::LayoutRight, MemSpace> baseNodes("baseNodes",4*dimOne*dimTwo);

			Kokkos::View<quad*, Kokkos::LayoutRight, MemSpace> refinedQuads("refinedQuads",4*dimOne*dimTwo);
			Kokkos::View<node*, Kokkos::LayoutRight, MemSpace> refinedNodes("refinedNodes",9*dimOne*dimTwo);
	
			baseGridConnector<> baseConnector(dimOne,dimTwo,baseQuads,baseNodes);
			baseConnector.connect();

			baseGridInitializer<> baseInitializer(dimOne,dimTwo,baseQuads,baseNodes, 2.0, 2.0);
			baseInitializer.initializeCoords();
	

			struct timeval begin, end;

			gettimeofday(&begin,NULL);
			refinedGridConnector<> refGridInit(dimOne,dimTwo, refinedQuads, refinedNodes);
			refGridInit.connect();

			UniformRefiner<> refer(dimOne,dimTwo, 2.0,2.0, baseQuads, baseNodes, refinedNodes);
	
			refer.refine();
			gettimeofday(&end,NULL);

			double time = 1.0*(end.tv_sec-begin.tv_sec) +
     			 1.0e-6*(end.tv_usec-begin.tv_usec);
			double newNodesPerSec = (9.0*dimOne*dimTwo)/time;

			//std::cout << "Base time for a " << dimOne << " by " << dimTwo << " grid." << time <<" with a nodes/second of: " << newNodesPerSec << std::endl;
			avgTime +=time;
			avgNPS += newNodesPerSec;
			}
			
			std::cout << "    Average run time over "<<numRepeats<<" repeats for a " << dimOne << " by " << dimTwo << " grid." << avgTime/(double)numRepeats <<" with a nodes/second of: " << avgNPS/100.0 << std::endl;
		}

		
		std::cout << "Data for refinement in serial" << std::endl;
		for(int i=1;i<=5;i++)//growing size of problem
		{

		double avgTime=0.0;
		double avgNPS=0.0;
		
		int dimOne = 2*i*multiplier;
		int dimTwo = i*multiplier;

		for(int j=0;j<numRepeats;j++){
			
			Kokkos::View<quad*, Kokkos::LayoutRight, sMem> baseQuads("baseQuads",dimOne*dimTwo);
			Kokkos::View<node*, Kokkos::LayoutRight, sMem> baseNodes("baseNodes",4*dimOne*dimTwo);

			Kokkos::View<quad*, Kokkos::LayoutRight, sMem> refinedQuads("refinedQuads",4*dimOne*dimTwo);
			Kokkos::View<node*, Kokkos::LayoutRight, sMem> refinedNodes("refinedNodes",9*dimOne*dimTwo);
	
			baseGridConnector<sMem, sExec> baseConnector(dimOne,dimTwo,baseQuads,baseNodes);
			baseConnector.connect();

			baseGridInitializer<sMem, sExec> baseInitializer(dimOne,dimTwo,baseQuads,baseNodes, 2.0, 2.0);
			baseInitializer.initializeCoords();
	

			struct timeval begin, end;

			gettimeofday(&begin,NULL);
			refinedGridConnector<sMem, sExec> refGridInit(dimOne,dimTwo, refinedQuads, refinedNodes);
			refGridInit.connect();

			UniformRefiner<sMem, sExec> refer(dimOne,dimTwo, 2.0,2.0, baseQuads, baseNodes, refinedNodes);
	
			refer.refine();
			gettimeofday(&end,NULL);

			double time = 1.0*(end.tv_sec-begin.tv_sec) +
     			 1.0e-6*(end.tv_usec-begin.tv_usec);
			double newNodesPerSec = (9.0*dimOne*dimTwo)/time;

			//std::cout << "Base time for a " << dimOne << " by " << dimTwo << " grid." << time <<" with a nodes/second of: " << newNodesPerSec << std::endl;
			avgTime +=time;
			avgNPS += newNodesPerSec;
			}
			
			std::cout << "    Average run time over "<<numRepeats<<" repeats for a " << dimOne << " by " 
					<< dimTwo << " grid." << avgTime/(double)numRepeats <<" with a nodes/second of: " << avgNPS/(double)numRepeats << std::endl;
		}	
	
}
