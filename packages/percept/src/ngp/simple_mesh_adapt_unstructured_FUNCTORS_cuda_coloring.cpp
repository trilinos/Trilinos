#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>


#include <iostream>
#include <cstdlib>
#include <sstream>
#include <Kokkos_Threads.hpp>

#ifdef KOKKOS_ENABLE_CUDA		
  typedef Kokkos::Cuda        ExecSpace ;
  typedef Kokkos::CudaSpace   MemSpace ;
  typedef Kokkos::LayoutRight DataLayout;
#elif defined(KOKKOS_ENABLE_OPENMP)
  typedef Kokkos::OpenMP     ExecSpace ;
  typedef Kokkos::OpenMP     MemSpace ;
  typedef Kokkos::LayoutRight DataLayout;
#else
  typedef Kokkos::Serial   ExecSpace ;
  typedef Kokkos::HostSpace   MemSpace ;
  typedef Kokkos::LayoutRight DataLayout;
#endif

typedef Kokkos::Serial	sExec;
typedef Kokkos::HostSpace sMem;

struct node {

	double x;
	double y;

	int color;
	
	node() {
		x = 0.0;
		y = 0.0;
		color =  -1;
	}//

	node(double _x, double _y) {
		x = _x;
		y = _y;
		color = -1; //initialize to invalid color
	}//

	~node() {
	}//
};

struct quad {

	int nodes[4]; //AN ARRAY OF INDICES 
	unsigned ID;
	int color; //color should be implicit in the view index
	bool isEmpty;


	quad() {
		isEmpty = true;
		ID = 0;
		color = -1; //initialize to invalid color NOTE: may not need, as element color will be implicit in the row index of the colored 2D view

	}

	~quad() {


	}
	//	  e0
	//	0 - - 1		 //assuming this ordering for edges. May not be able to do that in the future.
	//   e3 |     |	e1
	//	|     |
	//	3 - - 2
	//	  e2
	//
	//
	//back = nodes(0 to 1)
	//front = nodes(2 to 3)

	//left = (3 to 0)
	//right = (1 to 2)
};



//////////////////////////////////// Using the memory space set by compile flag


struct baseGridConnectorAndCoordInit
{	//predicated on nodes being shared between quads
	Kokkos::View<node*, DataLayout, MemSpace> nodes;
	Kokkos::View<quad*, DataLayout, MemSpace> quads;

	int mRows;
	int mColumns;
	int numNodes;
	int numQuads;

	bool connected;
	bool initialized;

	Kokkos::View<node*,DataLayout,MemSpace>::HostMirror nodeMirror;
	Kokkos::View<quad*, DataLayout, MemSpace>::HostMirror quadMirror;


	baseGridConnectorAndCoordInit(Kokkos::View<node*, DataLayout, MemSpace>  _nodes, Kokkos::View<quad*, DataLayout, MemSpace>  _quads,
			int _rows, int _cols)
	{
		nodes = _nodes;
		quads = _quads;
		
		mRows = _rows;
		mColumns=_cols;
		numNodes=(mRows+1)*(mColumns+1);
		numQuads = mRows*mColumns;

		connected =  false;
		initialized = false;

		nodeMirror =  Kokkos::create_mirror_view(nodes); //hard coding to the default typedefs for now because mirrors seem to not compile with templated data
		quadMirror =  Kokkos::create_mirror_view(quads);
	}

	~baseGridConnectorAndCoordInit(){}

	void connect(bool printDebug = true)
	{	
		//Kokkos::View<int*, DataLayout, MemSpace >::HostMirror theMirror; //kokkos bug with templating?
									//comment everything in thid function then uncomment the above. It gives a semicolon wanted complaint when compiling
									//If you switch from the template arguments to the typedefs defined above, it will compile
									//if it's just a view and not a mirror, it will compile with the template arguments

										
		//Kokkos::View<quad*, DataLayout, MemSpace>::HostMirror quadMirror; //hard coding to the default typedefs for now because mirrors seem to not compile with templated data
		//quadMirror = Kokkos::create_mirror_view(quads);
		Kokkos::deep_copy(quadMirror,quads);

		//connectivity setup in serial on host
		quad * prevHorizontal;
		quad * prevVertical;
		int nodesHorizontal = mColumns + 1;
		//int nodesVertical = mRows +1;
		for(int i=0;i<mRows;i++)
		{

			for(int j=0;j<mColumns;j++)
			{
				quad * currentEntityView = &(quadMirror(mColumns * i + j));
				quadMirror(mColumns * i + j).color = -1; //expicity coloring for now because of strange bug
				//if(printDebug)
				//	std::cout << "The color for quad " << mColumns * i + j << " is: " << quadMirror(mColumns * i + j).color <<std::endl;
				if(j==0)
					prevHorizontal = NULL;
				else
					prevHorizontal = &(quadMirror(mColumns * i + j - 1));

				if(i==0)
					prevVertical=NULL;
				else
					prevVertical =  &(quadMirror(mColumns * (i - 1) + j));



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

				currentEntityView->isEmpty = false;
				currentEntityView = NULL;
				prevHorizontal = NULL;
				prevVertical = NULL;
			}//endfor

		}//endfor
		//connnectivity done

		Kokkos::deep_copy(quads, quadMirror);
		connected =  true;
	}//end connect()
	
	void initCoords(double height, double width, bool printDebug=true)
	{
		if(!connected)
			return;
		double yInterval = height / (double) mRows; //how far apart each node will be in their given rows
		double xInterval = width / (double) mColumns;
		
		Kokkos::View<double*, DataLayout, sMem> intervalsHost("intervalsHost", 2);
		
		intervalsHost(0) = xInterval;
		intervalsHost(1) = yInterval;

		//Kokkos::View<node*,DataLayout,MemSpace>::HostMirror nodeMirror; 
		//nodeMirror = Kokkos::create_mirror_view(nodes);		
		Kokkos::deep_copy(nodeMirror,nodes);
		for(int i=0; i<(1 + mRows);i++)
		{
			for(int j=0; j<(1 + mColumns);j++)
			{
				nodeMirror(i*(1 + mColumns) + j).x = j*intervalsHost(0);
				nodeMirror(i*(1 + mColumns) + j).y = i*intervalsHost(1);
				nodeMirror(i*(1 + mColumns) + j).color = -1; //for now explicity coloring during initialization because colors seem to be getting overwritten with 0's 
										//despite constructor
			
				if(printDebug)
					std::cout<< "For node " << i*(1 + mColumns) + j << " the color is: " <<nodeMirror(i*(1 + mColumns) + j).color <<std::endl;
			}
		}

		Kokkos::deep_copy(nodes,nodeMirror);
		initialized = true;

	}//end initCoords()

	void coordDump() const
	{
		//Kokkos::View<node*,DataLayout,MemSpace>::HostMirror nodeMirror; //hard coding to the default typedefs for now because mirrors seem to not compile with templated data
		//nodeMirror = Kokkos::create_mirror_view(nodes);		
		Kokkos::deep_copy(nodeMirror,nodes);

		//Kokkos::View<quad*, DataLayout, MemSpace>::HostMirror quadMirror; //hard coding to the default typedefs for now because mirrors seem to not compile with templated data
		//quadMirror = Kokkos::create_mirror_view(quads);
		Kokkos::deep_copy(quadMirror,quads);
		
		std::cout << "Output from shared node view" << std::endl;
		for(int i = 0; i<mRows;i++ )
		{
			for(int j = 0; j< mColumns;j++)
			{	
				std::cout<<"This quad's color: "<< quadMirror(i*mColumns + j).color <<std::endl<<std::endl;

				std::cout << "X: " <<nodeMirror(quadMirror(i*mColumns + j).nodes[0]).x << " "
						<< "Y: "<<nodeMirror(quadMirror(i*mColumns + j).nodes[0]).y << " --- "
						<< "This node's color: " << nodeMirror(quadMirror(i*mColumns + j).nodes[0]).color << std::endl;

				std::cout << "X: " <<nodeMirror(quadMirror(i*mColumns + j).nodes[1]).x << " "
						<< "Y: "<<nodeMirror(quadMirror(i*mColumns + j).nodes[1]).y << " --- "
						<< "This node's color: " <<  nodeMirror(quadMirror(i*mColumns + j).nodes[1]).color << std::endl;

				std::cout << "X: " <<nodeMirror(quadMirror(i*mColumns + j).nodes[2]).x << " "
						<< "Y: "<<nodeMirror(quadMirror(i*mColumns + j).nodes[2]).y << " --- "
						<< "This node's color: " <<  nodeMirror(quadMirror(i*mColumns + j).nodes[2]).color << std::endl;
			
				std::cout << "X: " <<nodeMirror(quadMirror(i*mColumns + j).nodes[3]).x << " "
						<< "Y: "<<nodeMirror(quadMirror(i*mColumns + j).nodes[3]).y << " --- "
						<< "This node's color: " <<  nodeMirror(quadMirror(i*mColumns + j).nodes[3]).color << std::endl;

				std::cout<<"-----------------"<<std::endl;
			}
		}
	}//end coordDump()
	
}; //baseGridConnectorAndCoordInit


struct baseGridColorer
{
	//colors grid to maximize parallel performance
	//Takes a 1D grid of quads, sorts it into a 2-D grid of quads with each row of the 2-D grid corresponding to a thread color
	//This colors the datastructure but ultimately Kokkos decides the mapping of threads to work 

	//Kokkos::View<quad*,DataLayout, MemSpace> baseQuads;
	//Kokkos::View<node*,DataLayout, MemSpace> baseNodes;	//this just gets colored
	Kokkos::View<quad**,DataLayout, MemSpace> coloredQuads; //base gets sorted into this
	int baseQuadsSize;
	int sizeOfRefNodes; //as you're coloring the mesh, this number can be determined	
 
	Kokkos::View<quad**,DataLayout, MemSpace>::HostMirror hostColors;
	//Kokkos::View<quad*,DataLayout, MemSpace>::HostMirror hostQuads;
	//Kokkos::View<node*,DataLayout, MemSpace>::HostMirror hostNodes;

/*
	int num_threads;
	int num_colors;
	int cols;
	int leftOvers;
*/
	baseGridColorer(//Kokkos::View<quad*,DataLayout,MemSpace>  _baseQuads, 
		//Kokkos::View<node*,DataLayout, MemSpace>  _baseNodes, 
		int _baseQuadsSize, 
		int _baseNodesSize,
		//Kokkos::View<quad*,DataLayout, MemSpace>::HostMirror  _hostQuads,
		//Kokkos::View<node*,DataLayout, MemSpace>::HostMirror  _hostNodes,
		//Kokkos::View<quad**,DataLayout, MemSpace> _coloredQuads,
		bool printDebug = true)
	{
		//baseQuads =_baseQuads;
		//baseNodes =_baseNodes;
		baseQuadsSize =_baseQuadsSize;
		sizeOfRefNodes = _baseNodesSize; //will increase later
		
		//hostQuads = _hostQuads;//Kokkos::create_mirror_view(baseQuads);
		//hostNodes = _hostNodes;//Kokkos::create_mirror_view(baseNodes);

		//coloredQuads = _coloredQuads; //attempted to get rid of gllibc error by temporarily sizing coloredQuads before hand then resizing. Seems to like that less than my currentcode

		if(printDebug){
			std::cout << "Made it through the constructor of baseGridColorer" << std::endl;
			std::cout << "size of base quad grid is " << baseQuadsSize << std::endl;
			std::cout << "size of base node grid is " << sizeOfRefNodes << std::endl;
		}
	}

	~baseGridColorer(){}

	void color(Kokkos::View<node*,DataLayout, MemSpace>::HostMirror & hostNodes,
		Kokkos::View<quad*,DataLayout, MemSpace>::HostMirror & hostQuads,
		Kokkos::View<node*,DataLayout, MemSpace> & baseNodes,
		Kokkos::View<quad*,DataLayout, MemSpace> & baseQuads,
		bool printDebug=true)
	{			
#ifdef KOKKOS_ENABLE_CUDA
		int num_threads = 4992; //hard coding this to the number of threads available on the Telsa K80 GPU, the one found on the ASCIGPU machines
#elif defined(KOKKOS_ENABLE_OPENMP)
		const char* ompThrdCnt = std::getenv("OMP_NUM_THREADS");
		std::stringstream converter;
		converter << ompThrdCnt;
		int num_threads;
		converter >> num_threads;
#else
		int num_threads = 1;
#endif
		
		if(printDebug){
			std::cout << "Successfully determined the max number of possible threads that can be used as: " <<num_threads <<std::endl;
			std::cout<<"base quad size is: "<<baseQuadsSize <<std::endl; 
		}
		
		
		//calculate num colors needed
		int num_colors = 0;
		
		if(num_threads<baseQuadsSize)
			num_colors = num_threads;
		else
			num_colors=baseQuadsSize;
		
		
		if(printDebug)
			std::cout << "determined the number of colors (rows) needed to be: " <<num_colors <<std::endl;
		
		//calculate dimensions of colored view
		int cols = baseQuadsSize/num_colors;
		
		int leftOvers = baseQuadsSize%num_colors;
		
		if(leftOvers)
			cols++; //make room to append leftover elements
				//note this will mean our colored view may contain some "empty" elements, hence the field isEmpty on the quad DS
		
		if(printDebug)
			std::cout << "determined the number of columns in 2D view to be: " <<cols <<std::endl;
		
		
		coloredQuads = Kokkos::View<quad**,DataLayout, MemSpace>("colors",num_colors,cols); //doesn't mind this, it's the mirror creation in seems to have a problem with
		//Kokkos::Experimental::resize(coloredQuads, num_colors, cols);	//attempted to solve glibc error by having the coloredview allocated before then proceeding to resize it. Didn't work	
		hostColors = Kokkos::create_mirror_view(coloredQuads);		
		if(printDebug)
			std::cout << "Successfully allocated space for colored mesh view" <<std::endl;	
		//being the sorting process

		//Kokkos::View<quad**,DataLayout, MemSpace>::HostMirror hostColors; //it seems to like this approach even less
		//hostColors = Kokkos::create_mirror_view(coloredQuads); //doesn't like this mirror creation. Why? //perhaps it's getting reclaimed at the end of the function AND at the end of the program? I think there is a difference between creating this mirror in the constructor and creating within a function. I think it gets reclaimed within the scope of a function.
//That's why whenever I comment this function out in my main, my other code that creates mirrors in its constructor runs just finel
		Kokkos::deep_copy(hostColors, coloredQuads);
				
		Kokkos::deep_copy(hostQuads,baseQuads);

		Kokkos::deep_copy(hostNodes,baseNodes);
		
		if(printDebug)
			std::cout << "Successfully created mirrors" <<std::endl;
		
		int numFullRowsWritten = 0;
		for(int r=0; r < num_colors; r++)
		{
			if(printDebug)
				std::cout << "about to start coloring with color "<< r <<std::endl;

			sizeOfRefNodes++; //for each quad you color, you'll need a centroid for you refine nodes. This will change to 
						//sizeOfRefNodes += N^2 where N = number of nodes placed along an edge during refinement
			bool wFullRow = false;
			int iColsLim;
			if(r<leftOvers){
				iColsLim = cols;
				wFullRow = true;
				}
			else if(leftOvers)
				iColsLim = cols - -1;
			else //if not leftovers
				iColsLim = cols;
			
			for(int iCols = 0; iCols<iColsLim; iCols++)
			{
				int iBaseQuads = iCols+numFullRowsWritten*cols + (r-numFullRowsWritten)*iColsLim;
				hostQuads(iBaseQuads).color = r;
				hostColors(r,iCols) = hostQuads(iBaseQuads);
				//color nodes
				for(int iEdge = 0; iEdge<4;iEdge++)
				{
					
					int iNodeBegin = iEdge;
					int iNodeEnd = (iEdge < 3 ? iEdge+1: 0); //if you're edge 3, your ending node is not zero
					if (hostNodes(hostQuads(iBaseQuads).nodes[iNodeBegin]).color<0 || //if either node is not colored, you're coloring something on this edge
							hostNodes(hostQuads(iBaseQuads).nodes[iNodeEnd]).color<0){ //color less than zero are considered invalid/uncolored
						
						//you'll need to place a node along this edge later
						if( hostNodes(hostQuads(iBaseQuads).nodes[iNodeBegin]).color<0 ){
							hostNodes(hostQuads(iBaseQuads).nodes[iNodeBegin]).color = r;
							if(printDebug)
								std::cout << "coloring node " << hostQuads(iBaseQuads).nodes[iNodeBegin] << " color "<< r <<std::endl;
							}

						if( hostNodes(hostQuads(iBaseQuads).nodes[iNodeEnd]).color<0 ){
							hostNodes(hostQuads(iBaseQuads).nodes[iNodeEnd]).color =r;
							if(printDebug)
								std::cout << "coloring node " << hostQuads(iBaseQuads).nodes[iNodeEnd] << " color "<< r <<std::endl;
						}
					} //coloring
					if(hostNodes(hostQuads(iBaseQuads).nodes[iNodeBegin]).color==r ||
						hostNodes(hostQuads(iBaseQuads).nodes[iNodeEnd]).color==r){
						sizeOfRefNodes++;
					}//counting new nodes
					
					
				}//edges
			} //columns
					
			if(wFullRow)
				numFullRowsWritten++;
		}//colors
		
		//bring data back to device
		Kokkos::deep_copy(coloredQuads, hostColors);
		Kokkos::deep_copy(baseQuads,hostQuads);
		Kokkos::deep_copy(baseNodes,hostNodes);
		
		if (printDebug)
			std::cout << "the number of nodes to be allocated for our refine grid is: " << sizeOfRefNodes <<std::endl;
		
	}//color()

}; //baseGridColorer


void unshareNodesAfterColoring()
{ //Unshares nodes so threads don't create memory access conflicts when trying to perform calcuations
  //Not optimal memory wise, though, because even elements within the same spacial AND color domain will get their nodes unshared, which reall isn't necessary as only a singlee thread will ever need 
			//access to their field data


  
	
} //unshareNodesAfterColoring


////////////////////////////////////


void pf_serialized_test()
{
	for(int i=0;i<10;i++){
		Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(i,i+1),KOKKOS_LAMBDA (int j){
			//is it more efficient to do hostside things in a loop then copy them back over to the device or this?
			printf("%d\n",j);
		});}


}




int get_max_num_threads(int length)
{

	//Issue with trying to use team's: the max team size is 1024
	//Even if I could find the number of threads I have to work with on the device (gpu or cpu's) and dictate HOW MUCH work each team does with a range policy for such a team,
			//I'm not sure if I will get the expected mapping, which I need to conincide with how I color the mesh...How do I force color mapping?! 
	//Also, the max number of threads per team is 1024, what about the max number of teams? Is it just the max number of available threads?
		//Say I have a bunch of one thread teams that number the total number of threads, and I give each of them a team range, would that (physical) range be consistently distributed?
	Kokkos::View<int*, DataLayout, MemSpace > tempView("tempView",length);
	Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,length), KOKKOS_LAMBDA (int i) {
		//Kokkos::Impl::Threads bob;
			tempView(i) = 1; //Kokkos::league_size();//bob.league_rank() ;//OMP_GET_NUM_THREADS();
	});
	Kokkos::View<int*, DataLayout, MemSpace >::HostMirror theMirror;
	theMirror = Kokkos::create_mirror_view(tempView);
	Kokkos::deep_copy(theMirror,tempView);

	int sum =0;
	for(int i=0;i<length;i++)
		sum += theMirror(i);

	//if(const char* env_p = std::getenv("OMP_NUM_THREADS"))
        
	//I think this tells me the number of a available threads
	std::cout << "the sum is: "<<sum <<std::endl;
	return sum;//Kokkos::cuda_get_thread_num();

}

TEST(DISABLED_refine, unstructuredFunctorsColors)
{	
	#ifdef KOKKOS_ENABLE_CUDA
	std::cout << "Running Test with Kokkos::Cuda execution space and Kokkos::CudaSpace memory space" <<std::endl;
	#elif defined(KOKKOS_ENABLE_OPENMP)
	std::cout << "Running Test with Kokkos::OpenMP execution space and Kokkos::OpenMP memory space" <<std::endl;
	#else
	std::cout << "Running Test with Kokkos::Serial execution space and Kokkos::HostSpace memory space" <<std::endl;
	#endif
	/*(
	pf_serialized_test();

	
	node first(7.6,2.1);
	node second = first;
	std::cout << second.y << second.x <<std::endl;
	quad fQuad;
	fQuad.isEmpty = false;
	fQuad.color=7;
	quad sQuad;
	sQuad = fQuad;

	fQuad.isEmpty = true;
	fQuad.color=11;
	std::cout << sQuad.isEmpty <<  " " << sQuad.color << std::endl;l
	*/

	int rows = 2;
	int cols = 3;
	double width = 1.0;
	double height = 1.0;
	
	Kokkos::View<quad*,DataLayout,MemSpace> testQuads("testQuads",rows*cols);
	Kokkos::View<node*,DataLayout,MemSpace> testNodes("testNodes",(rows+1)*(cols+1));
	//Kokkos::View<quad**,DataLayout, MemSpace> coloredQuads ("colors",1,1); //this seemed to be causing glibc error...but only when I also run the connector code


	baseGridConnectorAndCoordInit init(testNodes,testQuads,rows, cols); // all thiis baseGridConnectorAndCoordInot code runs fine when I comment out coloredQuads from above!
	init.connect(false);
	init.initCoords(width,height,false);
	init.coordDump();
	
	baseGridColorer theColorer(//testQuads, //I get a glibc error when I run this, with out without the coloredQuads instantiation above
			//testNodes, 
			rows*cols,
			(rows+1)*(cols+1)//,
			//init.quadMirror,
			//init.nodeMirror//,
			//coloredQuads			
			);

	
	//init.quadMirror = 0;
	//init.nodeMirror = 0;
 	
	theColorer.color(init.nodeMirror, init.quadMirror, init.nodes,init.quads);
							
	init.coordDump(); //The above function appears to run correctly, I even make it to the end of the unit test. However, I'm still getting the glibc detected error at the end :/

	
	
	sleep(5); //I think memory is getting reclaimed twice as the error occurs after the sleep function

	

	
}//end test
