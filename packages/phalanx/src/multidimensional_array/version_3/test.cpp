
#include <iostream>
#include "Array.hpp"

struct Spatial : public phdmesh::ArrayDimTag {
  const char * name() const ;
  static const Spatial& descriptor();
};

struct Quadrature : public phdmesh::ArrayDimTag {
  const char * name() const ;
  static const Quadrature& descriptor();
};

struct Node : public phdmesh::ArrayDimTag {
  const char * name() const ;
  static const Node& descriptor();
};

struct Cell : public phdmesh::ArrayDimTag {
  const char * name() const ;
  static const Cell& descriptor();
};

const char * Spatial::name() const 
{ static const char n[] = "Spatial" ; return n ; }
const Spatial & Spatial::descriptor() 
{ static const Spatial myself ; return myself ; }

const char * Quadrature::name() const 
{ static const char n[] = "Quadrature" ; return n ; }
const Quadrature & Quadrature::descriptor() 
{ static const Quadrature myself ; return myself ; }

const char * Node::name() const 
{ static const char n[] = "Node" ; return n ; }
const Node & Node::descriptor() 
{ static const Node myself ; return myself ; }

const char * Cell::name() const 
{ static const char n[] = "Cell" ; return n ; }
const Cell & Cell::descriptor() 
{ static const Cell myself ; return myself ; }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace phdmesh;
  
  try {
    
    {
      cout << "\nStarting MultiDimensionalArray Example!\n" << endl;

      std::vector<double> memory(72);
      Array<double,NaturalOrder,Cell,Node,Spatial,Spatial> 
	a(&memory[0], 3, 4, 2, 2);

      std::vector<Array<double,NaturalOrder,Cell,Node,Spatial,Spatial>::size_type> dimensions;

      std::size_t cell_dim = a.dimension(0);
      std::size_t node_dim = a.dimension(1);
      std::size_t row_dim = a.dimension(2);
      std::size_t col_dim = a.dimension(3);

      for (std::size_t cell = 0; cell < cell_dim; ++cell)
	for (std::size_t node = 0; node < node_dim; ++node)
	  for (std::size_t row = 0; row < row_dim; ++row)
	    for (std::size_t col = 0; col < col_dim; ++col)
	      {
		a(cell,node,row,col) = 2.0;
		cout << "a[" << cell <<"," << node << "," 
		     << row << "," << col << "] = " 
		     << a(cell,node,row,col) << endl; 
	      }
      
      cout << endl;

      for (std::size_t i = 0; i < a.size(); ++i)
	{
	  a[i] = 3.0;
	  cout << "a[" << i << "] = " << a[i] << endl;
	}

      // check truncating down to matrix
      Array<double,NaturalOrder,Spatial,Spatial> m = (a.truncate(0)).truncate(0);


      cout << "\nFinished MultiDimensionalArray Example!\n" << endl;
    }

    // *********************************************************************
    // Finished all testing
    // *********************************************************************
    std::cout << "\nRun has completed successfully!\n" << std::endl; 
    // *********************************************************************
    // *********************************************************************
    
  }
  catch (const std::exception& e) {
    std::cout << "************************************************" << endl;
    std::cout << "************************************************" << endl;
    std::cout << "Exception Caught!" << endl;
    std::cout << "Error message is below\n " << e.what() << endl;
    std::cout << "************************************************" << endl;
  }
  catch (...) {
    std::cout << "************************************************" << endl;
    std::cout << "************************************************" << endl;
    std::cout << "Unknown Exception Caught!" << endl;
    std::cout << "************************************************" << endl;
  }

  return 0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
