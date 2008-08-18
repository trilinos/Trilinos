
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
      
      ArrayDimension<Cell,Node,Spatial,Spatial> n_dim;
      ArrayNatural<double*,Cell,Node,Spatial,Spatial> a;

      //cout << n_dim << endl;


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
