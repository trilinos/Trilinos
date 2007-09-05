#include "Intrepid_MultiCell.hpp"

using namespace std;
using namespace Intrepid;


int main(int argc, char *argv[]) {

  cout << "TEST 1: class MultiCell in 2D";

  double trinodes[] = {0.0, 0.0,            // nodes of the first triangle
                       1.0, 0.0,
                       1.0, 1.0,
                       0.0, 0.0,            // nodes of the second triangle
                       1.0, 1.0,
                       0.0, 1.0,
                       0.0, 1.0,            // nodes of the third triangle
                       1.5, 1.0,
                       0.5, 2.0};

  short triedgeorients[] = {1, 1, -1,       // constructor does not check edge orientations
                           -1, -1, 1,       // for consistency! It is user's responsibility
                            1, 1, -1};      // to provide a single orientation per edge

  short trifaceorients[] = {0,              // triangles don't have faces, this is dummy data
                            0,              // because we don't have separate 2D constructor
                            0};

  MultiCell<double> triMcell(3,                   // number of cells (triangles) in the multicell instance
                             2,                   // ambient dimension
                             TRI,                 // type of cells forming the multicell
                             trinodes,            // list of node coordinates
                             triedgeorients,      // edge orientations
                             trifaceorients);     // face orientations (dummy)

  cout << triMcell << endl;         // display the newly created multicell

  cout << "Testing multicell interface for the generating cell type...\n\n";

  cout << "\t type                   = " << triMcell.getMyType() << "\n";
  cout << "\t name                   = " << triMcell.getMyName() << "\n";
  cout << "\t ambient dimension      = " << triMcell.getMyAmbientDimension() <<"\n";
  cout << "\t topological dimension  = " << triMcell.getMyTopoDimension() << "\n";
  cout << "\t # of nodes             = " << triMcell.getMyNumNodes() << "\n"; 
  cout << "\t # of 0-subcells        = " << triMcell.getMyNumSubcells(0) << "\n";
  cout << "\t # of 1-subcells        = " << triMcell.getMyNumSubcells(1) << "\n";
  cout << "\t # of 2-subcells        = " << triMcell.getMyNumSubcells(2) << "\n";
  cout << "\t # of 3-subcells        = " << triMcell.getMyNumSubcells(3) << "\n";

  cout << "\t 1-subcell with index 0 = " << triMcell.getName(triMcell.getMySubcellType(1,0)) <<"\n";
  cout << "\t 1-subcell with index 1 = " << triMcell.getName(triMcell.getMySubcellType(1,1)) <<"\n";
  cout << "\t 1-subcell with index 2 = " << triMcell.getName(triMcell.getMySubcellType(1,2)) <<"\n";

  cout << "\t 2-subcell with index 0 = " << triMcell.getName(triMcell.getMySubcellType(2,0)) <<"\n\n";


  std::vector<int> subcell_node_conn;               // space for the node connectivities
  triMcell.getMySubcellNodes(2,0,subcell_node_conn);
  cout << " Node connectivity of the 0th 1-subcell -> { ";
  for(unsigned int i=0; i<subcell_node_conn.size(); i++) cout << subcell_node_conn[i]<<" ";
  cout << "}\n\n";

  cout << "Accessing node coordinates of the cell with cell_id_ = 1...\n";
  Point<double> Node_0 = triMcell[1][0];            // triMcell[i][j] = triMcell.getPoint(i,j)
  Point<double> Node_1 = triMcell[1][1];            // and returns the j-th node of the ith cell
  Point<double> Node_2 = triMcell.getPoint(1,2);    // as a Point object

  cout << "triMcell[1][0] =\n" << Node_0 << "\n";   // << is overloaded for Point objects
  cout << "triMcell[1][1] =\n" << Node_1 << "\n";
  cout << "triMcell[1][2] =\n" << Node_2 << "\n";
  
  cout << " testing [] operator for Node_0 \n";
  cout << "Node_0[0] = "<< Node_0[0] <<"\n";
  cout << "Node_0[1] = "<< Node_0[1] <<"\n";
  
  cout << " testing [] operator for Node_1 \n";
  cout << "Node_1[0] = "<< Node_1[0] <<"\n";
  cout << "Node_1[1] = "<< Node_1[1] <<"\n";
  
  cout << "END TEST 1: class MultiCell in 2D\n\n";


  cout << "TEST 2: class MultiCell in 3D";

  double prismnodes[] = {0.0, 0.0, 0.0,         // nodes of the first prism
                         1.0, 0.0, 0.0,
                         0.5, 0.5, 0.0,
                         0.0, 0.0, 1.0,
                         1.0, 0.0, 1.0,
                         0.5, 0.5, 1.0, 
                         0.0, 0.0, 1.0,         // nodes of the second prism
                         1.0, 0.0, 1.0,         // = nodes of the first prism
                         0.5, 0.5, 1.0,         // offset by 1 in the z-axes
                         0.0, 0.0, 2.0,
                         1.0, 0.0, 2.0,
                         0.5, 0.5, 2.0};

  short prismedgeorients[] = {1, 1, 1, 1, 1, 1, 1, 1, 1,        // orientations of prism 1 edges
                             -1, -1, -1, 1, 1, 1, -1, -1, -1};  // orientations of prism 2 edges

  short prismfaceorients[] = {1, 1, 1, 1, 1,                    // orientations of prism 1 faces
                             -1, -1, -1, -1, -1};               // orientations of prism 2 faces

  MultiCell<double> prismMcell(2,                         // number of cells (prisms) in the multicell 
                               3,                         // ambient dimension
                               PRISM,                     // type of cells forming the multicell
                               prismnodes,                // list of node coordinates
                               prismedgeorients,          // list of edge orentations
                               prismfaceorients);         // list of face orientations

  cout << prismMcell << endl;

  cout << "Accessing node coordinates ...\n";
  cout << "prismMcell[1][4] =\n" << prismMcell[1][4] << "\n\n";

  cout << "END TEST 2: class MultiCell in 3D\n\n";
  
  return 0;
}
