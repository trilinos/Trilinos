// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <string>
#include <fstream>

#include <boost/version.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  
  GlobalMPISession mpi_session(&argc, &argv);

  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    // *********************************************************************
    // Start of Graphviz Testing
    // *********************************************************************
#if (BOOST_VERSION>=104400)
    {
      
      typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
	boost::property<boost::vertex_name_t, std::string, 
	boost::property<boost::vertex_color_t, std::string> >,
	boost::property<boost::edge_name_t, std::string> > Graph;

      typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t;
      typedef boost::graph_traits<Graph>::edge_descriptor edge_t;
      
      Graph g;

      boost::add_edge(0,1,g);
      boost::add_edge(0,2,g);
      boost::add_edge(1,3,g);
      boost::add_edge(2,3,g);

      boost::dynamic_properties dp;
      dp.property("id",boost::get(boost::vertex_name, g));
      dp.property("fontcolor",boost::get(boost::vertex_color, g));
      
      boost::put("id",dp,(vertex_t) 0,std::string("EV0"));
      boost::put("id",dp,(vertex_t) 1,std::string("EV1"));
      boost::put("id",dp,(vertex_t) 2,std::string("EV2"));
      boost::put("id",dp,(vertex_t) 3,std::string("EV3"));

      boost::put("fontcolor",dp,(vertex_t) 0,std::string("black"));
      boost::put("fontcolor",dp,(vertex_t) 1,std::string("black"));
      boost::put("fontcolor",dp,(vertex_t) 2,std::string("red"));
      boost::put("fontcolor",dp,(vertex_t) 3,std::string("black"));

      const std::string filename = "Graphviz_file.dot";
      std::ofstream output_file(filename.c_str());
      write_graphviz_dp(std::cout, g, dp, std::string("id"));
      write_graphviz_dp(output_file, g, dp, std::string("id"));
    }
#endif

    // *********************************************************************
    // *********************************************************************
    std::cout << "\nTest passed!\n" << std::endl; 
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

  TimeMonitor::summarize();
    
  return 0;
}
