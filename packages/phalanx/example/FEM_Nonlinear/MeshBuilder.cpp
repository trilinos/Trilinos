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

#include <cmath>
#include "MeshBuilder.hpp"
#include "Element_Linear2D.hpp"
#include "Teuchos_TestForException.hpp"

//**********************************************************************
MeshBuilder::MeshBuilder(const Teuchos::RCP<Epetra_Comm>& comm,
			 int num_elements_x, int num_elements_y, 
			 double length_x, double length_y, 
			 int debug_level) :
  m_comm(comm),
  m_num_elements_x(num_elements_x),
  m_num_elements_y(num_elements_y),
  m_length_x(length_x),
  m_length_y(length_y),
  m_num_global_elements(m_num_elements_x * m_num_elements_y),
  m_debug_level(debug_level),
  m_print_process(m_comm->MyPID() == 0),
  m_my_elements(Teuchos::rcp(new std::vector<Element_Linear2D>(0)))
{
  
  int num_procs = m_comm->NumProc();
  int my_pid = m_comm->MyPID();
  int last_proc = num_procs - 1;

  // Divide elements across processors in x direction only.
  TEST_FOR_EXCEPTION(num_procs > num_elements_x,std::logic_error,
		     "Number of processors must be less than number of elements in the x direction.");

  // number of columns of elements on this processor
  int base_elements_per_process = num_elements_x / num_procs;
  int residual_elements = num_elements_x % num_procs;

  m_num_my_elements_x = base_elements_per_process;
  if (my_pid < residual_elements)
    m_num_my_elements_x += 1;

  // Determine the starting global id number
  int starting_global_element_id = 
    base_elements_per_process * my_pid * m_num_elements_y;
  
  starting_global_element_id += my_pid < residual_elements ?
    my_pid * m_num_elements_y : residual_elements * m_num_elements_y;
  
  int previous_num_elements_x = starting_global_element_id / m_num_elements_y;
  
  int starting_global_node_id = 
    previous_num_elements_x * (m_num_elements_y + 1);

  double dx = m_length_x / static_cast<double>(num_elements_x);
  double dy = m_length_y / static_cast<double>(num_elements_y);

  // Create elements
  std::vector<unsigned> global_node_ids(4);
  std::vector<double> x_coords(4);
  std::vector<double> y_coords(4);
  int global_element_id = starting_global_element_id;
  int local_element_id = 0;
  int my_offset_index_x = previous_num_elements_x;
    
  for (int xdir = 0; xdir < m_num_my_elements_x; ++ xdir) {

    for (int ydir = 0; ydir < m_num_elements_y; ++ ydir) {
    
      double x_left = static_cast<double>(my_offset_index_x + xdir) * dx;
      double x_right = x_left + dx;
      double y_bottom = static_cast<double>(ydir) * dy;
      double y_top = y_bottom + dy;

      x_coords[0] = x_left;
      x_coords[1] = x_right;
      x_coords[2] = x_right;
      x_coords[3] = x_left;
      y_coords[0] = y_bottom;
      y_coords[1] = y_bottom;
      y_coords[2] = y_top;
      y_coords[3] = y_top;

      global_node_ids[0] = 
	starting_global_node_id + (num_elements_y + 1) * xdir + ydir;
      global_node_ids[1] = global_node_ids[0] + num_elements_y + 1;
      global_node_ids[2] = global_node_ids[1] + 1;
      global_node_ids[3] = global_node_ids[0] + 1;
      
      Element_Linear2D e(global_node_ids, global_element_id, 
			 local_element_id, x_coords, y_coords);

      // Set the node ownership
      e.setOwnsNode(0, true);
      e.setOwnsNode(3, true);

      if ( (xdir == (m_num_my_elements_x - 1)) && (my_pid == last_proc) ) {
	e.setOwnsNode(1, true);
	e.setOwnsNode(2, true);
      }
      else if (xdir == (m_num_my_elements_x - 1)) {
	e.setOwnsNode(1, false);
	e.setOwnsNode(2, false);
      }
      else {
	e.setOwnsNode(1, true);
	e.setOwnsNode(2, true);
      }

      m_my_elements->push_back(e);

      ++global_element_id;
      ++local_element_id;
    }
  }
  
  // Create the boundary node lists
  m_left_node_set.clear();
  m_right_node_set.clear();
  m_top_node_set.clear();
  m_bottom_node_set.clear();
  int x_stop = 0;
  if (my_pid == last_proc) 
    x_stop = m_num_my_elements_x +1; // add the last column of nodes to proc
  else
    x_stop = m_num_my_elements_x;
    
  for (int xdir = 0; xdir < x_stop; ++xdir) {
    for (int ydir = 0; ydir < m_num_elements_y + 1; ++ydir) {
      int gid = starting_global_node_id + (m_num_elements_y + 1) * xdir + ydir;

      if (my_pid == 0 && xdir == 0)
	m_left_node_set.push_back(gid);
      
      if (my_pid == last_proc && xdir == (x_stop-1))
	m_right_node_set.push_back(gid);
      
      if (ydir == m_num_elements_y)
	m_top_node_set.push_back(gid);
      
      if (ydir == 0)
	m_bottom_node_set.push_back(gid);
    }
  }

}

//**********************************************************************
Teuchos::RCP< std::vector<Element_Linear2D> > 
MeshBuilder::myElements() const
{
  return m_my_elements;
}

//**********************************************************************
const std::vector<int>& MeshBuilder::leftNodeSetGlobalIds() const
{
  return m_left_node_set;
}

//**********************************************************************
const std::vector<int>& MeshBuilder::rightNodeSetGlobalIds() const
{
  return m_right_node_set;
}

//**********************************************************************
const std::vector<int>& MeshBuilder::topNodeSetGlobalIds() const
{
  return m_top_node_set;
}

//**********************************************************************
const std::vector<int>& MeshBuilder::bottomNodeSetGlobalIds() const
{
  return m_bottom_node_set;
}

//**********************************************************************
void MeshBuilder::print(std::ostream& os) const
{
  if (m_print_process)
    os << "MeshBuilder (debug_level: " << m_debug_level << ")" << std::endl;

  m_comm->Barrier();

  if (m_debug_level > 5)
    os  << "PID(" << m_comm->MyPID() << ") m_num_my_elements_x = " 
	<< m_num_my_elements_x << std::endl; 

  m_comm->Barrier();  
  
  if (m_debug_level > 6)
    for (std::size_t i = 0; i < m_my_elements->size(); ++i)
      os << (*m_my_elements)[i] << std::endl;

  m_comm->Barrier();  

  if (m_debug_level > 7)
    for (std::size_t i = 0; i < m_left_node_set.size(); ++i)
      os << "left node set (pid=" << m_comm->MyPID() << "):" 
	 << m_left_node_set[i] << std::endl;

  m_comm->Barrier(); 

  if (m_debug_level > 7)
    for (std::size_t i = 0; i < m_right_node_set.size(); ++i)
      os << "right node set (pid=" << m_comm->MyPID() << "):" 
	 << m_right_node_set[i] << std::endl;
 
  m_comm->Barrier(); 

  if (m_debug_level > 7)
    for (std::size_t i = 0; i < m_top_node_set.size(); ++i)
      os << "top node set (pid=" << m_comm->MyPID() << "):" 
	 << m_top_node_set[i] << std::endl;
 
  m_comm->Barrier(); 

  if (m_debug_level > 7)
    for (std::size_t i = 0; i < m_bottom_node_set.size(); ++i)
      os << "bottom node set (pid=" << m_comm->MyPID() << "):" 
	 << m_bottom_node_set[i] << std::endl;
 
}

//**********************************************************************
std::ostream& operator<<(std::ostream& os, const MeshBuilder& e)
{
  e.print(os);
  return os;
}

//**********************************************************************
