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

#ifndef PHX_EXAMPLE_MESH_BUILDER_HPP
#define PHX_EXAMPLE_MESH_BUILDER_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_Comm.h"
#include "Element_Linear2D.hpp"

/** \brief Creates a mesh and allocates the linear algebra objects.

*/
class MeshBuilder {
  
public:

  MeshBuilder(const Teuchos::RCP<Epetra_Comm>& comm, int num_elements_x, 
	      int num_elements_y, double length_x, double length_y, 
	      int debug_level = 0);

  //! Returns the local elements owned/used by this processor
  Teuchos::RCP< std::vector<Element_Linear2D> > myElements() const;

  const std::vector<int>& leftNodeSetGlobalIds() const;

  const std::vector<int>& rightNodeSetGlobalIds() const;

  const std::vector<int>& topNodeSetGlobalIds() const;

  const std::vector<int>& bottomNodeSetGlobalIds() const;

  void print(std::ostream& os) const;

private:
  
  Teuchos::RCP<Epetra_Comm> m_comm;
  int m_num_elements_x;
  int m_num_elements_y;
  int m_length_x;
  int m_length_y;
  int m_num_global_elements;
  int m_debug_level;

  int m_num_my_elements_x;
  int m_num_my_elements;
  bool m_print_process;
  
  Teuchos::RCP< std::vector<Element_Linear2D> > m_my_elements;

  std::vector<int> m_left_node_set;
  std::vector<int> m_right_node_set;
  std::vector<int> m_top_node_set;
  std::vector<int> m_bottom_node_set;

};

std::ostream& operator<<(std::ostream& os, const MeshBuilder& b);

#endif
