// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
