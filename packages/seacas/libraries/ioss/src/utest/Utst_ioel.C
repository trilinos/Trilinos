// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
//         
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <Ioss_CodeTypes.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#define OUTPUT std::cerr

#include <Ioss_ConcreteVariableType.h>
#include <Ioss_Initializer.h>
#include <Ioss_VariableType.h>
#include <Ioss_Utils.h>

#include <Ioss_ElementTopology.h>


using namespace Ioss;

// ========================================================================
static int  test_all_elements();
static void test_aliases(const NameList &elements);
static bool test_element(const std::string& type);
// ========================================================================

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  StorageInitializer initialize_storage;
  Ioss::Initializer  initialize_topologies;

  int err_count = test_all_elements();
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  OUTPUT << "\n" << argv[0];;
  if (err_count == 0) {
    OUTPUT << "\nSIERRA execution successful." << std::endl;
    return EXIT_SUCCESS;
  } else {
    OUTPUT << "\nSIERRA execution failed." << std::endl;
    return EXIT_FAILURE;
  }
}

int test_all_elements()
{
  int err_count = 0;

  NameList elements;
  int element_count = Ioss::ElementTopology::describe(&elements);

  OUTPUT.setf(std::ios::left);
  for (int i=0; i < element_count; i++) {
    OUTPUT << "Testing element: " << std::setw(20) << elements[i];
    bool result = test_element(elements[i]);
    if (result == true || elements[i] == "unknown" || elements[i] == "invalid_topology")
      OUTPUT << "OK" << '\n';
    else {
      OUTPUT << "\n        element: " << std::setw(20) << elements[i]
	   << "FAIL" << '\n';
      err_count++;
    }
  }

  test_aliases(elements);

  // Check that asking for invalid element returns NULL pointer.
  Ioss::ElementTopology *invalid = Ioss::ElementTopology::factory("Greg", true);
  if (invalid == NULL) {
    OUTPUT << "Testing request for invalid element: "
	 << std::setw(40) << "OK" << '\n';
  } else {
    OUTPUT << "Testing request for invalid element: "
	 << std::setw(40) << "FAIL" << '\n';
    err_count++;
  }
  return err_count;
}

bool test_element(const std::string& type)
{
  // NOTE: For true test, should run with purify checking enabled to
  //       ensure we are not running off the end of arrays...

  bool result = true;
  Ioss::ElementTopology *element = Ioss::ElementTopology::factory(type);
  if (element == NULL) {
    OUTPUT << "ERROR: Element type '" << type << "' could not be constructed.";
    // Must return since we have a NULL pointer and can't do further tests...
    return false;
  }

  // See if the name is an alias for another element (type != name())
  std::string name = element->name();
  OUTPUT << "(" << name << ")\t\t";

  // Check that name is alias for name...
  if (!element->is_alias(type)) {
    OUTPUT << "\n\tName is not valid alias";
    result = false;
  }

  // Check that master element name is an alias...
  if (!element->is_alias(element->master_element_name())) {
    if (element->name() == "edge2d2" || element->name() == "edge2d3") { // kluge
      OUTPUT << "\n\tMaster element name is not valid alias (ignore) ";
    } else {
      OUTPUT << "\n\tMaster element name is not valid alias";
      result = false;
    }
  }

  // Check that the hash id method of selecting the element returns the correct element.
  unsigned int hash_val = Ioss::ElementTopology::get_unique_id(name);
  if (Ioss::ElementTopology::factory(hash_val) != element) {
    OUTPUT << "\n\tElement to hash value conversion is not valid";
    result = false;
  }

  int order     = element->order();

  bool homo_edges = element->edges_similar();
  bool homo_faces = element->faces_similar();

  int nn  = element->number_nodes();
  if (nn <= 0) {
    OUTPUT << "\n\tInvalid node count";
    result = false;
  }

  int ncn = element->number_corner_nodes();
  if (ncn <= 0 || ncn > nn) {
    OUTPUT << "\n\tInvalid corner node count";
    result = false;
  }

  int ne  = element->number_edges();
  if (ne < 0) {
    OUTPUT << "\n\tInvalid edge count";
    result = false;
  }

  int nf  = element->number_faces();
  if (nf < 0) {
    OUTPUT << "\n\tInvalid face count";
    result = false;
  }

  // Verify Euler's Formula holds... V-E+F=2
  if (element->parametric_dimension() == 3) {
    int euler = ncn-ne+nf;
    if (euler != 2) {
      OUTPUT << "\n\tEuler's formula violated (V-E+F=2), value = " << euler << "\n";
      result = false;
    }
  }

  int nne = element->number_nodes_edge(0);
  if (nne == -1) {
    if (homo_edges != false) {
      OUTPUT << "\n\tInconsistent edge homogeneity...\n";
      result = false;
    } else {
      for (int edge = 1; edge <= ne; edge++) {
	int nnei = element->number_nodes_edge(edge);
	if (nnei < 0 || nnei > nn) {
	  OUTPUT << "\n\tInconsistent nodes per edge...\n";
	  result = false;
	}
      }
    }
  } else {
    if (nne < 0 || nne > nn) {
      OUTPUT << "\n\tInconsistent nodes per edge...\n";
      result = false;
    }
  }

  int nnf = element->number_nodes_face(0);
  if (nnf > nn || nnf < -1) {
    OUTPUT << "\n\tInvalid face node count";
    result = false;
  }

  // Check boundary and other topologies...
  if (nf > 0) {
    for (int i=0; i <= nf; i++) {
      Ioss::ElementTopology *face = element->face_type(i);
      if (face == NULL && i > 0) {
	OUTPUT << "\n\tBad face type for face " << i;
	result = false;
      } else if (face == NULL && i == 0 && homo_faces == true) {
	OUTPUT << "\n\tHomogenous faces, but null face_type";
	result = false;
      }	else if (face != NULL) {
	unsigned int nnfi = element->number_nodes_face(i);
	if (nnfi != (unsigned int) face->number_nodes()) {
	  OUTPUT << "\n\tNode count mismatch on face " << i;
	  result = false;
	}
	if (i != 0) {
	  std::vector<int> conn = element->face_connectivity(i);
	  if ((unsigned int)nnfi != conn.size()) {
	    OUTPUT << "\n\tNode count and face connectivity size "
	      "mismatch on face " << i;
	    result = false;
	  }
	}
      }
    }
  }
  // Edges...
  if (ne > 0) {
    for (int i=0; i <= ne; i++) {
      Ioss::ElementTopology *edge = element->edge_type(i);
      if (edge == NULL && i > 0) {
	OUTPUT << "\n\tBad edge type for edge " << i;
	result = false;
      } else if (edge == NULL && i == 0 && homo_edges == true) {
	OUTPUT << "\n\tHomogenous edges, but null edge_type";
	result = false;
      }	else if (edge != NULL) {
	unsigned int nnei = element->number_nodes_edge(i);
	if (nnei != (unsigned int) edge->number_nodes()) {
	  OUTPUT << "\n\tNode count mismatch on edge " << i;
	  result = false;
	}
	if (i != 0) {
	  std::vector<int> conn = element->edge_connectivity(i);
	  if ((unsigned int)nnei != conn.size()) {
	    OUTPUT << "\n\tNode count and edge connectivity size "
	      "mismatch on edge " << i;
	    result = false;
	  }
	}
      }
    }
  }

  // Variable types...
  const VariableType *vt = VariableType::factory(element->name());
  if (vt == NULL) {
    OUTPUT << "\n\tVariable Type does not exist for this name";
    result = false;
  } else {
    // See if component counts match...
    int vt_comp = vt->component_count();
    if (nn != vt_comp) {
      OUTPUT << "\n\tNode count does not match component count";
      result = false;
    }
  }

  // For elements with dimension == 3
  // Get face-edge-order
  // Get face-node-connectivity
  // Foreach edge in face, get nodal connectivity
  //   ensure that node is in face connectivity...
  //
  // For an edge:   1------3------2
  //
  // For a face: Corner nodes are first in connectivity
  //             Center nodes follow.
  // So:
  //
  //                        2		x
  //    3----6----2        / \		x
  //    |         |       /   \		x
  //    7         5      5     4	x
  //    |         |     /       \	x
  //    0----4----1    0----3----1	x
  //
  if (element->parametric_dimension() == 3) {
    for (int i=1; i <= nf; i++) {

      // Nodes defining face...
      std::vector<int> face_con = element->face_connectivity(i);

      unsigned int fncn = element->face_type(i)->number_corner_nodes();
      unsigned int num_edges_face = element->number_edges_face(i);
      if (fncn != num_edges_face) {
	OUTPUT << "\n\tFace corner node count should match edges/face for face "
	     << i;
	result = false;
      }

      // Nodes defining face...
      std::vector<int> face_conn = element->face_connectivity(i);

      // Edges defining face...
      std::vector<int> face_edge_conn = element->face_edge_connectivity(i);
      if ((unsigned int)num_edges_face != face_edge_conn.size()) {
	OUTPUT << "\n\tEdges per face mismatch for face " << i;
	result = false;
      } else {
	for (unsigned int j=0; j < num_edges_face; j++) {
	  // Not implemented in all elements yet...
	  std::vector<int> edge_conn = element->
	    edge_connectivity(face_edge_conn[j]+1);
	  // Check that first two nodes in 'edge_conn' match
	  // corresponding nodes in 'face_con'
	  if ( (edge_conn[0] != face_conn[j] &&
		edge_conn[1] != face_conn[j]) ||
	       (edge_conn[0] != face_conn[(j+1)%fncn] &&
		edge_conn[1] != face_conn[(j+1)%fncn]) ) {
	    OUTPUT << "\n\tEdge Connectivity does not match face "
	      "connectivity for edge " << j+1 << " on face " << i;
	    result = false;
	  }
	  if (order == 2) {
	    if (edge_conn.size() != 3) {
	      OUTPUT << "\n\tInvalid edge connectivity count.";
	      result = false;
	    }
	    if (face_conn.size() < (unsigned int)fncn + num_edges_face) {
	      OUTPUT << "\n\tInvalid face connectivity count.";
	      result = false;
	    }
	    if (edge_conn[2] != face_conn[fncn+j]) {
	      OUTPUT << "\n\tMid-Side Node Edge Connectivity does not match face "
		"connectivity for edge " << j+1 << " on face " << i;
	      result = false;
	    }
	  }
	}
      }
    }
  }
  return result;
}

void test_aliases(const NameList &elements)
{
  int count = elements.size();
  OUTPUT << "\n\nTesting Element Topology Aliases...\n";

  for (int i=0; i < count; i++) {
    Ioss::ElementTopology *el_top = Ioss::ElementTopology::factory(elements[i]);
    if (el_top->name() == elements[i]) {
      OUTPUT << "Element: " << std::setw(20) << elements[i]
	   << "(" << el_top->master_element_name() << ") has the following aliases:\n\t";
      for (int j=0; j < count; j++) {
	if (i != j && el_top->is_alias(elements[j])) {
	  OUTPUT << elements[j] << ", ";
	}
      }
      OUTPUT << std::endl;
    }
  }
}
