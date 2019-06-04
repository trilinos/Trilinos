/*
 * Copyright(C) 2010-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef __VTK_EXODUS_II_MULTI_BLOCK_DATA_SET_H
#define __VTK_EXODUS_II_MULTI_BLOCK_DATA_SET_H

#include "vtkMultiBlockDataSet.h"
#include <map>
#include <string>
#include <vector>

class vtkVariant;
class vtkPoints;
class Ve2mSideSetInfo;

class vtkExodusIIMultiBlockDataSet : public vtkMultiBlockDataSet
{
public:
  vtkTypeMacro(vtkExodusIIMultiBlockDataSet, vtkMultiBlockDataSet);
  void PrintSelf(ostream &os, vtkIndent indent);

  static vtkExodusIIMultiBlockDataSet *New();

  // Description:
  // Initializes the element blocks to NULL data sets with ids in element_block_id_list.
  // This method must be called first.
  void InitializeElementBlocks(const std::vector<int> &element_block_id_list);

  // Description:
  // Initializes the vtkMultiBlockDataSet with a global array of points
  // defined by num_points, dimension (2,3), and data.  Clears any existing data.
  void InitializeGlobalPoints(int num_points, int dimension, const double *data);

  // Description:
  // Clears all nodal and element variables from the vtkMultiBlockDataSet.
  // Clears the global vtkPoints.
  void ReleaseMemory();

  // Description:
  // Creates a global variable on the vtkExodusIIMultiBlockDataSet.
  // Input argument v must have the same type as the data
  // contained in input the argument array data. Creates
  // the global variable on all element blocks.
  void CreateGlobalVariable(std::vector<std::string> &component_names, vtkVariant &v,
                            const void *data);

  // Description:
  // Creates a vtkUnstructuredGrid on the vtkExodusIIMultiBlockDataSet
  // that represents and element block in the Exodus II data.  The global_points
  // array contains all of the points in the Exodus II file.
  void CreateElementBlock(const char *elem_block_name, int elem_block_id,
                          const std::string &elem_type, int nodes_per_elem, int num_elem,
                          vtkVariant &v, const int64_t *global_elem_ids, void *connectivity);

  // Description:
  // Creates an element variable the vtkExodusIIMultiBlockDataSet.
  void CreateElementVariable(std::vector<std::string> &component_names, int elem_block_id,
                             vtkVariant &v, const void *data);

  // Description:
  // Creates a nodal variable the vtkExodusIIMultiBlockDataSet.
  void CreateNodalVariable(std::vector<std::string> &component_names, vtkVariant &v,
                           const void *data);

  // Description:
  // Creates a vtkUnstructuredGrid representing the node set in the Exodus II
  // data. Node sets are arbitrary lists of mesh point ids.
  void CreateNodeSet(const char *node_set_name, int node_set_id, int num_ids, vtkVariant &v,
                     const void *ids);

  // Description:
  // Creates a vtkUnstructuredGrid representing the side set (also Side Block) in
  // the Exodus II data. Side sets are collections of element faces and edges.
  void CreateSideSet(/*const char* side_set_name,*/
                     const char *ss_owner_name, int side_set_id, int num_ids, vtkVariant &v,
                     const void *element_ids, const void *face_ids);

  // Description:
  // If true (the default), vector variables will contain a
  // trailing underscore in their name.  The default behavior
  // is consistent with the ParaView Exodus II file reader.
  vtkSetMacro(UnderscoreVectors, int);
  vtkGetMacro(UnderscoreVectors, int);
  vtkBooleanMacro(UnderscoreVectors, int);

  // Description:
  // If true (the default), displacements will be applied to the
  // mesh nodes before being sent to the in-situ pipeline.  The node
  // displacement variable is called either DISPL or displ.  The
  // default behavior is consistent with the ParaView Exodus II
  // file reader.
  vtkSetMacro(ApplyDisplacements, int);
  vtkGetMacro(ApplyDisplacements, int);
  vtkBooleanMacro(ApplyDisplacements, int);

protected:
  vtkExodusIIMultiBlockDataSet();
  ~vtkExodusIIMultiBlockDataSet();

private:
  vtkExodusIIMultiBlockDataSet(const vtkExodusIIMultiBlockDataSet &); // Not implemented.
  void operator=(const vtkExodusIIMultiBlockDataSet &);               // Not implemented.
  std::map<int, std::map<int, int>> ebmap;
  std::map<int, std::map<int, int>> ebmap_reverse;
  std::map<int, std::map<int, int>> global_elem_id_map;
  std::vector<int> global_point_id_to_global_elem_id;
  std::map<int, unsigned int> ebidmap;
  std::map<int, unsigned int> nsidmap;
  std::map<int, std::map<int, int>> nsmap;
  std::map<int, unsigned int> ssidmap;

  // ssinfomap is used to help track when we see a new sideset. CreateSideSet
  // is called once for each sideset for each block which the sideset spans,
  // and we combine those into a single sideset entity in the vtk
  // representation; also lets us do the necessary bookkeeping to combine
  // the data from the different blocks into the same sideset
  std::map<int, Ve2mSideSetInfo *> ssinfomap;

  std::map<int, std::map<int, int>> ssmap;
  void ContainsVector(std::vector<std::string> &component_names,
                      std::vector<std::string> &prefix_name);
  double GetArrayValue(vtkVariant &v, const void *data, int index);
  void       ReleaseGlobalPoints();
  vtkPoints *global_points;
  int        num_global_points;
  int        UnderscoreVectors;
  int        ApplyDisplacements;
  void CreateGlobalVariableInternal(std::vector<std::string> &component_names,
                                    vtkMultiBlockDataSet *eb, unsigned int bid, vtkVariant &v,
                                    const void *data);
  void CreateNodalVariableInternal(std::vector<std::string> &component_names,
                                   vtkMultiBlockDataSet *eb, std::map<int, unsigned int> &id_map,
                                   std::map<int, std::map<int, int>> &point_map, vtkVariant &v,
                                   const void *data);
  void CreateElementVariableInternal(std::vector<std::string> &component_names,
                                     vtkMultiBlockDataSet *eb, unsigned int bid, vtkVariant &v,
                                     const void *data);
  void ReleaseMemoryInternal(vtkMultiBlockDataSet *eb);
};

// see ssinfomap above, lets us combine sidesets which span multiple
// blocks int a single sideset entity
class Ve2mSideSetInfo
{
public:
  Ve2mSideSetInfo();
  ~Ve2mSideSetInfo();

  int bid;
  std::map<int, int> unique_points;
  std::vector<int> object_ids;
};

#endif /* __VTK_EXODUS_II_MULTI_BLOCK_DATA_SET_H */
