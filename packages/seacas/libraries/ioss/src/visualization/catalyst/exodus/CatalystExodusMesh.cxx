// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystExodusMesh.h"
#include "CatalystManager.h"
#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkCompositeDataIterator.h"
#include "vtkDataAssembly.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkInformation.h"
#include "vtkObjectFactory.h"
#include "vtkPartitionedDataSet.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"
#include "vtksys/SystemTools.hxx"
#include <unistd.h>

namespace Iovs_exodus {

  CatalystExodusMesh::CatalystExodusMesh(Iovs::CatalystManager *cm,
                                         CatalystPipelineInfo  &catalystPipelineInfo)
  {

    this->catalystPipelineInfo = catalystPipelineInfo;
    this->catManager           = cm;
    this->global_points        = nullptr;
    this->num_global_points    = 0;
    this->UnderscoreVectors    = true;
    this->ApplyDisplacements   = true;
    this->writeCatalystMesh    = false;

    vtkNew<vtkDataAssembly> assembly;
    assembly->SetRootNodeName(ASSEMBLY_ROOT_NAME.c_str());
    assembly->AddNode(ASSEMBLY_ELEMENT_BLOCKS.c_str());
    this->vpdc->SetDataAssembly(assembly);
  }

  CatalystExodusMesh::~CatalystExodusMesh()
  {
    this->ReleaseGlobalPoints();
    this->ebmap.clear();
    this->ebmap_reverse.clear();
    this->global_elem_id_map.clear();
    this->global_point_id_to_global_elem_id.clear();
    this->ebidmap.clear();
  }

  void CatalystExodusMesh::ReleaseGlobalPoints()
  {
    if (this->global_points) {
      this->global_points->Delete();
      this->global_points = nullptr;
    }
  }

  vtkPartitionedDataSetCollection *CatalystExodusMesh::getPartitionedDataSetCollection()
  {
    return this->vpdc.GetPointer();
  }

  bool CatalystExodusMesh::UnderscoreVectorsON() { return this->UnderscoreVectors; }

  void CatalystExodusMesh::SetUnderscoreVectors(bool status) { this->UnderscoreVectors = status; }

  bool CatalystExodusMesh::ApplyDisplacementsON() { return this->ApplyDisplacements; }

  void CatalystExodusMesh::SetApplyDisplacements(bool status) { this->ApplyDisplacements = status; }

  void CatalystExodusMesh::PerformCoProcessing(std::vector<int>         &error_and_warning_codes,
                                               std::vector<std::string> &error_and_warning_messages)
  {

    this->catManager->PerformCoProcessing(error_and_warning_codes, error_and_warning_messages,
                                          catalystPipelineInfo);
  }

  void CatalystExodusMesh::SetTimeData(double currentTime, int timeStep)
  {
    this->catManager->SetTimeData(currentTime, timeStep, catalystPipelineInfo);
  }

  void CatalystExodusMesh::CreateGlobalVariable(const std::string &variable_name, int num_comps,
                                                const double *data)
  {

    vtkVariant v((double)0.0);
    this->CreateGlobalVariableVariant(variable_name, num_comps, v, data);
  }

  void CatalystExodusMesh::CreateGlobalVariable(const std::string &variable_name, int num_comps,
                                                const int *data)
  {

    vtkVariant v((int)0);
    this->CreateGlobalVariableVariant(variable_name, num_comps, v, data);
  }

  void CatalystExodusMesh::CreateGlobalVariable(const std::string &variable_name, int num_comps,
                                                const int64_t *data)
  {

    vtkVariant v((int64_t)0);
    this->CreateGlobalVariableVariant(variable_name, num_comps, v, data);
  }

  void CatalystExodusMesh::CreateGlobalVariableVariant(const std::string &variable_name,
                                                       int num_comps, vtkVariant &v,
                                                       const void *data)
  {
    for (std::map<int, unsigned int>::iterator iter = this->ebidmap.begin();
         iter != this->ebidmap.end(); ++iter) {
      vtkUnstructuredGrid *ug = getUnstructuredGrid(iter->second);
      this->CreateGlobalVariableInternal(variable_name, num_comps, ug, v, data);
    }
  }

  void CatalystExodusMesh::CreateGlobalVariableInternal(const std::string &variable_name,
                                                        int num_comps, vtkUnstructuredGrid *ug,
                                                        vtkVariant &v, const void *data)
  {
    if (ug == nullptr) {
      return;
    }
    vtkFieldData *fieldData = ug->GetFieldData();
    vtkDataArray *da        = fieldData->GetArray(variable_name.c_str());
    if (da) {
      fieldData->RemoveArray(variable_name.c_str());
    }
    vtkDataArray *arr = vtkDataArray::CreateDataArray(v.GetType());
    arr->SetName(variable_name.c_str());
    arr->SetNumberOfComponents(num_comps);
    arr->SetNumberOfTuples(1);
    fieldData->AddArray(arr);
    arr->Delete();

    for (int cc = 0; cc < num_comps; cc++) {
      arr->SetComponent(0, cc, this->GetArrayValue(v, data, cc));
    }
  }

  void CatalystExodusMesh::InitializeGlobalPoints(int num_points, int dimension, const double *data)
  {

    this->global_point_id_to_global_elem_id.resize(num_points);
    this->num_global_points = num_points;
    this->global_points     = vtkPoints::New();
    this->global_points->SetNumberOfPoints(num_points);
    vtkDoubleArray *coords = vtkDoubleArray::New();
    coords->SetNumberOfComponents(3);
    coords->SetNumberOfTuples(num_points);
    int index = 0;
    for (int i = 0; i < num_points; i++) {
      coords->SetComponent(i, 0, data[index++]);
      coords->SetComponent(i, 1, data[index++]);
      if (dimension != 2)
        coords->SetComponent(i, 2, data[index++]);
      else
        coords->SetComponent(i, 2, 0.0);
    }
    this->global_points->SetData(coords);
    coords->Delete();
  }

  void
  CatalystExodusMesh::InitializeElementBlocks(const ElementBlockIdNameList &elemBlockNameIdList)
  {
    this->ReleaseGlobalPoints();
    this->ebmap.clear();
    this->ebmap_reverse.clear();
    this->global_elem_id_map.clear();
    this->global_point_id_to_global_elem_id.clear();
    this->ebidmap.clear();

    for (unsigned int i = 0; i < elemBlockNameIdList.size(); i++) {
      this->ebidmap[elemBlockNameIdList[i].first] = i;
      unsigned int                  pdsIdx        = vpdc->GetNumberOfPartitionedDataSets();
      vtkNew<vtkPartitionedDataSet> pds;
      vpdc->SetPartitionedDataSet(pdsIdx, pds);
      vpdc->GetMetaData(pdsIdx)->Set(vtkCompositeDataSet::NAME(), elemBlockNameIdList[i].second);
      auto assembly = vpdc->GetDataAssembly();
      auto node =
          assembly->AddNode(elemBlockNameIdList[i].second.c_str(), getElementBlocksAssemblyNode());
      assembly->SetAttribute(node, ASSEMBLY_LABEL.c_str(), elemBlockNameIdList[i].second.c_str());
      assembly->AddDataSetIndex(node, pdsIdx);
    }
  }

  void CatalystExodusMesh::CreateElementBlock(const char *elem_block_name, int elem_block_id,
                                              const std::string &elem_type, int nodes_per_elem,
                                              int num_elem, const int64_t *global_elem_ids,
                                              int *connectivity)
  {

    vtkVariant v((int)0);
    this->CreateElementBlockInternal(elem_block_name, elem_block_id, elem_type, nodes_per_elem,
                                     num_elem, v, global_elem_ids, connectivity);
  }

  void CatalystExodusMesh::CreateElementBlock(const char *elem_block_name, int elem_block_id,
                                              const std::string &elem_type, int nodes_per_elem,
                                              int num_elem, const int64_t *global_elem_ids,
                                              int64_t *connectivity)
  {

    vtkVariant v((int64_t)0);
    this->CreateElementBlockInternal(elem_block_name, elem_block_id, elem_type, nodes_per_elem,
                                     num_elem, v, global_elem_ids, connectivity);
  }

  void CatalystExodusMesh::CreateElementBlockInternal(const char        *elem_block_name,
                                                      int                elem_block_id,
                                                      const std::string &elem_type,
                                                      int nodes_per_elem, int num_elem,
                                                      vtkVariant &v, const int64_t *global_elem_ids,
                                                      void *connectivity)
  {

    std::string elemType(vtksys::SystemTools::UpperCase(elem_type));

    int point_count;
    int vtk_type;

    // Check for quadratic elements
    if ((elemType.substr(0, 3) == "TRI") && (nodes_per_elem == 6)) {
      vtk_type    = VTK_QUADRATIC_TRIANGLE;
      point_count = 6;
    }
    else if ((elemType.substr(0, 3) == "SHE") && (nodes_per_elem == 8)) {
      vtk_type    = VTK_QUADRATIC_QUAD;
      point_count = 8;
    }
    else if ((elemType.substr(0, 3) == "SHE") && (nodes_per_elem == 9)) {
      vtk_type    = VTK_QUADRATIC_QUAD;
      point_count = 8;
    }
    else if ((elemType.substr(0, 3) == "TET") && (nodes_per_elem == 10)) {
      vtk_type    = VTK_QUADRATIC_TETRA;
      point_count = 10;
    }
    else if ((elemType.substr(0, 3) == "TET") && (nodes_per_elem == 11)) {
      vtk_type    = VTK_QUADRATIC_TETRA;
      point_count = 10;
    }
    else if ((elemType.substr(0, 3) == "WED") && (nodes_per_elem == 15)) {
      vtk_type    = VTK_QUADRATIC_WEDGE;
      point_count = 15;
    }
    else if (((elemType.substr(0, 3) == "HEX") && (nodes_per_elem == 20)) ||
             ((elemType.substr(0, 3) == "HEX") && (nodes_per_elem == 21))) {
      vtk_type    = VTK_QUADRATIC_HEXAHEDRON;
      point_count = 20;
    }
    else if ((elemType.substr(0, 3) == "HEX") && (nodes_per_elem == 27)) {
      vtk_type    = VTK_TRIQUADRATIC_HEXAHEDRON;
      point_count = 27;
    }
    else if ((elemType.substr(0, 3) == "QUA") && (nodes_per_elem == 8)) {
      vtk_type    = VTK_QUADRATIC_QUAD;
      point_count = 8;
    }
    else if ((elemType.substr(0, 3) == "QUA") && (nodes_per_elem == 9)) {
      vtk_type    = VTK_BIQUADRATIC_QUAD;
      point_count = 9;
    }
    else if ((elemType.substr(0, 3) == "TRU") && (nodes_per_elem == 3)) {
      vtk_type    = VTK_QUADRATIC_EDGE;
      point_count = 3;
    }
    else if ((elemType.substr(0, 3) == "BEA") && (nodes_per_elem == 3)) {
      vtk_type    = VTK_QUADRATIC_EDGE;
      point_count = 3;
    }
    else if ((elemType.substr(0, 3) == "BAR") && (nodes_per_elem == 3)) {
      vtk_type    = VTK_QUADRATIC_EDGE;
      point_count = 3;
    }
    else if ((elemType.substr(0, 3) == "EDG") && (nodes_per_elem == 3)) {
      vtk_type    = VTK_QUADRATIC_EDGE;
      point_count = 3;
    }

    // Check for linear elements
    else if (elemType.substr(0, 3) == "CIR") {
      vtk_type    = VTK_VERTEX;
      point_count = 1;
    }
    else if (elemType.substr(0, 3) == "SPH") {
      vtk_type    = VTK_VERTEX;
      point_count = 1;
    }
    else if (elemType.substr(0, 3) == "BAR") {
      vtk_type    = VTK_LINE;
      point_count = 2;
    }
    else if (elemType.substr(0, 3) == "TRU") {
      vtk_type    = VTK_LINE;
      point_count = 2;
    }
    else if (elemType.substr(0, 3) == "BEA") {
      vtk_type    = VTK_LINE;
      point_count = 2;
    }
    else if (elemType.substr(0, 3) == "EDG") {
      vtk_type    = VTK_LINE;
      point_count = 2;
    }
    else if (elemType.substr(0, 3) == "TRI") {
      vtk_type    = VTK_TRIANGLE;
      point_count = 3;
    }
    else if (elemType.substr(0, 3) == "QUA") {
      vtk_type    = VTK_QUAD;
      point_count = 4;
    }
    else if (elemType.substr(0, 3) == "TET") {
      vtk_type    = VTK_TETRA;
      point_count = 4;
    }
    else if (elemType.substr(0, 3) == "PYR") {
      vtk_type    = VTK_PYRAMID;
      point_count = 5;
    }
    else if (elemType.substr(0, 3) == "WED") {
      vtk_type    = VTK_WEDGE;
      point_count = 6;
    }
    else if (elemType.substr(0, 3) == "HEX") {
      vtk_type    = VTK_HEXAHEDRON;
      point_count = 8;
    }
    else if (elemType.substr(0, 3) == "NSI") {
      vtk_type    = VTK_POLYGON;
      point_count = 0;
    }
    else if (elemType.substr(0, 3) == "NFA") {
      vtk_type    = VTK_POLYHEDRON;
      point_count = 0;
    }
    else if ((elemType.substr(0, 3) == "SHE") && (nodes_per_elem == 3)) {
      vtk_type    = VTK_TRIANGLE;
      point_count = 3;
    }
    else if ((elemType.substr(0, 3) == "SHE") && (nodes_per_elem == 4)) {
      vtk_type    = VTK_QUAD;
      point_count = 4;
    }
    else if ((elemType.substr(0, 8) == "STRAIGHT") && (nodes_per_elem == 2)) {
      vtk_type    = VTK_LINE;
      point_count = 2;
    }
    else if (elemType.substr(0, 3) == "SUP") {
      vtk_type    = VTK_POLY_VERTEX;
      point_count = nodes_per_elem;
    }
    else if ((elemType.substr(0, 4) == "NULL") && (num_elem == 0)) {
      std::cout << "NULL element block found";
      std::cout << "Unable to create element block " << elem_block_name;
      return; // silently ignore empty element blocks
    }
    else {
      std::cout << "Unsupported element type: " << elemType.c_str();
      std::cout << "Unable to create element block " << elem_block_name;
      return;
      // cell types not currently handled
      // quadratic wedge - 15,16 nodes
      // quadratic pyramid - 13 nodes
    }

    vtkIdType cell_vertex_order[point_count];
    for (int p = 0; p < point_count; p++) {
      cell_vertex_order[p] = p;
    }

    if (vtk_type == VTK_QUADRATIC_WEDGE) {
      cell_vertex_order[12] = 9;
      cell_vertex_order[13] = 10;
      cell_vertex_order[14] = 11;

      cell_vertex_order[9]  = 12;
      cell_vertex_order[10] = 13;
      cell_vertex_order[11] = 14;
    }
    else if (vtk_type == VTK_QUADRATIC_HEXAHEDRON) {
      cell_vertex_order[16] = 12;
      cell_vertex_order[17] = 13;
      cell_vertex_order[18] = 14;
      cell_vertex_order[19] = 15;

      cell_vertex_order[12] = 16;
      cell_vertex_order[13] = 17;
      cell_vertex_order[14] = 18;
      cell_vertex_order[15] = 19;
    }
    else if (vtk_type == VTK_TRIQUADRATIC_HEXAHEDRON) {
      cell_vertex_order[16] = 12;
      cell_vertex_order[17] = 13;
      cell_vertex_order[18] = 14;
      cell_vertex_order[19] = 15;

      cell_vertex_order[12] = 16;
      cell_vertex_order[13] = 17;
      cell_vertex_order[14] = 18;
      cell_vertex_order[15] = 19;

      cell_vertex_order[23] = 20;
      cell_vertex_order[24] = 21;
      cell_vertex_order[25] = 22;
      cell_vertex_order[26] = 23;

      cell_vertex_order[21] = 24;
      cell_vertex_order[22] = 25;
      cell_vertex_order[20] = 26;
    }

    vtkIdType pts[point_count];
    int       index = 0;
    while (index < num_elem * point_count) {
      for (int p = 0; p < point_count; p++) {
        int index_orig = index;

        int64_t conn_ind                     = this->GetArrayValue(v, connectivity, index_orig) - 1;
        this->ebmap[elem_block_id][conn_ind] = -1;
        index++;
      }
      if (point_count < nodes_per_elem) {
        for (int p = 0; p < (nodes_per_elem - point_count); p++) {
          index++;
        }
      }
    }

    vtkUnstructuredGrid *ug     = vtkUnstructuredGrid::New();
    vtkPoints           *points = vtkPoints::New();
    double               x[3];
    vtkIdType            i = 0;
    for (std::map<int, int>::iterator ii = this->ebmap[elem_block_id].begin();
         ii != this->ebmap[elem_block_id].end(); ++ii) {

      this->global_points->GetPoint((*ii).first, x);
      (*ii).second                          = i;
      this->ebmap_reverse[elem_block_id][i] = (*ii).first;
      points->InsertNextPoint(x);
      i++;
    }

    index = 0;
    std::vector<int> object_ids;
    while (index < num_elem * point_count) {
      for (int p = 0; p < point_count; p++) {
        int64_t conn_ind = this->GetArrayValue(v, connectivity, index) - 1;
        this->global_point_id_to_global_elem_id[conn_ind] = global_elem_ids[ug->GetNumberOfCells()];

        conn_ind                  = this->GetArrayValue(v, connectivity, index++) - 1;
        pts[cell_vertex_order[p]] = this->ebmap[elem_block_id][conn_ind];
      }
      if (point_count < nodes_per_elem) {
        for (int p = 0; p < (nodes_per_elem - point_count); p++) {
          int64_t conn_ind = this->GetArrayValue(v, connectivity, index) - 1;
          this->global_point_id_to_global_elem_id[conn_ind] =
              global_elem_ids[ug->GetNumberOfCells()];
          index++;
        }
      }

      this->global_elem_id_map[elem_block_id][global_elem_ids[ug->GetNumberOfCells()]] =
          ug->GetNumberOfCells();

      ug->InsertNextCell(vtk_type, point_count, pts);
      object_ids.push_back(this->ebidmap[elem_block_id]);
    }

    ug->SetPoints(points);

    vpdc->GetPartitionedDataSet(this->ebidmap[elem_block_id])
        ->SetPartition(PDS_UNSTRUCTURED_GRID_INDEX, ug);

    ug->Delete();
    points->Delete();
  }

  void CatalystExodusMesh::CreateElementVariable(const std::string &variable_name, int num_comps,
                                                 int elem_block_id, const double *data)
  {

    vtkVariant v((double)0.0);
    this->CreateElementVariableVariant(variable_name, num_comps, elem_block_id, v, data);
  }

  void CatalystExodusMesh::CreateElementVariable(const std::string &variable_name, int num_comps,
                                                 int elem_block_id, const int *data)
  {

    vtkVariant v((int)0);
    this->CreateElementVariableVariant(variable_name, num_comps, elem_block_id, v, data);
  }

  void CatalystExodusMesh::CreateElementVariable(const std::string &variable_name, int num_comps,
                                                 int elem_block_id, const int64_t *data)
  {

    vtkVariant v((int64_t)0);
    this->CreateElementVariableVariant(variable_name, num_comps, elem_block_id, v, data);
  }

  void CatalystExodusMesh::CreateElementVariableVariant(const std::string &variable_name,
                                                        int num_comps, int elem_block_id,
                                                        vtkVariant &v, const void *data)
  {
    vtkUnstructuredGrid *ug = getUnstructuredGrid(this->ebidmap[elem_block_id]);
    this->CreateElementVariableInternal(variable_name, num_comps, ug, v, data);
  }

  void CatalystExodusMesh::CreateElementVariableInternal(const std::string &variable_name,
                                                         int num_comps, vtkUnstructuredGrid *ug,
                                                         vtkVariant &v, const void *data)
  {
    if (ug == nullptr) {
      return;
    }
    vtkFieldData *cell_data = ug->GetCellData();
    vtkDataArray *da        = cell_data->GetArray(variable_name.c_str());
    if (da) {
      cell_data->RemoveArray(variable_name.c_str());
    }
    vtkDataArray *arr = vtkDataArray::CreateDataArray(v.GetType());
    arr->SetName(variable_name.c_str());
    arr->SetNumberOfComponents(num_comps);
    arr->SetNumberOfTuples(ug->GetNumberOfCells());
    cell_data->AddArray(arr);
    arr->Delete();

    for (int ii = 0; ii < ug->GetNumberOfCells(); ii++) {
      for (int cc = 0; cc < num_comps; cc++) {
        arr->SetComponent(ii, cc, this->GetArrayValue(v, data, num_comps * ii + cc));
      }
    }
  }

  void CatalystExodusMesh::CreateNodalVariable(const std::string &variable_name, int num_comps,
                                               const double *data)
  {

    vtkVariant v((double)0.0);
    this->CreateNodalVariableVariant(variable_name, num_comps, v, data);
  }

  void CatalystExodusMesh::CreateNodalVariable(const std::string &variable_name, int num_comps,
                                               const int *data)
  {

    vtkVariant v((int)0);
    this->CreateNodalVariableVariant(variable_name, num_comps, v, data);
  }

  void CatalystExodusMesh::CreateNodalVariable(const std::string &variable_name, int num_comps,
                                               const int64_t *data)
  {

    vtkVariant v((int64_t)0);
    this->CreateNodalVariableVariant(variable_name, num_comps, v, data);
  }

  void CatalystExodusMesh::CreateNodalVariableVariant(const std::string &variable_name,
                                                      int num_comps, vtkVariant &v,
                                                      const void *data)
  {
    for (std::map<int, unsigned int>::iterator iter = this->ebidmap.begin();
         iter != this->ebidmap.end(); ++iter) {
      vtkUnstructuredGrid *ug = getUnstructuredGrid(iter->second);
      this->CreateNodalVariableInternal(variable_name, num_comps, ug, iter->first, this->ebmap, v,
                                        data);
    }
  }

  void CatalystExodusMesh::CreateNodalVariableInternal(const std::string &variable_name,
                                                       int num_comps, vtkUnstructuredGrid *ug,
                                                       int element_block_id,
                                                       std::map<int, std::map<int, int>> &point_map,
                                                       vtkVariant &v, const void *data)
  {
    if (ug == nullptr) {
      return;
    }

    bool displace_nodes = false;
    if ((variable_name.substr(0, 3) == "DIS") || (variable_name.substr(0, 3) == "dis")) {
      displace_nodes = true;
    }

    if (ug->GetNumberOfCells() == 0) {
      return;
    }
    vtkFieldData *point_data = ug->GetPointData();
    vtkDataArray *da         = point_data->GetArray(variable_name.c_str());
    if (da) {
      point_data->RemoveArray(variable_name.c_str());
    }
    vtkDataArray *arr = vtkDataArray::CreateDataArray(v.GetType());
    arr->SetName(variable_name.c_str());
    arr->SetNumberOfComponents(num_comps);
    arr->SetNumberOfTuples(ug->GetPoints()->GetNumberOfPoints());
    point_data->AddArray(arr);
    arr->Delete();

    int index = 0;
    for (int i = 0; i < this->num_global_points; i++) {
      std::map<int, int>::iterator mit = point_map[element_block_id].find(i);

      for (int c = 0; c < num_comps; c++) {
        if (mit != point_map[element_block_id].end()) {
          arr->SetComponent(point_map[element_block_id][i], c,
                            this->GetArrayValue(v, data, index++));
        }
        else {
          index++;
        }
      }

      if (displace_nodes) {
        if (mit != point_map[element_block_id].end()) {
          vtkPoints *points = ug->GetPoints();
          double     x[3];
          this->global_points->GetPoint(i, x);
          for (int c = 0; c < num_comps; c++) {
            x[c] += arr->GetComponent(point_map[element_block_id][i], c);
          }
          points->SetPoint(point_map[element_block_id][i], x);
        }
      }
    }
  }

  void CatalystExodusMesh::ReleaseMemory()
  {
    for (std::map<int, unsigned int>::iterator iter = this->ebidmap.begin();
         iter != this->ebidmap.end(); ++iter) {
      vtkUnstructuredGrid *ug = getUnstructuredGrid(iter->second);
      this->ReleaseMemoryInternal(ug);
    }

    this->ebmap_reverse.clear();
    this->global_elem_id_map.clear();
    this->global_point_id_to_global_elem_id.clear();

    this->catManager->WriteToLogFile(catalystPipelineInfo);
  }

  void CatalystExodusMesh::ReleaseMemoryInternal(vtkUnstructuredGrid *ug)
  {
    if (ug == nullptr) {
      return;
    }
    vtkFieldData *data               = ug->GetCellData();
    vtkDataArray *global_element_ids = nullptr;
    vtkDataArray *object_id          = nullptr;
    while (data->GetNumberOfArrays() > 0) {
      if (std::string(data->GetArray(0)->GetName()) == "ids") {
        global_element_ids = data->GetArray(0);
        global_element_ids->Register(0);
      }
      if (std::string(data->GetArray(0)->GetName()) == "object_id") {
        object_id = data->GetArray(0);
        object_id->Register(0);
      }
      data->RemoveArray(data->GetArray(0)->GetName());
    }

    if (global_element_ids) {
      data->AddArray(global_element_ids);
      ug->GetCellData()->SetActiveGlobalIds("ids");
      global_element_ids->Delete();
    }

    if (object_id) {
      data->AddArray(object_id);
      object_id->Delete();
    }

    data                          = ug->GetPointData();
    vtkDataArray *global_node_ids = 0;
    while (data->GetNumberOfArrays() > 0) {
      if (std::string(data->GetArray(0)->GetName()) == "ids") {
        global_node_ids = data->GetArray(0);
        global_node_ids->Register(0);
      }
      data->RemoveArray(data->GetArray(0)->GetName());
    }

    if (global_node_ids) {
      data->AddArray(global_node_ids);
      ug->GetPointData()->SetActiveGlobalIds("ids");
      global_node_ids->Delete();
    }
  }

  void CatalystExodusMesh::logMemoryUsageAndTakeTimerReading()
  {
    this->catManager->logMemoryUsageAndTakeTimerReading(catalystPipelineInfo);
  }

  void CatalystExodusMesh::Delete() { this->catManager->DeletePipeline(catalystPipelineInfo); }

  double CatalystExodusMesh::GetArrayValue(vtkVariant &v, const void *data, int index)
  {
    if (v.IsDouble()) {
      return ((double)static_cast<const double *>(data)[index]);
    }
    else if (v.IsFloat()) {
      return ((double)static_cast<const float *>(data)[index]);
    }
    else if (v.IsInt()) {
      return ((double)static_cast<const int *>(data)[index]);
    }
    else if (v.IsLongLong()) {
      return ((double)static_cast<const long long *>(data)[index]);
    }
    else if (v.IsLong()) {
      return ((double)static_cast<const long *>(data)[index]);
    }
    else {
      std::cout << "Unhandled type found in GetArrayValue: " << v.GetTypeAsString();
      return 0.0;
    }
  }

  int CatalystExodusMesh::getElementBlocksAssemblyNode()
  {
    auto assembly = vpdc->GetDataAssembly();
    return assembly->GetFirstNodeByPath(
        ("/" + ASSEMBLY_ROOT_NAME + "/" + ASSEMBLY_ELEMENT_BLOCKS).c_str());
  }

  vtkUnstructuredGrid *CatalystExodusMesh::getUnstructuredGrid(int blockId)
  {
    vtkDataSet *ds =
        vpdc->GetPartitionedDataSet(blockId)->GetPartition(PDS_UNSTRUCTURED_GRID_INDEX);
    return vtkUnstructuredGrid::SafeDownCast(ds);
  }

} // namespace Iovs_exodus
