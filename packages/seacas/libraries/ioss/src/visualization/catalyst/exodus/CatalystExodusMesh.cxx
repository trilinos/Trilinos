// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystExodusMesh.h"
#include "CatalystManager.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkCompositeDataIterator.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkInformation.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"
#include "vtksys/SystemTools.hxx"
#include <unistd.h>

namespace Iovs_exodus {

CatalystExodusMesh::CatalystExodusMesh(Iovs::CatalystManager *cm,
    CatalystPipelineInfo& catalystPipelineInfo) {

    this->catalystPipelineInfo = catalystPipelineInfo;
    this->multiBlock = vtkMultiBlockDataSet::New();
    this->catManager = cm;
    this->global_points = nullptr;
    this->num_global_points = 0;
    this->UnderscoreVectors = true;
    this->ApplyDisplacements = true;
    this->writeCatalystMesh = false;
}

CatalystExodusMesh::~CatalystExodusMesh() {
    this->ReleaseGlobalPoints();
    this->ebmap.clear();
    this->ebmap_reverse.clear();
    this->global_elem_id_map.clear();
    this->global_point_id_to_global_elem_id.clear();
    this->ebidmap.clear();
    this->nsmap.clear();
    this->nsidmap.clear();
    this->ssmap.clear();
    this->ssidmap.clear();
    std::map<int, Ve2mSideSetInfo *>::iterator ssinfomapiter;
    for (ssinfomapiter = this->ssinfomap.begin(); ssinfomapiter != this->ssinfomap.end();
         ssinfomapiter++) {
        if (ssinfomapiter->second != NULL) {
            delete ssinfomapiter->second;
        }
    }
    this->ssinfomap.clear();
    this->multiBlock->Delete();
    this->multiBlock = nullptr;
}

void CatalystExodusMesh::ReleaseGlobalPoints() {
    if (this->global_points) {
        this->global_points->Delete();
        this->global_points = nullptr;
    }
}

vtkMultiBlockDataSet* CatalystExodusMesh::getMultiBlockDataSet() {
    return this->multiBlock;
}

bool CatalystExodusMesh::UnderscoreVectorsON() {
    return this->UnderscoreVectors;
}

void CatalystExodusMesh::SetUnderscoreVectors(bool status) {
    this->UnderscoreVectors = status;
}

bool CatalystExodusMesh::ApplyDisplacementsON() {
    return this->ApplyDisplacements;
}

void CatalystExodusMesh::SetApplyDisplacements(bool status) {
    this->ApplyDisplacements = status;
}

void CatalystExodusMesh::PerformCoProcessing(
    std::vector<int> &error_and_warning_codes,
        std::vector<std::string> &error_and_warning_messages) {

    this->catManager->PerformCoProcessing(error_and_warning_codes,
        error_and_warning_messages, catalystPipelineInfo);
}

void CatalystExodusMesh::SetTimeData(double currentTime, int timeStep) {
    this->catManager->SetTimeData(currentTime, timeStep,
        catalystPipelineInfo);
}

void CatalystExodusMesh::CreateGlobalVariable(
    std::vector<std::string> &component_names, const double *data) {

    vtkVariant v((double)0.0);
    this->CreateGlobalVariableVariant(component_names, v, data);
}

void CatalystExodusMesh::CreateGlobalVariable(
    std::vector<std::string> &component_names, const int *data) {

    vtkVariant v((int)0);
    this->CreateGlobalVariableVariant(component_names, v, data);
}

void CatalystExodusMesh::CreateGlobalVariableVariant(
    std::vector<std::string> &component_names, vtkVariant &v,
        const void *data) {

    vtkMultiBlockDataSet *eb =
        vtkMultiBlockDataSet::SafeDownCast(\
            this->multiBlock->GetBlock(ELEMENT_BLOCK_MBDS_ID));

    for (std::map<int, unsigned int>::iterator iter = this->ebidmap.begin();
         iter != this->ebidmap.end(); ++iter) {
        this->CreateGlobalVariableInternal(component_names, eb,
            iter->second, v, data);
    }

/*
    eb = vtkMultiBlockDataSet::SafeDownCast(\
        this->multiBlock->GetBlock(NODE_SETS_MBDS_ID));

    for (std::map<int, unsigned int>::iterator iter = this->nsidmap.begin();
         iter != this->nsidmap.end(); ++iter) {
        this->CreateGlobalVariableInternal(component_names, eb,
            iter->second, v, data);
    }

    eb = vtkMultiBlockDataSet::SafeDownCast(\
        this->multiBlock->GetBlock(SIDE_SETS_MBDS_ID));

    for (std::map<int, unsigned int>::iterator iter = this->ssidmap.begin();
         iter != this->ssidmap.end(); ++iter) {
        this->CreateGlobalVariableInternal(component_names, eb,
            iter->second, v, data);
    }
*/
}

void CatalystExodusMesh::CreateGlobalVariableInternal(
    std::vector<std::string> &component_names, vtkMultiBlockDataSet *eb,
        unsigned int bid, vtkVariant &v, const void *data) {

    int number_data_components = 1;
    std::vector<std::string> prefix_name;
    this->ContainsVector(component_names, prefix_name);
    std::vector<std::string> component_names_buffer = component_names;
    if (prefix_name.size() == 1) {
        number_data_components = component_names.size();
        component_names_buffer = prefix_name;
    }

    vtkUnstructuredGrid *ug = vtkUnstructuredGrid::SafeDownCast(eb->GetBlock(bid));
    vtkFieldData *field_data = ug->GetFieldData();
    std::vector<vtkDataArray *> data_arrays;
    for (std::vector<std::string>::iterator it = component_names_buffer.begin();
        it != component_names_buffer.end(); ++it) {

        vtkDataArray *da = field_data->GetArray((*it).c_str());
        if (da) {
            field_data->RemoveArray((*it).c_str());
        }
        vtkDataArray *arr = vtkDataArray::CreateDataArray(v.GetType());
        arr->SetName((*it).c_str());
        arr->SetNumberOfComponents(number_data_components);
        arr->SetNumberOfTuples(1);
        field_data->AddArray(arr);
        data_arrays.push_back(arr);
        arr->Delete();
    }

    int index = 0;
    for (std::vector<vtkDataArray *>::iterator it = data_arrays.begin();
        it != data_arrays.end(); ++it) {

        for (int c = 0; c < number_data_components; c++) {
            (*it)->SetComponent(0, c, this->GetArrayValue(v, data, index++));
        }
    }
}

void CatalystExodusMesh::InitializeGlobalPoints(int num_points, int dimension,
    const double *data) {

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

void CatalystExodusMesh::InitializeElementBlocks(
    const std::vector<int> &element_block_id_list) {

    this->multiBlock->Initialize();
    this->ReleaseGlobalPoints();
    this->ebmap.clear();
    this->ebmap_reverse.clear();
    this->global_elem_id_map.clear();
    this->global_point_id_to_global_elem_id.clear();
    this->ebidmap.clear();
    this->nsmap.clear();
    this->nsidmap.clear();
    this->ssmap.clear();
    this->ssidmap.clear();

    vtkMultiBlockDataSet *eb = vtkMultiBlockDataSet::New();
    this->multiBlock->SetBlock(ELEMENT_BLOCK_MBDS_ID, eb);
    this->multiBlock->GetMetaData(ELEMENT_BLOCK_MBDS_ID)->\
        Set(vtkCompositeDataSet::NAME(), ELEMENT_BLOCK_MBDS_NAME);

    for (unsigned int i = 0; i < element_block_id_list.size(); i++) {
        this->ebidmap[element_block_id_list[i]] = i;
        vtkUnstructuredGrid *ug = vtkUnstructuredGrid::New();
        eb->SetBlock(i, ug);
        ug->Delete();
    }
    eb->Delete();

/*
    eb = vtkMultiBlockDataSet::New();
    this->multiBlock->SetBlock(SIDE_SETS_MBDS_ID, eb);
    this->multiBlock->GetMetaData(SIDE_SETS_MBDS_ID)->\
        Set(vtkCompositeDataSet::NAME(), SIDE_SETS_MBDS_NAME);
    eb->Delete();

    eb = vtkMultiBlockDataSet::New();
    this->multiBlock->SetBlock(NODE_SETS_MBDS_ID, eb);
    this->multiBlock->GetMetaData(NODE_SETS_MBDS_ID)->\
        Set(vtkCompositeDataSet::NAME(), NODE_SETS_MBDS_NAME);
    eb->Delete();
*/
}

void CatalystExodusMesh::CreateElementBlock(const char *elem_block_name,
    int elem_block_id, const std::string &elem_type, int nodes_per_elem,
        int num_elem, const int64_t *global_elem_ids, int *connectivity) {

    vtkVariant v((int)0);
    this->CreateElementBlockInternal(elem_block_name, elem_block_id,
        elem_type, nodes_per_elem, num_elem, v, global_elem_ids, connectivity);
}

void CatalystExodusMesh::CreateElementBlock(const char *elem_block_name,
    int elem_block_id, const std::string &elem_type, int nodes_per_elem,
        int num_elem, const int64_t *global_elem_ids, int64_t *connectivity) {

    vtkVariant v((int64_t)0);
    this->CreateElementBlockInternal(elem_block_name, elem_block_id,
        elem_type, nodes_per_elem, num_elem, v, global_elem_ids, connectivity);
}

void CatalystExodusMesh::CreateElementBlockInternal(
    const char * elem_block_name, int elem_block_id,
        const std::string &elem_type, int nodes_per_elem, int num_elem,
            vtkVariant &v, const int64_t *global_elem_ids,
                void *connectivity) {

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
        this->multiBlock->RemoveBlock(ELEMENT_BLOCK_MBDS_ID);
        return; // silently ignore empty element blocks
    }
    else {
        std::cout << "Unsupported element type: " << elemType.c_str();
        std::cout << "Unable to create element block " << elem_block_name;
        this->multiBlock->RemoveBlock(ELEMENT_BLOCK_MBDS_ID);
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

            int64_t conn_ind = this->GetArrayValue(v, connectivity, index_orig) - 1;
            this->ebmap[elem_block_id][conn_ind] = -1;
            index++;
        }
        if (point_count < nodes_per_elem) {
            for (int p = 0; p < (nodes_per_elem - point_count); p++) {
                index++;
            }
        }
    }

    vtkMultiBlockDataSet *eb =
        vtkMultiBlockDataSet::SafeDownCast(\
            this->multiBlock->GetBlock(ELEMENT_BLOCK_MBDS_ID));
    vtkUnstructuredGrid *ug =
        vtkUnstructuredGrid::SafeDownCast(eb->GetBlock(\
            this->ebidmap[elem_block_id]));
    vtkPoints *points = vtkPoints::New();
    double x[3];
    vtkIdType i = 0;
    for (std::map<int, int>::iterator ii = this->ebmap[elem_block_id].begin();
        ii != this->ebmap[elem_block_id].end(); ++ii) {

        this->global_points->GetPoint((*ii).first, x);
        (*ii).second = i;
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

            conn_ind = this->GetArrayValue(v, connectivity, index++) - 1;
            pts[cell_vertex_order[p]] = this->ebmap[elem_block_id][conn_ind];
        }
        if (point_count < nodes_per_elem) {
            for (int p = 0; p < (nodes_per_elem - point_count); p++) {
                int64_t conn_ind = this->GetArrayValue(v, connectivity, index) - 1;
                this->global_point_id_to_global_elem_id[conn_ind] = global_elem_ids[ug->GetNumberOfCells()];
                index++;
            }
        }

        this->global_elem_id_map[elem_block_id][global_elem_ids[ug->GetNumberOfCells()]] =
            ug->GetNumberOfCells();

        ug->InsertNextCell(vtk_type, point_count, pts);
        object_ids.push_back(this->ebidmap[elem_block_id]);
    }

    ug->SetPoints(points);
    eb->GetMetaData(this->ebidmap[elem_block_id])->\
        Set(vtkCompositeDataSet::NAME(), elem_block_name);
    points->Delete();

    std::vector<std::string> element_block_name;
    element_block_name.push_back("ElementBlockIds");
    vtkVariant vb((int)this->ebidmap[elem_block_id]);
    int bid = this->ebidmap[elem_block_id];
    this->CreateGlobalVariableInternal(element_block_name, eb,
        this->ebidmap[elem_block_id], vb, &bid);

    std::vector<std::string> component_names;
    component_names.push_back("ObjectId");
    this->CreateElementVariableInternal(component_names, eb, bid,
        vb, &object_ids[0]);
}

void CatalystExodusMesh::CreateNodeSet(const char *node_set_name,
    int node_set_id, int num_ids, const int *data) {

    vtkVariant v((int)0);
    this->CreateNodeSetVariant(node_set_name, node_set_id, num_ids, v, data);
}

void CatalystExodusMesh::CreateNodeSet(const char *node_set_name,
    int node_set_id, int num_ids, const int64_t *data) {

    vtkVariant v((int64_t)0);
    this->CreateNodeSetVariant(node_set_name, node_set_id, num_ids, v, data);
}

void CatalystExodusMesh::CreateNodeSetVariant(const char *node_set_name,
    int node_set_id, int num_ids, vtkVariant &v, const void *ids) {

    vtkMultiBlockDataSet *eb = vtkMultiBlockDataSet::SafeDownCast(\
        this->multiBlock->GetBlock(NODE_SETS_MBDS_ID));
    unsigned int bid = eb->GetNumberOfBlocks();
    vtkUnstructuredGrid *ug = vtkUnstructuredGrid::New();
    vtkPoints *points = vtkPoints::New();
    points->SetNumberOfPoints(num_ids);
    this->nsidmap[node_set_id] = bid;
    std::vector<int> object_ids;
    std::vector<int> global_element_id;

    for (int i = 0; i < num_ids; i++) {
        double x[3];
        vtkIdType ptids[1];
        int id = (int)this->GetArrayValue(v, ids, i) - 1;
        this->global_points->GetPoint(id, x);
        ptids[0] = (vtkIdType)i;
        points->SetPoint(i, x);
        ug->InsertNextCell(VTK_VERTEX, 1, ptids);
        this->nsmap[node_set_id][id] = i;
        object_ids.push_back(node_set_id);

        global_element_id.push_back(this->global_point_id_to_global_elem_id[id]);
    }

    ug->SetPoints(points);
    eb->SetBlock(bid, ug);
    eb->GetMetaData(bid)->Set(vtkCompositeDataSet::NAME(), node_set_name);
    ug->Delete();
    points->Delete();

    if (num_ids > 0) {
        std::vector<std::string> component_names;
        component_names.push_back("ObjectId");
        vtkVariant val((int)0);
        this->CreateElementVariableInternal(component_names, eb, bid,
            val, &object_ids[0]);

        component_names.clear();
        component_names.push_back("GlobalElementId");
        this->CreateElementVariableInternal(component_names, eb, bid,
            val, &global_element_id[0]);
    }
}

void CatalystExodusMesh::CreateSideSet(const char *ss_owner_name,
    int side_set_id, int num_ids, const int *element_ids,
        const int *face_ids) {
    /*NOTE: Jeff Mauldin JAM 2015Oct8
    CreateSideSet is called once for each block which the sideset
    spans, and the side_set_name for the side set is the ss_owner_name
    with additional characters to indicate which block we are doing.
    The current implementation of the sierra sideset construction
    creates a single independent sideset and collects all the
    nodes and elements from the side set from each block spanned by
    the sideset into that single sideset.  It needs to have the
    ss_owner_name, not the side_set_name, because that is the name
    in the input deck for the sideset for reference for things like
    extractblock.  It may become necessary at a later date to
    pass in both, but for now we
    are just passing in ss_owner_name to give us correct
    functionality while not changing the function interface*/

    vtkVariant v((int)0);
    this->CreateSideSetVariant(ss_owner_name, side_set_id, num_ids,
        v, element_ids, face_ids);
}

void CatalystExodusMesh::CreateSideSet(const char *ss_owner_name,
    int side_set_id, int num_ids, const int64_t *element_ids,
        const int64_t *face_ids) {

    vtkVariant v((int64_t)0);
    this->CreateSideSetVariant(ss_owner_name, side_set_id, num_ids,
        v, element_ids, face_ids);
}

void CatalystExodusMesh::CreateSideSetVariant(const char *ss_owner_name,
    int side_set_id, int num_ids, vtkVariant &v, const void *element_ids,
        const void *face_ids) {

    vtkMultiBlockDataSet *ssb = vtkMultiBlockDataSet::SafeDownCast(\
        this->multiBlock->GetBlock(SIDE_SETS_MBDS_ID));
    vtkMultiBlockDataSet *eb = vtkMultiBlockDataSet::SafeDownCast(\
        this->multiBlock->GetBlock(ELEMENT_BLOCK_MBDS_ID));

    // side set can span multiple blocks; we appear to get callback for same
    // side set once per block in which the side set exists
    // we need to track this situation and only construct one unstructured grid
    // per sideset, but put info from every block the sideset spans into the
    // unstructured grid (and the map of the point ids)
    int newSideSetFlag;
    vtkPoints * points = nullptr;
    vtkUnstructuredGrid *ug = nullptr;
    unsigned int bid = -1;

    Ve2mSideSetInfo * ssinfo = nullptr;
    std::map<int, Ve2mSideSetInfo *>::iterator ssinfomapiter;
    ssinfomapiter = this->ssinfomap.find(side_set_id);
    if (ssinfomapiter != ssinfomap.end()) {
        // we have already seen this side_set_id, so we add this side set rather
        // than constructing a new one
        newSideSetFlag = 0;
        ssinfo = ssinfomapiter->second;
        bid = ssinfo->bid;
        ug = vtkUnstructuredGrid::SafeDownCast(ssb->GetBlock(bid));
        points = ug->GetPoints();
    }
    else {
        // we have not seen this side_set_id, so we need to create it and add
        // it to the map of information about existing side sets
        newSideSetFlag = 1;
        bid = ssb->GetNumberOfBlocks();
        ug = vtkUnstructuredGrid::New();
        points = vtkPoints::New();
        ssinfo = new Ve2mSideSetInfo();
        ssinfomap[side_set_id] = ssinfo;
        ssinfo->bid = bid;
        this->ssidmap[side_set_id] = bid;
    }

    vtkIdType point_index = ssinfo->unique_points.size();

    for (int i = 0; i < num_ids; i++) {
        int e_id;
        int e_bid;
        int face_id;
        vtkUnstructuredGrid *eb_ug = 0;
        vtkCell *cell = 0;

        for (std::map<int, unsigned int>::iterator iter = this->ebidmap.begin();
            iter != this->ebidmap.end(); ++iter) {

            eb_ug = vtkUnstructuredGrid::SafeDownCast(\
                eb->GetBlock(this->ebidmap[iter->first]));
            e_id = this->GetArrayValue(v, element_ids, i);
            if (this->global_elem_id_map[iter->first].find(e_id) !=
                this->global_elem_id_map[iter->first].end()) {

                face_id = this->GetArrayValue(v, face_ids, i) - 1;
                e_id = this->global_elem_id_map[iter->first][e_id];
                e_bid = iter->first;
                cell = eb_ug->GetCell(e_id);
                break;
            }
        }

        if (cell != 0 && cell->GetCellType() != VTK_EMPTY_CELL) {
            if (cell->GetCellType() == VTK_HEXAHEDRON ||
                cell->GetCellType() == VTK_QUADRATIC_HEXAHEDRON ||
                cell->GetCellType() == VTK_TRIQUADRATIC_HEXAHEDRON) {

                face_id = HEXAHEDRON_FACE_MAP[face_id];
            }
            else if (cell->GetCellType() == VTK_WEDGE ||
                     cell->GetCellType() == VTK_QUADRATIC_WEDGE) {
                face_id = WEDGE_FACE_MAP[face_id];
            }
            else {
            }

            vtkCell *face = eb_ug->GetCell(e_id)->GetFace(face_id);
            if (!face) {
                face = eb_ug->GetCell(e_id)->GetEdge(face_id);
            }
            vtkPoints *face_points = face->GetPoints();
            double x[3];
            vtkIdType pts[face_points->GetNumberOfPoints()];

            for (int j = 0; j < face_points->GetNumberOfPoints(); j++) {
                int global_point_id = this->\
                    ebmap_reverse[e_bid][face->GetPointId(j)];
                if (ssinfo->unique_points.find(global_point_id) !=
                    ssinfo->unique_points.end()) {

                    pts[j] = ssinfo->unique_points[global_point_id];
                    if (pts[j] >= points->GetNumberOfPoints()) {
                        std::cerr << "ERROR!!!  BAD INDEX\n";
                    }
                }
                else {
                    if (point_index > points->GetNumberOfPoints()) {
                        std::cerr << "ERROR!!!  BAD INDEX\n";
                    }
                    ssinfo->unique_points[global_point_id] = point_index;
                    this->ssmap[side_set_id][global_point_id] = point_index;
                    pts[j] = point_index++;
                    face_points->GetPoint(j, x);
                    points->InsertNextPoint(x);
                }
            }
            ug->InsertNextCell(face->GetCellType(),
                face->GetNumberOfPoints(), pts);
            ssinfo->object_ids.push_back(side_set_id);
        }
    }

    if (newSideSetFlag) {
        ug->SetPoints(points);
        ssb->SetBlock(bid, ug);
        ssb->GetMetaData(bid)->Set(vtkCompositeDataSet::NAME(), ss_owner_name);
        ug->Delete();
        points->Delete();
    }

    if (num_ids > 0) {
        std::vector<std::string> component_names;
        component_names.push_back("ObjectId");
        vtkVariant val((int)0);
        this->CreateElementVariableInternal(component_names, ssb, bid, val,
            &(ssinfo->object_ids)[0]);
    }
}

void CatalystExodusMesh::CreateElementVariable(
    std::vector<std::string> &component_names, int elem_block_id,
        const double *data) {

    vtkVariant v((double)0.0);
    this->CreateElementVariableVariant(component_names,
        elem_block_id, v, data);
}

void CatalystExodusMesh::CreateElementVariable(
    std::vector<std::string> &component_names, int elem_block_id,
        const int *data) {

    vtkVariant v((int)0);
    this->CreateElementVariableVariant(component_names,
        elem_block_id, v, data);
}

void CatalystExodusMesh::CreateElementVariable(
    std::vector<std::string> &component_names, int elem_block_id,
        const int64_t *data) {

    vtkVariant v((int64_t)0);
    this->CreateElementVariableVariant(component_names,
        elem_block_id, v, data);
}

void CatalystExodusMesh::CreateElementVariableVariant(
    std::vector<std::string> &component_names, int elem_block_id, vtkVariant &v,
        const void *data) {

  vtkMultiBlockDataSet *eb =
      vtkMultiBlockDataSet::SafeDownCast(this->multiBlock->\
          GetBlock(ELEMENT_BLOCK_MBDS_ID));
  this->CreateElementVariableInternal(component_names, eb,
      this->ebidmap[elem_block_id], v, data);
}

void CatalystExodusMesh::CreateElementVariableInternal(
    std::vector<std::string> &component_names, vtkMultiBlockDataSet *eb,
        unsigned int bid, vtkVariant &v, const void *data) {

    vtkUnstructuredGrid *ug = vtkUnstructuredGrid::SafeDownCast(\
        eb->GetBlock(bid));

    int number_data_components = 1;
    std::vector<std::string> prefix_name;
    this->ContainsVector(component_names, prefix_name);
    std::vector<std::string> component_names_buffer = component_names;
    if (prefix_name.size() == 1) {
        number_data_components = component_names.size();
        component_names_buffer = prefix_name;
    }

    vtkFieldData* cell_data = ug->GetCellData();
    std::vector<vtkDataArray *> data_arrays;
    for (std::vector<std::string>::iterator it = component_names_buffer.begin();
         it != component_names_buffer.end(); ++it) {
        vtkDataArray *da = cell_data->GetArray((*it).c_str());
        if (da) {
            cell_data->RemoveArray((*it).c_str());
        }
        vtkDataArray *arr = vtkDataArray::CreateDataArray(v.GetType());
        arr->SetName((*it).c_str());
        arr->SetNumberOfComponents(number_data_components);
        arr->SetNumberOfTuples(ug->GetNumberOfCells());
        cell_data->AddArray(arr);
        data_arrays.push_back(arr);
        arr->Delete();
    }

    const double *dp = static_cast<const double *>(data);
    int index = 0;
    for (int i = 0; i < ug->GetNumberOfCells(); i++) {
        for (std::vector<vtkDataArray *>::iterator it = data_arrays.begin();
            it != data_arrays.end(); ++it) {
            for (int c = 0; c < number_data_components; c++) {
                (*it)->SetComponent(i, c,
                    this->GetArrayValue(v, data, index++));
            }
        }
    }
}

void CatalystExodusMesh::CreateNodalVariable(
    std::vector<std::string> &component_names, const double *data) {

    vtkVariant v((double)0.0);
    this->CreateNodalVariableVariant(component_names, v, data);
}

void CatalystExodusMesh::CreateNodalVariable(
    std::vector<std::string> &component_names, const int *data) {

    vtkVariant v((int)0);
    this->CreateNodalVariableVariant(component_names, v, data);
}

void CatalystExodusMesh::CreateNodalVariable(
    std::vector<std::string> &component_names, const int64_t *data) {

    vtkVariant v((int64_t)0);
    this->CreateNodalVariableVariant(component_names, v, data);
}

void CatalystExodusMesh::CreateNodalVariableVariant(
    std::vector<std::string> &component_names, vtkVariant &v,
        const void *data) {

    vtkMultiBlockDataSet *eb =
        vtkMultiBlockDataSet::SafeDownCast(\
            this->multiBlock->GetBlock(ELEMENT_BLOCK_MBDS_ID));

    this->CreateNodalVariableInternal(component_names, eb, this->ebidmap,
        this->ebmap, v, data);

/*
    eb = vtkMultiBlockDataSet::SafeDownCast(\
        this->multiBlock->GetBlock(NODE_SETS_MBDS_ID));

    this->CreateNodalVariableInternal(component_names, eb, this->nsidmap,
        this->nsmap, v, data);

    eb = vtkMultiBlockDataSet::SafeDownCast(\
        this->multiBlock->GetBlock(SIDE_SETS_MBDS_ID));

    this->CreateNodalVariableInternal(component_names, eb, this->ssidmap,
        this->ssmap, v, data);
*/
}

void CatalystExodusMesh::CreateNodalVariableInternal(
    std::vector<std::string> &component_names, vtkMultiBlockDataSet *eb,
        std::map<int, unsigned int> &id_map, std::map<int,
            std::map<int, int>> &point_map, vtkVariant &v, const void *data) {

    int number_data_components = 1;
    std::vector<std::string> prefix_name;
    bool displace_nodes = false;
    this->ContainsVector(component_names, prefix_name);
    std::vector<std::string> component_names_buffer = component_names;
    if (prefix_name.size() == 1) {
        number_data_components = component_names.size();
        component_names_buffer = prefix_name;
        if (this->ApplyDisplacementsON() && (number_data_components <= 3) &&
            (prefix_name[0].length() >= 3)) {
            if ((prefix_name[0].substr(0, 3) == "DIS") ||\
                (prefix_name[0].substr(0, 3) == "dis")) {
                displace_nodes = true;
            }
        }
    }

    for (std::map<int, unsigned int>::iterator iter = id_map.begin();
        iter != id_map.end(); ++iter) {

        vtkUnstructuredGrid *ug = vtkUnstructuredGrid::SafeDownCast(eb->GetBlock(iter->second));
        if (!ug) {
            continue;
        }
        if (ug->GetNumberOfCells() == 0) {
            continue;
        }
        vtkFieldData* point_data = ug->GetPointData();
        std::vector<vtkDataArray *> data_arrays;

        for (std::vector<std::string>::iterator it = component_names_buffer.begin();
             it != component_names_buffer.end(); ++it) {

            vtkDataArray *da = point_data->GetArray((*it).c_str());
            if (da) {
                point_data->RemoveArray((*it).c_str());
            }
            vtkDataArray *arr = vtkDataArray::CreateDataArray(v.GetType());
            arr->SetName((*it).c_str());
            arr->SetNumberOfComponents(number_data_components);
            arr->SetNumberOfTuples(ug->GetPoints()->GetNumberOfPoints());
            point_data->AddArray(arr);
            data_arrays.push_back(arr);
            arr->Delete();
        }

        int index = 0;
        for (int i = 0; i < this->num_global_points; i++) {
            std::map<int, int>::iterator mit = point_map[iter->first].find(i);
            for (std::vector<vtkDataArray *>::iterator it = data_arrays.begin();
                it != data_arrays.end(); ++it) {

                for (int c = 0; c < number_data_components; c++) {
                    if (mit != point_map[iter->first].end()) {
                        (*it)->SetComponent(point_map[iter->first][i], c,
                            this->GetArrayValue(v, data, index++));
                    }
                    else {
                        index++;
                    }
                }

                if (displace_nodes) {
                    if (mit != point_map[iter->first].end()) {
                        vtkPoints *points = ug->GetPoints();
                        double x[3];
                        this->global_points->GetPoint(i, x);
                        for (int c = 0; c < number_data_components; c++) {
                            x[c] += (*it)->GetComponent(\
                                point_map[iter->first][i], c);
                        }
                        points->SetPoint(point_map[iter->first][i], x);
                    }
                }
            }
        }
    }
}

void CatalystExodusMesh::ReleaseMemory() {

    vtkMultiBlockDataSet *eb =
        vtkMultiBlockDataSet::SafeDownCast(\
            this->multiBlock->GetBlock(ELEMENT_BLOCK_MBDS_ID));
    this->ReleaseMemoryInternal(eb);

/*
    eb = vtkMultiBlockDataSet::SafeDownCast(\
        this->multiBlock->GetBlock(NODE_SETS_MBDS_ID));
    this->ReleaseMemoryInternal(eb);

    eb = vtkMultiBlockDataSet::SafeDownCast(\
        this->multiBlock->GetBlock(SIDE_SETS_MBDS_ID));
    this->ReleaseMemoryInternal(eb);
*/

    this->ebmap_reverse.clear();
    this->global_elem_id_map.clear();
    this->global_point_id_to_global_elem_id.clear();

    this->catManager->WriteToLogFile(catalystPipelineInfo);
}

void CatalystExodusMesh::ReleaseMemoryInternal(vtkMultiBlockDataSet *eb) {

    vtkCompositeDataIterator *iter = eb->NewIterator();
    for (iter->InitTraversal(); !iter->IsDoneWithTraversal();
        iter->GoToNextItem()) {

        vtkUnstructuredGrid *ug = vtkUnstructuredGrid::SafeDownCast(\
            iter->GetCurrentDataObject());
        vtkFieldData * data = ug->GetCellData();
        vtkDataArray * global_element_ids = nullptr;
        vtkDataArray * object_id = nullptr;
        while (data->GetNumberOfArrays() > 0) {
            if (std::string(data->GetArray(0)->GetName()) ==
                "GlobalElementId") {

                global_element_ids = data->GetArray(0);
                global_element_ids->Register(0);
            }
            if (std::string(data->GetArray(0)->GetName()) == "ObjectId") {
                object_id = data->GetArray(0);
                object_id->Register(0);
            }
            data->RemoveArray(data->GetArray(0)->GetName());
        }

        if (global_element_ids) {
            data->AddArray(global_element_ids);
            global_element_ids->Delete();
        }

        if (object_id) {
            data->AddArray(object_id);
            object_id->Delete();
        }

        data = ug->GetPointData();
        vtkDataArray *global_node_ids = 0;
        while (data->GetNumberOfArrays() > 0) {
            if (std::string(data->GetArray(0)->GetName()) == "GlobalNodeId") {
                global_node_ids = data->GetArray(0);
                global_node_ids->Register(0);
            }
            data->RemoveArray(data->GetArray(0)->GetName());
        }

        if (global_node_ids) {
            data->AddArray(global_node_ids);
            global_node_ids->Delete();
        }
    }
    iter->Delete();
}

void CatalystExodusMesh::logMemoryUsageAndTakeTimerReading() {
    this->catManager->logMemoryUsageAndTakeTimerReading(
        catalystPipelineInfo);
}

void CatalystExodusMesh::Delete() {
    this->catManager->DeletePipeline(catalystPipelineInfo);
}

void CatalystExodusMesh::ContainsVector(std::vector<std::string> &component_names,
                                        std::vector<std::string> &prefix_name) {

  if (component_names.size() == 3) {
    if ((*component_names[0].rbegin() == 'X' || *component_names[0].rbegin() == 'x') &&
        (*component_names[1].rbegin() == 'Y' || *component_names[1].rbegin() == 'y') &&
        (*component_names[2].rbegin() == 'Z' || *component_names[2].rbegin() == 'z')) {
      prefix_name.push_back(component_names[0].substr(0, component_names[0].size() - 1));
    }
  }
  else if (component_names.size() == 2) {
    if ((*component_names[0].rbegin() == 'X' || *component_names[0].rbegin() == 'x') &&
        (*component_names[1].rbegin() == 'Y' || *component_names[1].rbegin() == 'y')) {
      prefix_name.push_back(component_names[0].substr(0, component_names[0].size() - 1));
    }
  }
  if (!this->UnderscoreVectorsON() && prefix_name.size() == 1) {
    if (*prefix_name[0].rbegin() == '_') {
      prefix_name[0] = prefix_name[0].substr(0, prefix_name[0].size() - 1);
    }
  }
}

double CatalystExodusMesh::GetArrayValue(vtkVariant &v, const void *data, int index) {
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

} // namespace Iovs_exodus
