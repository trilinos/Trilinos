#include <generated/Iogn_DashSurfaceMesh.h>

#include <vector>

namespace Iogn
{

int64_t DashSurfaceMesh::node_count() const
{
    return mDashSurfaceData.globalNumberOfNodes;
}

int64_t DashSurfaceMesh::node_count_proc() const
{
    return mDashSurfaceData.coordinates.size()/SPATIAL_DIMENSION;
}

int64_t DashSurfaceMesh::element_count() const
{
    return mDashSurfaceData.globalNumberOfElements;
}

int64_t DashSurfaceMesh::element_count(int64_t surfaceNumber) const
{
    if(surfaceNumber == 1)
    {
        return mDashSurfaceData.globalNumberOfElementsSurface1;
    }
    if(surfaceNumber == 2)
    {
        return mDashSurfaceData.globalNumberOfElementsSurface2;
    }
    throw std::exception();

    return INVALID;
}

int64_t DashSurfaceMesh::block_count() const
{
    return NUMBER_OF_SURFACES;
}

int64_t DashSurfaceMesh::nodeset_count() const
{
    return 0;
}

int64_t DashSurfaceMesh::sideset_count() const
{
    return NUMBER_OF_SURFACES;
}

int64_t DashSurfaceMesh::element_count_proc() const
{
    return (mDashSurfaceData.surface1Connectivity.size() + mDashSurfaceData.surface2Connectivity.size())/NUM_NODES_PER_QUAD_FACE;
}

int64_t DashSurfaceMesh::element_count_proc(int64_t block_number) const
{
    if(block_number == 1)
    {
        return mDashSurfaceData.surface1Connectivity.size()/NUM_NODES_PER_QUAD_FACE;
    }
    else if(block_number == 2)
    {
        return mDashSurfaceData.surface2Connectivity.size()/NUM_NODES_PER_QUAD_FACE;
    }
    throw std::exception();

    return INVALID;
}

int64_t DashSurfaceMesh::nodeset_node_count_proc(int64_t id) const
{
    return 0;
}

int64_t DashSurfaceMesh::sideset_side_count_proc(int64_t id) const
{
    return element_count_proc(id);
}

int64_t DashSurfaceMesh::communication_node_count_proc() const
{
    if ( mDashSurfaceData.sharedNodes )
        return mDashSurfaceData.sharedNodes->size();
    else
       return 0;
}

void DashSurfaceMesh::coordinates(double *coord) const
{
    std::copy(mDashSurfaceData.coordinates.begin(),mDashSurfaceData.coordinates.end(), coord);
}

void DashSurfaceMesh::coordinates(std::vector<double> &coord) const
{
    throw std::exception();
}

void DashSurfaceMesh::coordinates(int component, std::vector<double> &xyz) const
{
    throw std::exception();
}

void DashSurfaceMesh::coordinates(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z) const
{
    throw std::exception();
}

void DashSurfaceMesh::connectivity(int64_t block_number, int* connect) const
{
    switch(block_number)
    {
        case 1:
            std::copy(mDashSurfaceData.surface1Connectivity.begin(),mDashSurfaceData.surface1Connectivity.end(), connect);
            return;
        case 2:
            std::copy(mDashSurfaceData.surface2Connectivity.begin(),mDashSurfaceData.surface2Connectivity.end(), connect);
            return;
        default:
            throw std::exception();
    }
}

std::pair<std::string, int> DashSurfaceMesh::topology_type(int64_t block_number) const
{
    const int numNodesPerElement = 4;
    return std::make_pair(std::string("shell4"), numNodesPerElement);
}

void DashSurfaceMesh::sideset_elem_sides(int64_t setId, Int64Vector &elem_sides) const
{
    elem_sides.clear();
    size_t numElementsInSurface1 = element_count_proc(1);
    size_t numElementsInSurface2 = element_count_proc(2);
    switch(setId)
    {
        case 1:
            for(size_t i=0; i<numElementsInSurface1; ++i)
            {
                elem_sides.push_back(mDashSurfaceData.elementGlobalIds[i]);
                elem_sides.push_back(0);
            }
            return;
        case 2:
            for(size_t i=0; i<numElementsInSurface2; ++i)
            {
                elem_sides.push_back(mDashSurfaceData.elementGlobalIds[numElementsInSurface1+i]);
                elem_sides.push_back(0);
            }
            return;
        default:
            throw std::exception();
    }
}

void DashSurfaceMesh::nodeset_nodes(int64_t nset_id, Int64Vector &nodes) const
{
    return;
}

void DashSurfaceMesh::node_communication_map(MapVector &map, std::vector<int> &proc)
{
    if ( mDashSurfaceData.sharedNodes == 0 ) return;

    for (unsigned int i=0;i<mDashSurfaceData.sharedNodes->size();i++)
    {
        map[i] = (*mDashSurfaceData.sharedNodes)[i].nodeId;
        proc[i] = (*mDashSurfaceData.sharedNodes)[i].procId;
    }
    return;
}

void DashSurfaceMesh::node_map(IntVector &map)
{
    int size = node_count_proc();
    map.resize(size);

    for(int i = 0; i < size; i++)
    {
        map[i] = mDashSurfaceData.nodeGlobalIds[i];
    }
}

void DashSurfaceMesh::node_map(MapVector &map)
{
    int size = node_count_proc();
    map.resize(size);

    for(int i = 0; i < size; i++)
    {
        map[i] = mDashSurfaceData.nodeGlobalIds[i];
    }
}

void DashSurfaceMesh::element_map(int block_number, IntVector &map) const
{
    size_t numElementsInSurface1 = element_count_proc(1);
    size_t numElementsInSurface2 = element_count_proc(2);
    switch(block_number)
    {
        case 1:
            for(size_t i = 0; i < numElementsInSurface1; ++i)
            {
                map[i] = mDashSurfaceData.elementGlobalIds[i];
            }
            return;
        case 2:
            for(size_t i = 0; i < numElementsInSurface2; ++i)
            {
                map[numElementsInSurface1 + i] = mDashSurfaceData.elementGlobalIds[numElementsInSurface1 + i];
            }
            return;
        default:
            throw std::exception();
    }
}

void DashSurfaceMesh::element_map(int64_t block_number, MapVector &map) const
{
    size_t numElementsInSurface1 = element_count_proc(1);
    size_t numElementsInSurface2 = element_count_proc(2);
    switch(block_number)
    {
        case 1:
            for(size_t i = 0; i < numElementsInSurface1; ++i)
            {
                map[i] = mDashSurfaceData.elementGlobalIds[i];
            }
            return;
        case 2:
            for(size_t i = 0; i < numElementsInSurface2; ++i)
            {
                map[numElementsInSurface1 + i] = mDashSurfaceData.elementGlobalIds[numElementsInSurface1 + i];
            }
            return;
        default:
            throw std::exception();
    }
}

void DashSurfaceMesh::element_map(MapVector &map) const
{
    int count = element_count_proc();
    map.resize(count);

    for(int i = 0; i < count; i++)
    {
        map[i] = mDashSurfaceData.elementGlobalIds[i];
    }
}

void DashSurfaceMesh::element_map(IntVector &map) const
{
    int count = element_count_proc();
    map.resize(count);

    for(int i = 0; i < count; i++)
    {
        map[i] = mDashSurfaceData.elementGlobalIds[i];
    }
}

}
