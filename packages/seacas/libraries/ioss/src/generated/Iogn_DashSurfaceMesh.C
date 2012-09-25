#include <generated/Iogn_DashSurfaceMesh.h>

#include <vector>

namespace Iogn
{

int64_t DashSurfaceMesh::node_count() const
{
    return mCoordinates.size()/SPATIAL_DIMENSION;
}

int64_t DashSurfaceMesh::node_count_proc() const
{
    return node_count();
}

int64_t DashSurfaceMesh::element_count() const
{
    return (mQuadSurface1.size()+mQuadSurface2.size())/NUM_NODES_PER_QUAD_FACE;
}

int64_t DashSurfaceMesh::element_count(int64_t block_number) const
{
    if(block_number == 1)
    {
        return mQuadSurface1.size()/NUM_NODES_PER_QUAD_FACE;
    }
    else if(block_number == 2)
    {
        return mQuadSurface2.size()/NUM_NODES_PER_QUAD_FACE;
    }
    throw std::exception();

    return INVALID;
}

int64_t DashSurfaceMesh::block_count() const
{
    return 2;
}

int64_t DashSurfaceMesh::nodeset_count() const
{
    return 0;
}

int64_t DashSurfaceMesh::sideset_count() const
{
    return 2;
}

int64_t DashSurfaceMesh::element_count_proc() const
{
    return element_count();
}

int64_t DashSurfaceMesh::element_count_proc(int64_t block_number) const
{
    return element_count(block_number);
}

int64_t DashSurfaceMesh::nodeset_node_count_proc(int64_t id) const
{
    return 0;
}
int64_t DashSurfaceMesh::sideset_side_count_proc(int64_t id) const
{
    return element_count(id);
}
int64_t DashSurfaceMesh::communication_node_count_proc() const
{
    return 0;
}

void DashSurfaceMesh::coordinates(double *coord) const
{
    std::copy(mCoordinates.begin(),mCoordinates.end(), coord);
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
            std::copy(mQuadSurface1.begin(),mQuadSurface1.end(), connect);
            return;
        case 2:
            std::copy(mQuadSurface2.begin(),mQuadSurface2.end(), connect);
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
    size_t numElementsInSurface1 = mQuadSurface1.size()/NUM_NODES_PER_QUAD_FACE;
    switch(setId)
    {
        case 1:
            for(size_t i=0; i<numElementsInSurface1; ++i)
            {
                elem_sides.push_back(i+1);
                elem_sides.push_back(0);
            }
            return;
        case 2:
            for(size_t i=0; i<mQuadSurface2.size()/NUM_NODES_PER_QUAD_FACE; ++i)
            {
                elem_sides.push_back(numElementsInSurface1+i+1);
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
    return;
}

void DashSurfaceMesh::node_map(IntVector &map)
{
    map.resize(node_count());
    for(int i = 0; i < node_count(); i++)
    {
        map[i] = i + 1;
    }
}

void DashSurfaceMesh::node_map(MapVector &map)
{
    map.resize(node_count());
    for(int i = 0; i < node_count(); i++)
    {
        map[i] = i + 1;
    }
}

void DashSurfaceMesh::element_map(int block_number, IntVector &map) const
{
    size_t numElementsInSurface1 = mQuadSurface1.size() / NUM_NODES_PER_QUAD_FACE;
    size_t numElementsInSurface2 = mQuadSurface2.size() / NUM_NODES_PER_QUAD_FACE;
    switch(block_number)
    {
        case 1:
            for(size_t i = 0; i < numElementsInSurface1; ++i)
            {
                map[i] = i + 1;
            }
            return;
        case 2:
            for(size_t i = 0; i < numElementsInSurface2; ++i)
            {
                map[numElementsInSurface1 + i] = numElementsInSurface1 + i + 1;
            }
            return;
        default:
            throw std::exception();
    }
}

void DashSurfaceMesh::element_map(int64_t block_number, MapVector &map) const
{
    size_t numElementsInSurface1 = mQuadSurface1.size() / NUM_NODES_PER_QUAD_FACE;
    size_t numElementsInSurface2 = mQuadSurface2.size() / NUM_NODES_PER_QUAD_FACE;
    switch(block_number)
    {
        case 1:
            for(size_t i = 0; i < numElementsInSurface1; ++i)
            {
                map[i] = i + 1;
            }
            return;
        case 2:
            for(size_t i = 0; i < numElementsInSurface2; ++i)
            {
                map[numElementsInSurface1 + i] = numElementsInSurface1 + i + 1;
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
        map[i] = i + 1;
    }
}

void DashSurfaceMesh::element_map(IntVector &map) const
{
    int count = element_count_proc();
    map.resize(count);

    for(int i = 0; i < count; i++)
    {
        map[i] = i + 1;
    }
}

}
