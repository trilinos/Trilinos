#ifndef IOSS_Iogn_DashSurfaceMesh_h
#define IOSS_Iogn_DashSurfaceMesh_h

#include <generated/Iogn_GeneratedMesh.h>

namespace Iogn
{

enum
{
    INVALID                 = -1,
    NUMBER_OF_SURFACES      =  2,
    SPATIAL_DIMENSION       =  3,
    NUM_NODES_PER_QUAD_FACE =  4
};

struct SharedNode
{
    SharedNode() : nodeId(-1), procId(-1) {}
    int nodeId;
    int procId;
};

struct DashSurfaceData
{
    const std::vector<double> &coordinates;
    const std::vector<int> &surface1Connectivity;
    const std::vector<int> &surface2Connectivity;

    int globalNumberOfNodes;
    int globalNumberOfElements;

    int globalNumberOfElementsSurface1;
    int globalNumberOfElementsSurface2;

    std::vector<int> elementGlobalIds;
    std::vector<int> nodeGlobalIds;


    std::vector<SharedNode> *sharedNodes;

    DashSurfaceData(const std::vector<double> &coords, const std::vector<int> &connectivity1, const std::vector<int> &connectivity2)
    : coordinates(coords), surface1Connectivity(connectivity1), surface2Connectivity(connectivity2),
      sharedNodes(0)
    {
        this->setSerialDefaults();
    }

  private:

    void setSerialDefaults()
    {
        globalNumberOfNodes = coordinates.size()/SPATIAL_DIMENSION;

        globalNumberOfElementsSurface1 = surface1Connectivity.size()/NUM_NODES_PER_QUAD_FACE;
        globalNumberOfElementsSurface2 = surface2Connectivity.size()/NUM_NODES_PER_QUAD_FACE;
        globalNumberOfElements = globalNumberOfElementsSurface1 + globalNumberOfElementsSurface2;

        elementGlobalIds.resize(globalNumberOfElements);
        nodeGlobalIds.resize(globalNumberOfNodes);

        for (size_t i=0; i<elementGlobalIds.size();i++)
        {
            elementGlobalIds[i] = i+1;
        }

        for (size_t i=0; i<nodeGlobalIds.size();i++)
        {
            nodeGlobalIds[i] = i+1;
        }
    }
};

class DashSurfaceMesh : public GeneratedMesh
{
public:

    explicit DashSurfaceMesh(DashSurfaceData &dashSurfaceData) : mDashSurfaceData(dashSurfaceData) { }

    virtual ~DashSurfaceMesh() { }

    virtual int64_t node_count() const;
    virtual int64_t node_count_proc() const;

    virtual int64_t element_count() const;
    virtual int64_t element_count(int64_t block_number) const;
    virtual int64_t element_count_proc() const;
    virtual int64_t element_count_proc(int64_t block_number) const;

    virtual int64_t block_count() const;

    virtual int64_t nodeset_count() const;
    virtual int64_t nodeset_node_count_proc(int64_t id) const;

    virtual int64_t sideset_count() const;
    virtual int64_t sideset_side_count_proc(int64_t id) const;

    virtual int64_t communication_node_count_proc() const;

    virtual void coordinates(double *coord) const;
    virtual void coordinates(std::vector<double> &coord) const;
    virtual void coordinates(int component, std::vector<double> &xyz) const;
    virtual void coordinates(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z) const;

    virtual void connectivity(int64_t block_number, int* connect) const;

    virtual std::pair<std::string, int> topology_type(int64_t block_number) const;

    virtual void sideset_elem_sides(int64_t setId, std::vector<int64_t> &elem_sides) const;

    virtual void nodeset_nodes(int64_t nset_id, std::vector<int64_t> &nodes) const;

    virtual void node_communication_map(std::vector<int64_t> &map, std::vector<int> &proc);

    virtual void node_map(std::vector<int> &map);
    virtual void node_map(std::vector<int64_t> &map);

    virtual void element_map(int block_number, std::vector<int> &map) const;
    virtual void element_map(int64_t block_number, std::vector<int64_t> &map) const;
    virtual void element_map(std::vector<int64_t> &map) const;
    virtual void element_map(std::vector<int> &map) const;

private:
    DashSurfaceData & mDashSurfaceData;

};

}

#endif
