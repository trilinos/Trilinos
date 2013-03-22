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

enum Topology {
    Shell4 = 4,
    Hex8   = 8
};

inline std::string getTopologyName(Topology topology)
{
    switch(topology)
    {
        case Shell4:
            return std::string("shell4");
        case Hex8:
            return std::string("hex8");
    }
    throw std::exception();
}

struct ExodusData {
    std::vector<double> coordinates;
    const std::vector< std::vector<int> > elementBlockConnectivity;
    const std::vector<int> globalNumberOfElementsInBlock;
    const std::vector<int> localNumberOfElementsInBlock;
    const std::vector<Topology> blockTopologicalData;

    const int globalNumberOfNodes;

    const std::vector<int> globalIdsOfLocalElements;
    const std::vector<int> globalIdsOfLocalNodes;

    std::vector<SharedNode> *sharedNodes;
    ExodusData(const std::vector<double> coords, const std::vector< std::vector<int> > elemBlockConnectivity,
               const std::vector<int> globalNumOfElemsInBlock, const std::vector<int> localNumOfElemsInBlock,
               const std::vector<Topology> blockTopoData, int globalNumNodes, const std::vector<int> globalIdsOfLocalElems,
               const std::vector<int> globalIdsLocalNodes)
    : coordinates(coords), elementBlockConnectivity(elemBlockConnectivity),
      globalNumberOfElementsInBlock(globalNumOfElemsInBlock), localNumberOfElementsInBlock(localNumOfElemsInBlock),
      blockTopologicalData(blockTopoData), globalNumberOfNodes(globalNumNodes), globalIdsOfLocalElements(globalIdsOfLocalElems),
      globalIdsOfLocalNodes(globalIdsLocalNodes), sharedNodes(0)
    {

    }
};

struct DashSurfaceData
{
    const std::vector<double> &coordinates;
    const std::vector<int> &surfaceAConnectivity;
    const std::vector<int> &surfaceBConnectivity;

    int globalNumberOfNodes;
    int globalNumberOfElements;

    int globalNumberOfElementsSurface1;
    int globalNumberOfElementsSurface2;

    std::vector<int> globalIdsOfLocalElements;
    std::vector<int> globalIdsOfLocalNodes;

    std::vector<SharedNode> *sharedNodes;

    DashSurfaceData(const std::vector<double> &coords, const std::vector<int> &connectivity1, const std::vector<int> &connectivity2)
    : coordinates(coords), surfaceAConnectivity(connectivity1), surfaceBConnectivity(connectivity2),
      sharedNodes(0)
    {
        this->setSerialDefaults();
    }

  private:

    void setSerialDefaults()
    {
        globalNumberOfNodes = coordinates.size()/SPATIAL_DIMENSION;

        globalNumberOfElementsSurface1 = surfaceBConnectivity.size()/NUM_NODES_PER_QUAD_FACE;
        globalNumberOfElementsSurface2 = surfaceAConnectivity.size()/NUM_NODES_PER_QUAD_FACE;
        globalNumberOfElements = globalNumberOfElementsSurface1 + globalNumberOfElementsSurface2;

        globalIdsOfLocalElements.resize(globalNumberOfElements);
        globalIdsOfLocalNodes.resize(globalNumberOfNodes);

        for (size_t i=0; i<globalIdsOfLocalElements.size();i++)
        {
            globalIdsOfLocalElements[i] = i+1;
        }

        for (size_t i=0; i<globalIdsOfLocalNodes.size();i++)
        {
            globalIdsOfLocalNodes[i] = i+1;
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

class ExodusMesh : public GeneratedMesh
{
public:

    explicit ExodusMesh(const ExodusData &exodusData);

    virtual ~ExodusMesh() { }

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
    int64_t mGlobalNumberOfElements;
    int64_t mLocalNumberOfElements;

    const ExodusData & mExodusData;
    std::vector<int64_t> mElementOffsetForBlock;
};

}

#endif
