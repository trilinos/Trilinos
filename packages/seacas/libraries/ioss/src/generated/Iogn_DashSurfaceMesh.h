#ifndef IOSS_Iogn_DashSurfaceMesh_h
#define IOSS_Iogn_DashSurfaceMesh_h

#include <generated/Iogn_GeneratedMesh.h>

namespace Iogn {

class DashSurfaceMesh : public GeneratedMesh
{
public:
    enum {
        INVALID                 = -1,
        SPATIAL_DIMENSION       =  3,
        NUM_NODES_PER_QUAD_FACE =  4
    };

    explicit DashSurfaceMesh(const std::vector<double>    &coords,
                        const std::vector< int64_t > &quadSurface1,
                        const std::vector< int64_t > &quadSurface2)
    : mCoordinates(coords),
      mQuadSurface1(quadSurface1),
      mQuadSurface2(quadSurface2)
    {}

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
    const std::vector<double>    &mCoordinates;
    const std::vector< int64_t > &mQuadSurface1;
    const std::vector< int64_t > &mQuadSurface2;

};

}

#endif
