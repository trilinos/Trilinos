#ifndef GEOMETRYFACTORY_HPP
#define GEOMETRYFACTORY_HPP

#include <string>
#include "GeometryKernel.hpp"
#include "MeshGeometry.hpp"
class PerceptMesh;


class GeometryFactory
{
public:
    GeometryFactory(GeometryKernel* kernel, MeshGeometry* geometry);
    ~GeometryFactory();

    bool read_file(const std::string& filename);
    void create_evaluators(PerceptMesh* mesh_data);

protected:
    GeometryKernel* geomKernel;
    MeshGeometry* geomDatabase;

    std::vector<GeometryHandle> geometryHandles;

    PointSet& get_nodes_in_sideset(PerceptMesh* mesh_data, int sideset);
};

#endif // GEOMETRYFACTORY_HPP
