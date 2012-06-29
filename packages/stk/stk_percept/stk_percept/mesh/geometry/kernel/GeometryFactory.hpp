#ifndef GEOMETRYFACTORY_HPP
#define GEOMETRYFACTORY_HPP

#include <string>
#include "GeometryKernel.hpp"
#include "MeshGeometry.hpp"
#include <stk_percept/PerceptMesh.hpp>
using namespace stk;
using namespace mesh;
using namespace percept;

class GeometryFactory
{
public:
    GeometryFactory(GeometryKernel* kernel, MeshGeometry* geometry);
    ~GeometryFactory();

    bool read_file(const std::string& filename, PerceptMesh* mesh);

protected:
    GeometryKernel* geomKernel;
    MeshGeometry* geomDatabase;
};

#endif // GEOMETRYFACTORY_HPP
