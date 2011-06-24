#ifndef MESHGEOMETRY_HPP
#define MESHGEOMETRY_HPP

#include <stk_percept/PerceptMesh.hpp>
#include "GeometryKernel.hpp"

typedef std::vector<double> PointSet;
typedef int GeometryHandle;
using namespace stk;
using namespace mesh;
using namespace percept;

struct GeometryEvaluator
{
    GeometryEvaluator(Part* part) : mMesh(Selector(*part)) {}
    GeometryHandle mGeometry;
    Selector mMesh;
};

class MeshGeometry
{
public:
    MeshGeometry(GeometryKernel* geom);
    ~MeshGeometry();

    void add_evaluator(GeometryEvaluator* evaluator);
    void add_evaluators(std::vector<GeometryEvaluator*> evaluators);

    void snap_points_to_geometry(PerceptMesh* mesh_data);

    const std::vector<GeometryEvaluator*>& getGeomEvaluators();

protected:
    std::vector<GeometryEvaluator*> geomEvaluators;
    GeometryKernel* geomKernel;
};

#endif // MESHGEOMETRY_HPP
