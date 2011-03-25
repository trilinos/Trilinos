#ifndef MESHGEOMETRY_HPP
#define MESHGEOMETRY_HPP

#include "GeometryKernel.hpp"

typedef std::vector<double> PointSet;
typedef int GeometryHandle;

struct GeometryEvaluator
{
    GeometryHandle geometry;
    PointSet mesh;
};

class MeshGeometry
{
public:
    MeshGeometry(GeometryKernel* geom);
    ~MeshGeometry();

    void add_evaluator(GeometryEvaluator* evaluator);
    void add_evaluators(std::vector<GeometryEvaluator*> evaluators);

    void snap_points_to_geometry();

protected:
    std::vector<GeometryEvaluator*> geomEvaluators;
    GeometryKernel* geomKernel;
};

#endif // MESHGEOMETRY_HPP
