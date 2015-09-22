#ifndef GEOMETRYKERNEL_OPENNURBS_HPP
#define GEOMETRYKERNEL_OPENNURBS_HPP

//#include <stk_percept/mesh/geometry/opennurbs/opennurbs.h>
#include <opennurbs.h>
#include "GeometryKernel.hpp"

class GeometryKernelOpenNURBS : public GeometryKernel
{
public:
    GeometryKernelOpenNURBS();
    virtual ~GeometryKernelOpenNURBS();

    virtual bool read_file(const std::string& file_name,
                           std::vector<GeometryHandle>& geometry_entities);

    virtual std::string get_attribute(GeometryHandle geom);

    virtual void snap_to(KernelPoint& point, GeometryHandle geom,
                    double *converged_tolerance = NULL,
                    double *uvw_computed = NULL,
                    double *uvw_hint = NULL);
    virtual void normal_at(KernelPoint& point, GeometryHandle geom, std::vector<double>& normal);

    virtual bool is_curve(GeometryHandle geom) const;

    virtual bool is_surface(GeometryHandle geom) const;

private:
    ONX_Model onModel;
};

#endif // GEOMETRYKERNEL_OPENNURBS_HPP
