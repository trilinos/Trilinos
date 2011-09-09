#ifndef GEOMETRYKERNEL_HPP
#define GEOMETRYKERNEL_HPP

#include <vector>
#include <string>

typedef double* KernelPoint;
typedef int GeometryHandle;

class GeometryKernel
{
public:
    GeometryKernel() {}
    virtual ~GeometryKernel() {};

    virtual bool read_file(const std::string& file_name,
                           std::vector<GeometryHandle>& geometry_entities ) = 0;

    virtual std::string get_attribute(GeometryHandle geom) = 0;

    virtual void snap_to(KernelPoint& point, GeometryHandle geom) = 0;

    virtual bool is_curve(GeometryHandle geom) const = 0;

    virtual bool is_surface(GeometryHandle geom) const = 0;

};

#endif // GEOMETRYKERNEL_HPP
