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
    virtual ~GeometryKernel() = 0;

    virtual bool read_file(const std::string& file_name,
                           std::vector<GeometryHandle>& geometry_entities) = 0;

    virtual int get_attribute(const std::string name, GeometryHandle geom) = 0;

    virtual KernelPoint snap_to(KernelPoint point, GeometryHandle geom) = 0;
};

#endif // GEOMETRYKERNEL_HPP
