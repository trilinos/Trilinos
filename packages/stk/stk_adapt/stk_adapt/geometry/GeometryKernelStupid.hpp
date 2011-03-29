#ifndef GEOMETRYKERNELSTUPID_HPP
#define GEOMETRYKERNELSTUPID_HPP

#include "GeometryKernel.hpp"

class GeometryKernelStupid : public GeometryKernel
{
public:
    GeometryKernelStupid() : GeometryKernel() {}
    virtual ~GeometryKernelStupid() {}

    virtual bool read_file(const std::string& file_name,
                           std::vector<GeometryHandle>& geometry_entities)
    {for (int i=1; i<=6; i++) geometry_entities.push_back(i);
     return true;}

    virtual int get_attribute(const std::string name, GeometryHandle geom)
    {return -1;}

    virtual Point snap_to(KernelPoint point, GeometryHandle geom)
    {return point;}
};

#endif // GEOMETRYKERNELSTUPID_HPP
