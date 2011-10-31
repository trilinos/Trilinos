#ifndef GEOMETRYKERNELSTUPID_HPP
#define GEOMETRYKERNELSTUPID_HPP

#include "GeometryKernel.hpp"
#include <fstream>

class GeometryKernelStupid : public GeometryKernel
{
public:
    GeometryKernelStupid() : GeometryKernel() {}
    virtual ~GeometryKernelStupid() {}

    virtual bool read_file(const std::string& file_name,
                           std::vector<GeometryHandle>& geometry_entities)
    {
        std::ifstream file (file_name.c_str());
        if (!file.is_open())
            return false;
        int i=1;
        while (!file.eof())
        {
            geometry_entities.push_back(i);
            char string[256];
            file.getline(string, 255);
            geometryAttribute[i]=string;
            i++;
        }
        return true;
    }

    virtual std::string get_attribute(GeometryHandle geom)
    {
      return geometryAttribute[geom];
    }

    virtual void snap_to(KernelPoint& point, GeometryHandle geom)
    { }

    virtual bool is_curve(GeometryHandle geom) const
    {
      return true;
    }

    virtual bool is_surface(GeometryHandle geom) const
    {
      return true;
    }

private:
    std::map<int, std::string> geometryAttribute;
};

#endif // GEOMETRYKERNELSTUPID_HPP
