#ifndef KrinoMeshAdaptInputData_hpp
#define KrinoMeshAdaptInputData_hpp

#include <KrinoMeshAdaptAlgorithmParameters.hpp>
#include <string>

namespace krino
{

struct MeshAdaptInputData
{
    static constexpr const char * mDefaultString{"None"};
    std::string meshIn{mDefaultString};
    std::string meshOut{mDefaultString};
    MeshAdaptAlgorithmParameters algorithmParams;

    MeshAdaptInputData() {}
    MeshAdaptInputData(const std::string& iMeshIn,
                           const std::string& iMeshOut) :
        meshIn(iMeshIn),
        meshOut(iMeshOut)
    {}
};

}

#endif
