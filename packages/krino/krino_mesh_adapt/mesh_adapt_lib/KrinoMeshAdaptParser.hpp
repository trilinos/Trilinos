#ifndef KrinoMeshAdaptParser_hpp
#define KrinoMeshAdaptParser_hpp

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <KrinoMeshAdaptInputData.hpp>
#include <stk_util/parallel/Parallel.hpp>
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace krino
{

class MeshAdaptParser
{
public:
    const MeshAdaptInputData& get_input_data() const { return mInputData; }
    bool read_command_line(int argc, char *argv[], stk::ParallelMachine comm);

private:
    MeshAdaptInputData mInputData;
};

}


#endif
