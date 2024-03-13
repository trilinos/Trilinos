#ifndef KRINOMESHADAPTALGORITHMPARAMETERS_HPP_
#define KRINOMESHADAPTALGORITHMPARAMETERS_HPP_

namespace krino
{

struct MeshAdaptAlgorithmParameters
{
    bool force64Bit{false};
    bool assert32Bit{false};
    bool autoCompose{true};
    unsigned numUniformRefinementLevels{1};
};

}

#endif /* KRINOMESHADAPTALGORITHMPARAMETERS_HPP_ */
