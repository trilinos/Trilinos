#ifndef _ZOLTAN2_MACHINEREPRESENTATION_HPP_
#define _ZOLTAN2_MACHINEREPRESENTATION_HPP_

namespace Zoltan2{
template <typename pcoord_t>
class MachineRepresentation{

    int networkDim;
    int numProcs;
    pcoord_t **procCoords;

public:
    MachineRepresentation(int networkDim_, int numProcs_, pcoord_t **procCoords_){
        this->networkDim = networkDim_;
        this->numProcs = numProcs_;
        this->procCoords = procCoords_;
    }

    MachineRepresentation():
        networkDim(0), numProcs(0), procCoords(0){

    }
    ~MachineRepresentation(){}

    int getProcDim(){
        return networkDim;
    }
    pcoord_t** getProcCoords(){
        return procCoords;
    }
    int getNumProcs(){
        return numProcs;
    }

};
}
#endif
