#ifndef _ZOLTAN2_MACHINEREPRESENTATION_HPP_
#define _ZOLTAN2_MACHINEREPRESENTATION_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>
namespace Zoltan2{

/*! \brief MachineRepresentation Class
 * Finds the coordinate of the processor.
 * Used to find the processor coordinates or graph.
 */
template <typename pcoord_t>
class MachineRepresentation{

public:
    int networkDim;
    int numProcs;

    pcoord_t **procCoords;
    RCP<Comm<int> > comm;


    /*! \brief Constructor MachineRepresentation Class
     *  \param comm_ Communication object.
     */
    MachineRepresentation(RCP<Comm<int> > comm_):
        networkDim(0), numProcs(comm_->getSize()), procCoords(0), comm(comm_){

        //get network dimension.
        //TODO change.
        networkDim = 3;
        //allocate memory for processor coordinates.
        procCoords = new pcoord_t *[networkDim];
        for (int i = 0; i < networkDim; ++i){
            procCoords[i] = new pcoord_t [numProcs];
            memset (procCoords[i], 0, sizeof(pcoord_t) * numProcs);
        }
        //obtain the coordinate of the processor.
        this->getMyCoordinate();
        //reduceAll the coordinates of each processor.
        this->gatherMachineCoordinates();
    }


    /*! \brief getMyCoordinate function
     *  stores the coordinate of the current processor in procCoords[*][rank]
     */
    void getMyCoordinate(){
        int myRank = comm->getRank();

        int slice = pow( double(numProcs), double(1.0 / networkDim)) + 0.5 ;

        int m = myRank;
        for (int i = 0; i < networkDim; ++i){
            procCoords[i][myRank] = m / int(pow(slice, networkDim - i - 1));
            m = m % int(pow(slice, networkDim - i - 1));
        }
    }

    /*! \brief gatherMachineCoordinates function
     *  reduces and stores all machine coordinates.
     */
    void gatherMachineCoordinates(){
        pcoord_t *tmpVect = new pcoord_t [numProcs];

        for (int i = 0; i < networkDim; ++i){
            reduceAll<int, pcoord_t>(
                    *comm,
                    Teuchos::REDUCE_SUM,
                    numProcs,
                    procCoords[i],
                    tmpVect);
            pcoord_t *tmp = tmpVect;
            tmpVect = procCoords[i];
            procCoords[i] = tmp;
        }
        delete [] tmpVect;
    }

    /*! \brief destructor of the class
     * free memory in procCoords.
     */
    ~MachineRepresentation() {
        for (int i = 0; i < networkDim; ++i){
            delete [] procCoords[i];
        }
        delete []procCoords;
    }

    /*! \brief getProcDim function
     * returns the dimension of the physical processor layout.
     */
    int getProcDim() const{
        return networkDim;
    }

    /*! \brief getProcDim function
     * returns the coordinates of processors in two dimensional array.
     */
    pcoord_t** getProcCoords() const{
        return procCoords;
    }

    /*! \brief getNumProcs function
     * returns the number of processors.
     */
    int getNumProcs() const{
        return numProcs;
    }

};
}
#endif
