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

private:
    int networkDim;
    int numProcs;

    pcoord_t **procCoords;
    RCP<const Comm<int> > comm;

public:
    /*! \brief Constructor MachineRepresentation Class
     *  \param comm_ Communication object.
     */
    MachineRepresentation(RCP<const Comm<int> > comm_):
        networkDim(0), numProcs(comm_->getSize()), procCoords(0), 
        comm(comm_)
    {
        // WIll need this constructor to be specific to RAAMP (MD).
        // Will need a default constructor using, e.g., GeometricGenerator
        // or nothing at all, for when RAAMP is not available as TPL.
        //
        // (AG) In addition, need to be able to run without special
        // privileges in system (e.g., on hopper).
        // Notes:  For now, all cores connected to same NIC will get the
        // same coordinates; later, we could add extra coordinate dimensions
        // to represent nodes or dies (using hwloc info through RAAMP
        // data object).

        // (MD) will modify mapping test to use machine representation
        // #ifdef HAVE_ZOLTAN2_OVIS

        // Call initializer for RAAMP data object (AG)

        //get network dimension.
        //TODO change.
        // Call RAAMP Data Object to get the network dimension (AG)
        networkDim = 3;

        //allocate memory for processor coordinates.
        procCoords = new pcoord_t *[networkDim];
        for (int i = 0; i < networkDim; ++i){
            procCoords[i] = new pcoord_t [numProcs];
            memset (procCoords[i], 0, sizeof(pcoord_t) * numProcs);
        }
        //obtain the coordinate of the processor.
        this->getMyCoordinate(/*pcoord_t &xyz[networkDim]*/);
        // copy xyz into appropriate spot in procCoords. (MD)

        //reduceAll the coordinates of each processor.
        this->gatherMachineCoordinates();
    }


    /*! \brief Constructor MachineRepresentation Class
     *  \param comm_ Communication object.
     */
    MachineRepresentation(RCP<Comm<int> > comm_):
        networkDim(0), numProcs(comm_->getSize()), procCoords(0), comm(comm_){
        // WIll need this constructor to be specific to RAAMP (MD).
        // Will need a default constructor using, e.g., GeometricGenerator
        // or nothing at all, for when RAAMP is not available as TPL.
        //
        // (AG) In addition, need to be able to run without special
        // privileges in system (e.g., on hopper).  
        // Notes:  For now, all cores connected to same NIC will get the
        // same coordinates; later, we could add extra coordinate dimensions
        // to represent nodes or dies (using hwloc info through RAAMP
        // data object).

        // (MD) will modify mapping test to use machine representation
        // #ifdef HAVE_ZOLTAN2_OVIS

        // Call initializer for RAAMP data object (AG)

        //get network dimension.
        //TODO change.
        // Call RAAMP Data Object to get the network dimension (AG)
        networkDim = 3;

        //allocate memory for processor coordinates.
        procCoords = new pcoord_t *[networkDim];
        for (int i = 0; i < networkDim; ++i){
            procCoords[i] = new pcoord_t [numProcs];
            memset (procCoords[i], 0, sizeof(pcoord_t) * numProcs);
        }
        //obtain the coordinate of the processor.
        this->getMyCoordinate(/*pcoord_t &xyz[networkDim]*/);
        // copy xyz into appropriate spot in procCoords. (MD)

        //reduceAll the coordinates of each processor.
        this->gatherMachineCoordinates();
    }


    /*! \brief getMyCoordinate function
     *  stores the coordinate of the current processor in procCoords[*][rank]
     */
    void getMyCoordinate(/* pcoord_t &xyz[networkDim]*/){

        // Call RAAMP system to get coordinates and store in xyz (MD)
        // What is the RAAMP call?  (AG)
        // AG will return a view (pointer) to RAAMP's data.
        // We will copy it into xyz.

        // The code below may be good for the default constructor, perhaps,
        // but it should copy the data into xyz instead of the procCoords.
        int myRank = comm->getRank();

        int slice = int (pow( double(numProcs), double(1.0 / networkDim)) + 0.5 );

        int m = myRank;
        for (int i = 0; i < networkDim; ++i){
            procCoords[i][myRank] = m / int(pow(slice, double(networkDim - i - 1)));
            m = m % int(pow(double(slice), double(networkDim - i - 1)));
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
        // Free/release THE RAAMP Data Object.
        // Deinitialize/finalize/whatever (AG)
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
