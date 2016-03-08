
#include "Zoltan2_TaskMapping.hpp"
#include <Zoltan2_TestHelpers.hpp>
#include "Tpetra_MultiVector_decl.hpp"
#include <string>
#include "Teuchos_XMLParameterListHelpers.hpp"
typedef int part_t;
int main(int argc, char *argv[]){


  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Teuchos::Comm<int> > tcomm = Teuchos::DefaultComm<int>::getComm();
  int nx = 2, ny = 2, nz = 2;
  for ( int i = 1 ; i < argc ; ++i ) {
    if ( 0 == strcasecmp( argv[i] , "NX" ) ) {
      nx = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "NY" ) ) {
      ny = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "NZ" ) ) {
      nz = atoi( argv[++i] );
    }
    else{
      std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl ;
      return 1;
    }
  }


  try{
    int rank = tcomm->getRank();
    part_t numProcs = tcomm->getSize();


    Zoltan2::MachineRepresentation<zscalar_t, part_t> mach(*tcomm);

    int procDim = mach.getMachineDim();
    zscalar_t **procCoordinates = NULL;
    mach.getAllMachineCoordinatesView(procCoordinates);




    int coordDim = 3;
    part_t numParts = nx*ny*nz;
    zscalar_t **partCenters = NULL;

    partCenters = new zscalar_t * [coordDim];
    for(int i = 0; i < coordDim; ++i){
      partCenters[i] = new zscalar_t[numParts];
    }


    part_t *task_communication_xadj_ = new part_t [numParts+1];
    part_t *task_communication_adj_ = new part_t [numParts * 6];

    int prevNCount = 0;
    task_communication_xadj_[0] = 0;
    for (part_t i = 0; i < numParts; ++i) {
      int x = i % nx;
      int y = (i / (nx)) % ny;
      int z = (i / (nx)) / ny;
      partCenters[0][i] = x;
      partCenters[1][i] = y;
      partCenters[2][i] = z;

      if (x > 0){
        task_communication_adj_[prevNCount++] = i - 1;
      }
      if (x < nx - 1){
        task_communication_adj_[prevNCount++] = i + 1;
      }
      if (y > 0){
        task_communication_adj_[prevNCount++] = i - nx;
      }
      if (y < ny - 1){
        task_communication_adj_[prevNCount++] = i + nx;
      }
      if (z > 0){
        task_communication_adj_[prevNCount++] = i - nx * ny;
      }
      if (z < nz - 1){
        task_communication_adj_[prevNCount++] = i + nx * ny;
      }
      task_communication_xadj_[i+1] = prevNCount;
    }


    part_t *proc_to_task_xadj_ = new part_t[numProcs+1];
    part_t *proc_to_task_adj_ = new part_t[numParts];
    part_t *partArray = NULL;
    int partArraysize = -1;

    part_t *machineDimensions = NULL;
    //machineDimensions = hopper;
    Zoltan2::coordinateTaskMapperInterface<part_t, zscalar_t, zscalar_t>(
        tcomm,
        procDim,
        numProcs,
        procCoordinates,

        coordDim,
        numParts,
        partCenters,

        task_communication_xadj_,
        task_communication_adj_,
        NULL,

        proc_to_task_xadj_, /*output*/
        proc_to_task_adj_, /*output*/

        partArraysize,
        partArray,
        machineDimensions
    );

    if (tcomm->getRank() == 0){
      cout << "PASS" << endl;
    }

    if (tcomm->getRank() == 0){
      for (part_t i = 0; i < numProcs; ++i){
        std::cout << "\nProc i:" << i << " ";
        for (part_t j = proc_to_task_xadj_[i]; j < proc_to_task_xadj_[i+1]; ++j){
          std::cout << " " << proc_to_task_adj_[j];
        }
      }
      std::cout << std::endl;
    }
    delete [] proc_to_task_xadj_;
    delete [] proc_to_task_adj_;
    delete [] task_communication_xadj_;
    delete [] task_communication_adj_;

    for (int i = 0; i < coordDim; i++) delete [] partCenters[i];
    delete [] partCenters;

  }
  catch(std::string &s){
    cerr << s << endl;
  }

  catch(char * s){
    cerr << s << endl;
  }
  catch(char const * s){
    cerr << s << endl;
  }
}

