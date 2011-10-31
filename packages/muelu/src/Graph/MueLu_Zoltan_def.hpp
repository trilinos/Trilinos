#ifndef MUELU_ZOLTAN_HPP
#define MUELU_ZOLTAN_HPP

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)

#include "zoltan_cpp.h"

#include "Xpetra_Operator.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_VectorFactory.hpp"
#include "Xpetra_Map.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

/*!
  @class ZoltanInterface class
  @brief Interface to Zoltan library.

  Currently, this interface provides access only to the RCB partitioning in Zoltan.
*/

template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType,
          class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
class ZoltanInterface : public BaseClass {

#include "MueLu_UseShortNames.hpp"

  public:

   //! @name Constructors/Destructors
   //@{

   //! Constructor
   ZoltanInterface(RCP<const Teuchos::Comm<int> > const &comm, RCP<const FactoryBase> AFact = Teuchos::null) : AFact_(AFact) {

     Zoltan_Initialize(0, NULL, &zoltanVersion_);
     //TODO define zoltanComm_ as a subcommunicator?!;
     comm_ = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
     zoltanComm_ = comm_->getRawMpiComm();
     zoltanObj_ = rcp( new Zoltan( (*zoltanComm_)() ) );  //extract the underlying MPI_Comm handle and create a Zoltan object
     if (zoltanObj_==Teuchos::null) throw(Exceptions::RuntimeError("MueLu::Zoltan : Unable to create Zoltan data structure"));
   }

   //! Destructor
   virtual ~ZoltanInterface() {}
   //@}

   //! @name Input
   //@{
   void DeclareInput(Level & level) const {
     level.Request("A", AFact_.get());
     
   }
   //@}

   //! @name Set / Get methods.
   //@{
   void SetNumberOfPartitions(GO const numPartitions) {
     std::stringstream ss;
     ss << numPartitions;
     zoltanObj_->Set_Param("num_global_partitions",ss.str());
   }
   //@}

   //! @name Build methods.
   //@{
   void Build(Level &level) {
     // Tell Zoltan what kind of local/global IDs we will use.
     // In our case, each GID is two ints and there are no local ids.
     // One can skip this step if the IDs are just single ints.
     int rv;
     if ((rv=zoltanObj_->Set_Param("num_gid_entries", "1")) != ZOLTAN_OK )
       throw(Exceptions::RuntimeError("MueLu::Zoltan::Setup : setting parameter 'num_gid_entries' returned error code " + toString(rv)));
     if ( (rv=zoltanObj_->Set_Param("num_lid_entries", "0") ) != ZOLTAN_OK )
       throw(Exceptions::RuntimeError("MueLu::Zoltan::Setup : setting parameter 'num_lid_entries' returned error code " + toString(rv)));
     if ( (rv=zoltanObj_->Set_Param("obj_weight_dim", "1") ) != ZOLTAN_OK )
       throw(Exceptions::RuntimeError("MueLu::Zoltan::Setup : setting parameter 'obj_weight_dim' returned error code " + toString(rv)));

     RCP<Operator> A = level.Get< RCP<Operator> >("A",AFact_.get());
     RCP<MultiVector> XYZ = level.Get< RCP<MultiVector> >("coordinates",MueLu::NoFactory::get());
     problemDimension_ = XYZ->getNumVectors();

     zoltanObj_->Set_Num_Obj_Fn(GetLocalNumberOfRows,(void *) &*A);
     zoltanObj_->Set_Obj_List_Fn(GetLocalNumberOfNonzeros,(void *) &*A);
     zoltanObj_->Set_Num_Geom_Fn(GetProblemDimension, (void *) &problemDimension_);
     zoltanObj_->Set_Geom_Multi_Fn(GetProblemGeometry, (void *) &*XYZ);

     // Data pointers that Zoltan requires.
     ZOLTAN_ID_PTR import_gids = NULL;  // Global nums of objs to be imported   
     ZOLTAN_ID_PTR import_lids = NULL;  // Local indices to objs to be imported 
     int   *import_procs = NULL;        // Proc IDs of procs owning objs to be imported.
     int   *import_to_part = NULL;      // Partition #s to which imported objs should be assigned.
     ZOLTAN_ID_PTR export_gids = NULL;  // Global nums of objs to be exported   
     ZOLTAN_ID_PTR export_lids = NULL;  // local indices to objs to be exported 
     int   *export_procs = NULL;        // Proc IDs of destination procs for objs to be exported.
     int   *export_to_part = NULL;      // Partition #s for objs to be exported.
     int   num_imported;                // Number of objs to be imported.          
     int   num_exported;                // Number of objs to be exported.          
     int   newDecomp;                   // Flag indicating whether the decomposition has changed 
     int   num_gid_entries;             // Number of array entries in a global ID.  
     int   num_lid_entries;

     rv = zoltanObj_->LB_Partition(newDecomp, num_gid_entries, num_lid_entries,
                  num_imported, import_gids, import_lids, import_procs, import_to_part,
                  num_exported, export_gids, export_lids, export_procs, export_to_part);
     if (rv == ZOLTAN_FATAL) {
       throw(Exceptions::RuntimeError("Zoltan::LB_Partition() returned error code"));
     }

     //TODO check that A's row map is 1-1.  Zoltan requires this.
     RCP<const Map> rowMap = A->getRowMap();
     int mypid = rowMap->getComm()->getRank();
     RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(rowMap,false);
     ArrayRCP<GO> decompEntries = decomposition->getDataNonConst(0);
     if (newDecomp) {
       //TODO Right now, the new partition is numbered 0 ... n.  Will Zoltan let us used different partitioning
       //numbers, e.g., the processor ID?  I can imagine cases where we'd like migrate to a nonconsecutive subset
       //of processors.
       //Answer -- no.   So, we'll have to look at make decisions here about which PIDs should be active.
       for (typename ArrayRCP<GO>::iterator i = decompEntries.begin(); i != decompEntries.end(); ++i)
         *i = mypid;
       for (int i=0; i< num_exported; ++i)
         decompEntries[ rowMap->getLocalElement(export_gids[i]) ] = export_to_part[i];
     }
     level.Set<RCP<Xpetra::Vector<GO,LO,GO,NO> > >("partition",decomposition);

     zoltanObj_->LB_Free_Part(&import_gids, &import_lids, &import_procs, &import_to_part);
     zoltanObj_->LB_Free_Part(&export_gids, &export_lids, &export_procs, &export_to_part);

   } //Build()

   //@}

   //! @name Query methods (really functions) required by Zoltan.
   //@{

  /*! Callback function that returns the local number of objects. Required by Zoltan.

     In this case, the number of objects is the number of local rows.
     
     @param data (in) void pointer to an Xpetra::Operator.
     @param ierr (out) error code.
  */
  static int GetLocalNumberOfRows(void *data, int *ierr) {
    if (data == NULL) {
      *ierr = ZOLTAN_FATAL;
      return -1;
    }
    *ierr = ZOLTAN_OK;
    //TODO is there a safer way to cast?
    //Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> *A = (Operator*) data;
    Operator *A = (Operator*) data;
    return A->getRowMap()->getNodeNumElements();
  }

  /*! Callback function that returns the local number of nonzeros in the matrix. Required by Zoltan.

    FIXME: Note that this will not work properly for non-point matrices.

     @param data (in) void pointer to an Xpetra::Operator
     @param weights (out) array whose <tt>i</tt><sup>th</sup> entry is the number of nonzeros in local row \c i.
     @param ierr (out) error code
  */
  static void GetLocalNumberOfNonzeros(void *data, int NumGidEntries, int NumLidEntries, ZOLTAN_ID_PTR gids,
                                ZOLTAN_ID_PTR lids, int wgtDim, float *weights, int *ierr)
  {
    if (data == NULL || NumGidEntries < 1) {
      *ierr = ZOLTAN_FATAL;
      return;
    } else {
      *ierr = ZOLTAN_OK;
    }

    //TODO is there a safer way to cast?
    Operator *A = (Operator*) data;
    RCP<const Map> map = A->getRowMap();
    Teuchos::ArrayView<const LO> cols;
    Teuchos::ArrayView<const SC> vals;
    for (size_t i=0; i<map->getNodeNumElements(); ++i) {
      gids[i] = (ZOLTAN_ID_TYPE) map->getGlobalElement(i);
      A->getLocalRowView(i,cols,vals);
      weights[i] = cols.size();
    }

  } //GetLocalNumberOfNonzeros()

  /*! Callback function that returns the problem dimension. Required by Zoltan.

     @param data (in) void pointer to integer dimension
     @param ierr (out) error code
  */
  static int GetProblemDimension(void *data, int *ierr) {
    //TODO is there a safer way to cast?
    int dim = *((int*)data);
    *ierr = ZOLTAN_OK; /* set error flag */
    return(dim);
  }


  /*! Callback function that returns the problem dimension. Required by Zoltan.

     @param data (in) void pointer to Xpetra::MultiVector.
     @param coordinates (out) array of double coordinates, arranged like so: [x1 y1 z1 x2 y2 z2 ...].
     @param ierr (out) error code

     TODO -- should I return a view of the coordinates instead of copying them?
  */
  static void GetProblemGeometry(void *data, int numGIDEntries, int numLIDEntries, int numObjectIDs, 
                          ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int dim, double *coordinates, int *ierr)
  {
    if (data == NULL) {
      *ierr = ZOLTAN_FATAL;
      return;
    }

   //TODO is there a safer way to cast?
    MultiVector *XYZ = (MultiVector*) data;
    if (dim != (int) XYZ->getNumVectors()) {
      //FIXME I'm assuming dim should be 1,2,or 3 coming in?!
      *ierr = ZOLTAN_FATAL;
      return;
    }

    ArrayRCP<ArrayRCP<const SC> > XYZdata(dim);
    for (int i=0; i<dim; ++i) XYZdata[i] = XYZ->getData(i);
    for (size_t i=0; i<XYZ->getLocalLength(); ++i) {
      for (int j=0; j<dim; ++j) {
        coordinates[i*dim+j] = (double) XYZdata[j][i];
      }
    }

    *ierr = ZOLTAN_OK;

  } //GetProblemGeometry

  //@}

/*
   void Build() {

  if (Zoltan_LB_Partition(&*zoltanObj_, &new_decomp, &num_gid_entries, &num_lid_entries,
               &num_imported, &import_gids,
               &import_lids, &import_procs, &import_to_part,
               &num_exported, &export_gids,
               &export_lids, &export_procs, &export_to_part) == ZOLTAN_FATAL){
    printf("fatal(8)  error returned from Zoltan_LB_Partition()\n");
    return 0;
  }

  //cleanup
  Zoltan_LB_Free_Part(&import_gids, &import_lids,
                      &import_procs, &import_to_part);
  Zoltan_LB_Free_Part(&export_gids, &export_lids,
                      &export_procs, &export_to_part);

   } //Build
*/

  private:

    RCP<const Teuchos::MpiComm<int> > comm_;
    RCP<const FactoryBase> AFact_;
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > zoltanComm_;
    RCP<Zoltan> zoltanObj_;
    float zoltanVersion_;
    size_t problemDimension_;

};  //class ZoltanInterface

} //namespace MueLu

#endif //if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)

#define MUELU_ZOLTANINTERFACE_SHORT
#endif //ifndef MUELU_ZOLTAN_HPP
