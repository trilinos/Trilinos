#ifndef MUELU_ZOLTAN_DEF_HPP
#define MUELU_ZOLTAN_DEF_HPP

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)

#include <Teuchos_Utils.hpp>

#include "MueLu_Zoltan_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  ZoltanInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  ZoltanInterface(GO numPartitions, RCP<const FactoryBase> AFact) : numPartitions_(numPartitions), AFact_(AFact)
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  void ZoltanInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  DeclareInput(Level & level) const
  {
    level.DeclareInput("A", AFact_.get()); // use DeclareInput instead of Request!!!
  } //DeclareInput()

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  void ZoltanInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  SetNumberOfPartitions(GO const numPartitions)
  {
    numPartitions_ = numPartitions;
  } //SetNumberOfPartitions()

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  void ZoltanInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  Build(Level &level) const
  {
    RCP<Operator> A = level.Get< RCP<Operator> >("A",AFact_.get());
    // Tell Zoltan what kind of local/global IDs we will use.
    // In our case, each GID is two ints and there are no local ids.
    // One can skip this step if the IDs are just single ints.
    RCP<const Teuchos::MpiComm<int> > mpiComm = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(A->getRowMap()->getComm());
    float zoltanVersion_;
    Zoltan_Initialize(0, NULL, &zoltanVersion_);
    //TODO define zoltanComm_ as a subcommunicator?!;
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > zoltanComm_ = mpiComm->getRawMpiComm();
    RCP<Zoltan> zoltanObj_ = rcp( new Zoltan( (*zoltanComm_)() ) );  //extract the underlying MPI_Comm handle and create a Zoltan object
    if (zoltanObj_==Teuchos::null) throw(Exceptions::RuntimeError("MueLu::Zoltan : Unable to create Zoltan data structure"));
    int rv;
    if ((rv=zoltanObj_->Set_Param("num_gid_entries", "1")) != ZOLTAN_OK )
      throw(Exceptions::RuntimeError("MueLu::Zoltan::Setup : setting parameter 'num_gid_entries' returned error code " + Teuchos::toString(rv)));
    if ( (rv=zoltanObj_->Set_Param("num_lid_entries", "0") ) != ZOLTAN_OK )
      throw(Exceptions::RuntimeError("MueLu::Zoltan::Setup : setting parameter 'num_lid_entries' returned error code " + Teuchos::toString(rv)));
    if ( (rv=zoltanObj_->Set_Param("obj_weight_dim", "1") ) != ZOLTAN_OK )
      throw(Exceptions::RuntimeError("MueLu::Zoltan::Setup : setting parameter 'obj_weight_dim' returned error code " + Teuchos::toString(rv)));

    std::stringstream ss;
    ss << numPartitions_;
    zoltanObj_->Set_Param("num_global_partitions",ss.str());

    RCP<MultiVector> XYZ = level.Get< RCP<MultiVector> >("coordinates",MueLu::NoFactory::get());
    size_t problemDimension_ = XYZ->getNumVectors();

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
    RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition;
    if (newDecomp) {
      decomposition = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(rowMap,false); //Don't bother initializing, as this will just be overwritten.
      ArrayRCP<GO> decompEntries = decomposition->getDataNonConst(0);
      for (typename ArrayRCP<GO>::iterator i = decompEntries.begin(); i != decompEntries.end(); ++i)
        *i = mypid;
      for (int i=0; i< num_exported; ++i)
        decompEntries[ rowMap->getLocalElement(export_gids[i]) ] = export_to_part[i];
    } else {
      //Running on one processor, so decomposition is the trivial one, all zeros.
      decomposition = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(rowMap,true);
    }
    level.Set<RCP<Xpetra::Vector<GO,LO,GO,NO> > >("partition",decomposition);
    level.Set<GO>("number of partitions",numPartitions_);

    zoltanObj_->LB_Free_Part(&import_gids, &import_lids, &import_procs, &import_to_part);
    zoltanObj_->LB_Free_Part(&export_gids, &export_lids, &export_procs, &export_to_part);

  } //Build()

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  int ZoltanInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  GetLocalNumberOfRows(void *data, int *ierr)
  {
    if (data == NULL) {
      *ierr = ZOLTAN_FATAL;
      return -1;
    }
    *ierr = ZOLTAN_OK;
    //TODO is there a safer way to cast?
    //Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> *A = (Operator*) data;
    Operator *A = (Operator*) data;
    return A->getRowMap()->getNodeNumElements();
  } //GetLocalNumberOfRows()

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  void ZoltanInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  GetLocalNumberOfNonzeros(void *data, int NumGidEntries, int NumLidEntries, ZOLTAN_ID_PTR gids,
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

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  int ZoltanInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  GetProblemDimension(void *data, int *ierr)
  {
    //TODO is there a safer way to cast?
    int dim = *((int*)data);
    *ierr = ZOLTAN_OK; /* set error flag */
    return(dim);
  } //GetProblemDimension

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  void ZoltanInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  GetProblemGeometry(void *data, int numGIDEntries, int numLIDEntries, int numObjectIDs, 
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

} //namespace MueLu

#endif //if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)

#endif // MUELU_ZOLTAN_DEF_HPP
