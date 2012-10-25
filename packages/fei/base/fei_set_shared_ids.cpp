
#include <fei_set_shared_ids.hpp>
#include <fei_CommMap.hpp>
#include <fei_CommUtils.hpp>
#include <fei_TemplateUtils.hpp>
#include <snl_fei_RecordCollection.hpp>
#include <fei_SharedIDs.hpp>
#include <fei_LinearDecomposition.hpp>

namespace fei {

void copy_into_shared_ids(const fei::CommMap<int>::Type& procs_to_ids_and_sharing_procs,
                          const snl_fei::RecordCollection& records,
                          fei::SharedIDs<int>& sharedIDs)
{
  //Given a CommMap object which maps procs to packed vectors of ids and sharing
  //procs, and a RecordCollection of records that appear in the local subdomain:
  //Loop through the CommMap and for each id that exists in 'records', copy that
  //id along with the associated sharing procs into the sharedIDs object.

  fei::CommMap<int>::Type::const_iterator
    ip_iter = procs_to_ids_and_sharing_procs.begin(),
    ip_end = procs_to_ids_and_sharing_procs.end();

  for(; ip_iter != ip_end; ++ip_iter) {
    const std::vector<int>& ids_and_procs = ip_iter->second;
    size_t vsize = ids_and_procs.size();
    size_t offset = 0;
    while(offset < vsize) {
      int id = ids_and_procs[offset++];
      int num_procs = ids_and_procs[offset++];
      if (records.getRecordWithID(id) != NULL) {
        sharedIDs.addSharedID(id, num_procs, &ids_and_procs[offset]);
      }
      offset += num_procs;
    }
  }
}

void copy_remotelyowned_ids_into_CommMap(int myProc,
                                         const fei::LinearDecomposition<int>& lindecomp,
                                         const snl_fei::RecordCollection& records,
                                         fei::CommMap<int>::Type& procs_to_shared_ids)
{
  for(size_t i=0; i<records.getNumRecords(); ++i) {
    int ID = records.getRecordWithLocalID(i)->getID();
    int proc = lindecomp.which_proc(ID);
    if (proc != myProc) {
      addItemsToCommMap(proc, 1, &ID, procs_to_shared_ids);
    }
  }
}

void set_shared_ids(MPI_Comm comm,
                    const snl_fei::RecordCollection& records,
                    fei::SharedIDs<int>& sharedIDs,
                    int lowest_global_id, int highest_global_id)
{
  sharedIDs.getSharedIDs().clear();
  sharedIDs.getOwningProcs().clear();

  int numProcs = fei::numProcs(comm);
  if (numProcs < 2) return;
  int myProc = fei::localProc(comm);

  fei::LinearDecomposition<int> lindecomp(myProc,numProcs,
                                     lowest_global_id, highest_global_id);

  //Fill a CommMap (procs_to_shared_ids) that maps other procs to ids which we hold.
  //These are ids that appear locally, but which are *not* in our portion
  //of the linear-decomposition.
  fei::CommMap<int>::Type procs_to_shared_ids;
  copy_remotelyowned_ids_into_CommMap(myProc, lindecomp, records, procs_to_shared_ids);

  //Do a global-exchange where we send ids we share to procs that own them,
  //and receive IDs that we own from other procs that share them.
  fei::CommMap<int>::Type procs_to_owned_ids;
  fei::exchangeCommMapData<int>(comm, procs_to_shared_ids, procs_to_owned_ids);

  //transpose procs_to_owned_ids:
  fei::CommMap<int>::Type owned_ids_to_procs;

  fei::CommMap<int>::Type::iterator
    o_iter = procs_to_owned_ids.begin(), o_end = procs_to_owned_ids.end();

  for(; o_iter != o_end; ++o_iter) {
    int proc = o_iter->first;
    std::vector<int>& ids = o_iter->second;
    for(size_t i=0; i<ids.size(); ++i) {
      addItemsToCommMap(ids[i], 1, &proc, owned_ids_to_procs);
      if (records.getRecordWithID(ids[i]) != NULL) {
        addItemsToCommMap(ids[i], 1, &myProc, owned_ids_to_procs);
      }
    }
  }

  fei::CommMap<int>::Type procs_to_owned_ids_and_sharing_procs;

  for(o_iter=procs_to_owned_ids.begin(); o_iter!=procs_to_owned_ids.end(); ++o_iter) {
    int proc = o_iter->first;
    std::vector<int>& ids = o_iter->second;
    for(size_t i=0; i<ids.size(); ++i) {
      std::vector<int>& sharing_procs = owned_ids_to_procs[ids[i]];
      addItemsToCommMap(proc, 1, &ids[i], procs_to_owned_ids_and_sharing_procs, false);
      int num_sharing_procs = sharing_procs.size();
      addItemsToCommMap(proc, 1, &num_sharing_procs, procs_to_owned_ids_and_sharing_procs, false);
      addItemsToCommMap(proc, num_sharing_procs, &sharing_procs[0], procs_to_owned_ids_and_sharing_procs, false);
    }
  }

  fei::CommMap<int>::Type procs_to_shared_ids_and_sharing_procs;
  fei::exchangeCommMapData<int>(comm, procs_to_owned_ids_and_sharing_procs,
                                 procs_to_shared_ids_and_sharing_procs);

  copy_into_shared_ids(procs_to_owned_ids_and_sharing_procs, records, sharedIDs);
  copy_into_shared_ids(procs_to_shared_ids_and_sharing_procs, records, sharedIDs);
}

}//namespace fei

