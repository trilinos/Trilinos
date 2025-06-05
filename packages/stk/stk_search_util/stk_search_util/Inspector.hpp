/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_INSPECTOR_HPP_
#define STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_INSPECTOR_HPP_

#include <stk_mesh/base/Types.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_search_util/MeshUtility.hpp>

#include <vector>
#include <string>
#include <memory>
#include <fstream>

namespace stk {
namespace search {

template <typename SENDMESH, typename RECVMESH>
struct InspectorInfo {
  using EntityKeyA = typename SENDMESH::EntityKey;
  using EntityKeyB = typename RECVMESH::EntityKey;

  EntityKeyA domainEntityKey;
  EntityKeyB rangeEntityKey;

  int domainProc;
  int rangeProc;

  std::vector<std::string> domainParts;
};

template <typename SENDMESH, typename RECVMESH>
struct InspectorOutputFormatterBase
{
public:
  virtual void write_info(std::ostream& os, const std::vector< InspectorInfo<SENDMESH, RECVMESH> >& infoVec) const
  {
    if(!m_isFirstWrite) {
      os << std::endl;
    }
    write_delimiter(os);
    for( const InspectorInfo<SENDMESH, RECVMESH>& info : infoVec) {
      internal_write_info(os, info);
    }
    m_isFirstWrite = false;
  }

  virtual ~InspectorOutputFormatterBase() {}

protected:
  virtual void internal_write_info(std::ostream& os, const InspectorInfo<SENDMESH, RECVMESH>& info) const = 0;
  virtual void write_delimiter(std::ostream& os) const { os << get_time_stamp() << std::endl; }
  virtual std::string get_time_stamp() const
  {
    static const std::string defaultStamp("------------------------------------");
    return defaultStamp;
  }

  mutable bool m_isFirstWrite = true;
};

template <typename SENDMESH, typename RECVMESH>
struct InspectorOutputFormatter : public InspectorOutputFormatterBase<SENDMESH, RECVMESH>
{
public:
  InspectorOutputFormatter(std::shared_ptr<SENDMESH>& mesha, std::shared_ptr<RECVMESH>& meshb)
   : m_sendMesh(mesha)
   , m_recvMesh(meshb)
   {}

  const std::shared_ptr<SENDMESH> send_mesh() const { return m_sendMesh; }
  const std::shared_ptr<RECVMESH> recv_mesh() const { return m_recvMesh; }

protected:
  void write_delimiter(std::ostream& os) const override
  {
    os << m_recvMesh->get_inspector_delimiter() << " -> "
       << m_sendMesh->get_inspector_delimiter() << std::endl;
  }

  void internal_write_info(std::ostream& os, const InspectorInfo<SENDMESH, RECVMESH>& info) const override
  {
    os << "[" << info.rangeEntityKey  << ", PROC(" << info.rangeProc  << ")]"
       << " -> "
       << "[" << info.domainEntityKey << ", PROC(" << info.domainProc << ")]"
       << " : ";
    for(const std::string& partName : info.domainParts) {
      os << partName << " ";
    }
    os << std::endl;
  }

private:
  std::shared_ptr<SENDMESH> m_sendMesh;
  std::shared_ptr<RECVMESH> m_recvMesh;

  InspectorOutputFormatter() {}
};

template <typename SENDMESH, typename RECVMESH>
class Inspector
{
public:
  using EntityKeyA = typename SENDMESH::EntityKey;
  using EntityKeyB = typename RECVMESH::EntityKey;
  using EntityProcA = typename SENDMESH::EntityProc;
  using EntityProcB = typename RECVMESH::EntityProc;

  using EntityProcRelation = std::pair<EntityProcB, EntityProcA>;
  using EntityProcRelationVec = std::vector<EntityProcRelation>;
  using EntityKeyVecA = std::vector<EntityKeyA>;
  using EntityKeyVecB = std::vector<EntityKeyB>;

  Inspector(std::shared_ptr<SENDMESH>& mesha, std::shared_ptr<RECVMESH>& meshb,
        const std::string& fileName, const EntityKeyVecB& rangeEntities)
  : m_sendMesh(mesha)
  , m_recvMesh(meshb)
  , m_fileName(fileName)
  , m_rangeEntities(rangeEntities)
  , m_outputFormatter(std::make_shared<InspectorOutputFormatter<SENDMESH,RECVMESH>> (m_sendMesh, m_recvMesh))
  {
    stk::ParallelMachine comm = m_sendMesh->comm();
    if(stk::parallel_machine_rank(comm) == m_rootProc) {
      m_outputStream.open(m_fileName, std::ios::trunc);
    }
  }

  Inspector()
  : m_sendMesh(nullptr)
  , m_recvMesh(nullptr)
  , m_outputFormatter(std::make_shared<InspectorOutputFormatter<SENDMESH,RECVMESH>> (m_sendMesh, m_recvMesh))
  { }

  Inspector(const Inspector& original)
  : m_sendMesh(original.m_sendMesh)
  , m_recvMesh(original.m_recvMesh)
  , m_fileName(original.m_fileName)
  , m_outputFormatter(original.m_outputFormatter)
  {
    stk::ParallelMachine comm = m_sendMesh->comm();
    if(stk::parallel_machine_rank(comm) == m_rootProc) {
      m_outputStream.open(m_fileName, std::ios::app);
    }
  }

  void register_output_formatter(std::shared_ptr< InspectorOutputFormatterBase<SENDMESH, RECVMESH> > outputFormatter)
  {
    m_outputFormatter = outputFormatter;
  }

  ~Inspector() = default;

  void execute(const EntityProcRelationVec& rangeToDomain);

  const std::string& get_file_name() const {return m_fileName;}

private:
  std::shared_ptr<SENDMESH> m_sendMesh;
  std::shared_ptr<RECVMESH> m_recvMesh;

  std::string m_fileName;
  EntityKeyVecB m_rangeEntities;

  std::ofstream m_outputStream;

  int m_rootProc = 0;

  std::shared_ptr< InspectorOutputFormatterBase<SENDMESH, RECVMESH> > m_outputFormatter;

  InspectorInfo<SENDMESH, RECVMESH> extract_inspector_info(const EntityProcRelation& relation);
  std::vector< InspectorInfo<SENDMESH, RECVMESH> > get_info_vector(const EntityProcRelationVec& rangeToDomain);
  void communicate_inspector_info(const std::vector< InspectorInfo<SENDMESH, RECVMESH> >& infoVec,
                                  std::vector< InspectorInfo<SENDMESH, RECVMESH> >& accumulatedInfoVec);

  bool is_output_proc() const
  {
    stk::ParallelMachine comm = m_sendMesh->comm();
    int myProc = stk::parallel_machine_rank(comm);

    return myProc == m_rootProc;
  }
};

template <typename SENDMESH, typename RECVMESH>
InspectorInfo<SENDMESH, RECVMESH> Inspector<SENDMESH, RECVMESH>::extract_inspector_info(const EntityProcRelation& relation)
{
  InspectorInfo<SENDMESH, RECVMESH> info;

  info.domainEntityKey = relation.second.id();
  info.rangeEntityKey = relation.first.id();
  info.domainProc = relation.second.proc();
  info.rangeProc = relation.first.proc();
  info.domainParts = m_sendMesh->get_part_membership(info.domainEntityKey);

  return info;
}

template <typename SENDMESH, typename RECVMESH>
std::vector< InspectorInfo<SENDMESH, RECVMESH> >
Inspector<SENDMESH, RECVMESH>::get_info_vector(const EntityProcRelationVec& rangeToDomain)
{
  std::vector< InspectorInfo<SENDMESH, RECVMESH> > infoVec;

  auto comparator = [] (const EntityProcRelation& relation, const EntityProcB& entry) {
    return relation.first < entry;
  };

  for(const auto& object : m_rangeEntities) {
    EntityProcB minEntry, maxEntry;

    minEntry.set_proc(0);
    minEntry.set_id(object);

    maxEntry.set_proc(std::numeric_limits<int>::max());
    maxEntry.set_id(object);

    auto minIter = std::lower_bound(rangeToDomain.begin(), rangeToDomain.end(), minEntry, comparator);
    auto maxIter = std::lower_bound(rangeToDomain.begin(), rangeToDomain.end(), maxEntry, comparator);

    auto extractor = [&infoVec, &object, this] (const EntityProcRelation& entry) {
      InspectorInfo<SENDMESH, RECVMESH> info = this->extract_inspector_info(entry);

      if(object == info.rangeEntityKey) {
        infoVec.push_back(info);
      }
    };

    std::for_each(minIter, maxIter, extractor);
  }

  return infoVec;
}

template <typename SENDMESH, typename RECVMESH>
void Inspector<SENDMESH, RECVMESH>::communicate_inspector_info(const std::vector< InspectorInfo<SENDMESH, RECVMESH> >& infoVec,
                                                               std::vector< InspectorInfo<SENDMESH, RECVMESH> >& accumulatedInfoVec)
{
  stk::ParallelMachine comm = m_sendMesh->comm();
  stk::CommSparse commSparse(comm);
  int myProc = stk::parallel_machine_rank(comm);

  if(myProc == m_rootProc) {
    for(const InspectorInfo<SENDMESH, RECVMESH>& info : infoVec) {
      if((info.rangeProc == myProc) && (myProc == m_rootProc)) {
        accumulatedInfoVec.push_back(info);
      }
    }
  }

  int rootProc = m_rootProc;
  stk::pack_and_communicate(commSparse, [&commSparse, &infoVec, &myProc, &rootProc]()
                           {
                              stk::CommBuffer& buf = commSparse.send_buffer(rootProc);

                              for(const InspectorInfo<SENDMESH, RECVMESH>& info : infoVec) {
                                if((info.rangeProc == myProc) && (myProc != rootProc)) {
                                  buf.pack<EntityKeyA>(info.domainEntityKey);
                                  buf.pack<EntityKeyB>(info.rangeEntityKey);
                                  buf.pack<int>(info.domainProc);
                                  buf.pack<int>(info.rangeProc);
                                  buf.pack<size_t>(info.domainParts.size());

                                  for(const std::string& domainPart : info.domainParts) {
                                    buf.pack(domainPart);
                                  }
                                }
                              }
                           });

  stk::unpack_communications(commSparse, [&commSparse, &accumulatedInfoVec](int procId)
                            {
                               InspectorInfo<SENDMESH, RECVMESH> info;
                               stk::CommBuffer& buf = commSparse.recv_buffer(procId);

                               buf.unpack<EntityKeyA>(info.domainEntityKey);
                               buf.unpack<EntityKeyB>(info.rangeEntityKey);
                               buf.unpack<int>(info.domainProc);
                               buf.unpack<int>(info.rangeProc);

                               size_t numParts;
                               buf.unpack<size_t>(numParts);
                               info.domainParts.resize(numParts);

                               for(std::string& domainPart : info.domainParts) {
                                 buf.unpack(domainPart);
                               }

                               accumulatedInfoVec.push_back(info);
                            });

}

template <typename SENDMESH, typename RECVMESH>
void Inspector<SENDMESH, RECVMESH>::execute(const EntityProcRelationVec& rangeToDomain)
{
  std::vector< InspectorInfo<SENDMESH, RECVMESH> > infoVec = get_info_vector(rangeToDomain);
  std::vector< InspectorInfo<SENDMESH, RECVMESH> > accumulatedInfoVec;
  communicate_inspector_info(infoVec, accumulatedInfoVec);

  if(is_output_proc()) {
    m_outputFormatter->write_info(m_outputStream, accumulatedInfoVec);
  }
}

}}



#endif /* STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_INSPECTOR_HPP_ */
