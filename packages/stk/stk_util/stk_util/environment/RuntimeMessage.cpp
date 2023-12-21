// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#include "stk_util/environment/RuntimeMessage.hpp"
#include "stk_util/environment/Env.hpp"
#include "stk_util/parallel/Parallel.hpp"   // for parallel_machine_rank, MPI_Gather, MPI_Gatherv
#include "stk_util/stk_config.h"            // for STK_HAS_MPI
#include "stk_util/util/Bootstrap.hpp"      // for Bootstrap
#include "stk_util/util/Marshal.hpp"        // for operator<<, operator>>, Marshal
#include "stk_util/util/ReportHandler.hpp"  // for report
#include <algorithm>                        // for stable_sort
#include <functional>                       // for hash
#include <sstream>                          // for operator<<, basic_ostream, basic_ostream::ope...
#include <stdexcept>                        // for runtime_error
#include <string>                           // for string, operator<<, char_traits, operator==
#include <system_error>
#include <unordered_map>                    // for unordered_map<>::iterator, _Node_iterator
#include <utility>                          // for pair
#include <vector>                           // for vector, vector<>::const_iterator


namespace stk {
  typedef std::pair<MessageId, std::string> MessageKey;
}

namespace std {

template<>
struct hash<stk::MessageKey>
{
    size_t operator()(const stk::MessageKey& msgkey) const
    {
        return hash_messagekey(msgkey);
    }
    size_t operator()(const stk::MessageKey& msgkey1, const stk::MessageKey& msgkey2) const
    {
        return hash_messagekey(msgkey1)^hash_messagekey(msgkey2);
    }

private:
    size_t hash_messagekey(const stk::MessageKey& msgkey) const
    {
        return hash<stk::MessageId>()(msgkey.first + hash<std::string>()(msgkey.second));
    }
};

}

namespace stk {

MessageCode
MessageCode::s_defaultMessageCode(100000000);

namespace {

void bootstrap()
{
  register_message_type(MSG_WARNING, 10000000, "warning");
  register_message_type(MSG_DOOMED, 10000000, "error");
  register_message_type(MSG_EXCEPTION, 1000000, "exception");
  register_message_type(MSG_INFORMATION, 1000000, "informational");
}

stk::Bootstrap x(bootstrap);

typedef std::unordered_map<MessageKey, Throttle> MessageIdMap;

MessageIdMap s_messageIdMap;

MessageIdMap s_deferredMessageIdMap;

struct DeferredMessage 
{
  DeferredMessage()
  {}
  
  DeferredMessage(
    size_t              type,
    MessageId           message_id,
    size_t              throttle_cutoff,
    int                 throttle_group,
    const std::string & header,
    const std::string & aggregate)
    : m_type(type),
      m_messageId(message_id),
      m_rank(0),
      m_throttleCutoff(throttle_cutoff),
      m_throttleGroup(throttle_group),
      m_header(header),
      m_aggregate(aggregate)
  {}
  
  size_t                m_type;
  MessageId             m_messageId;
  int                   m_rank;
  size_t                m_throttleCutoff;
  int                   m_throttleGroup;
  std::string           m_header;
  std::string           m_aggregate;
};

typedef std::vector<DeferredMessage> DeferredMessageVector;

struct DeferredMessageLess
{
  bool operator()(const DeferredMessage &key_1, const DeferredMessage &key_2) const {
    return (key_1.m_type < key_2.m_type)
      || (!(key_2.m_type < key_1.m_type) && key_1.m_messageId < key_2.m_messageId)
      || (!(key_2.m_type < key_1.m_type) && !(key_2.m_messageId < key_1.m_messageId) && key_1.m_header < key_2.m_header);
  }
};
      
DeferredMessageVector s_deferredMessageVector;

struct MessageTypeInfo
{
  MessageTypeInfo()
    : m_count(0),
      m_maxCount(10000000),
      m_name("unknown")
  {}
  
  unsigned              m_count;
  unsigned              m_maxCount;
  std::string           m_name;
};

typedef std::unordered_map<unsigned, MessageTypeInfo> MessageTypeInfoMap;

MessageTypeInfoMap s_messageTypeInfo;

MessageTypeInfo &
get_message_type_info(
  unsigned              type)
{
  MessageTypeInfoMap::iterator it = s_messageTypeInfo.find(type & MSG_TYPE_MASK);
  if (it != s_messageTypeInfo.end()) {
    return (*it).second;
  }
  return s_messageTypeInfo[type & MSG_TYPE_MASK];
}

#ifdef STK_HAS_MPI
Marshal &operator<<(Marshal &mout, const DeferredMessage &s)  {
  mout << s.m_type << s.m_messageId << s.m_rank << s.m_throttleGroup << s.m_throttleCutoff << s.m_header << s.m_aggregate;
  return mout;
}

Marshal &operator>>(Marshal &min, DeferredMessage &s)  {
  min >> s.m_type >> s.m_messageId >> s.m_rank >> s.m_throttleGroup >> s.m_throttleCutoff >> s.m_header >> s.m_aggregate;
  return min;
}
#endif // STK_HAS_MPI

} // namespace <empty>

bool should_increment_message_count(unsigned messageType)
{
  return !(messageType & MSG_SYMMETRIC) || (stk::parallel_machine_rank(sierra::Env::parallel_comm()) == 0);
}

CutoffStatus
count_message(unsigned messageType, const MessageCode &   message_code)
{

    MessageId message_id     = message_code.m_id;
    const Throttle &throttle = message_code.m_throttle;

  std::pair<MessageIdMap::iterator, bool> res = s_messageIdMap.insert(MessageIdMap::value_type(MessageIdMap::key_type(message_id, std::string("")), throttle));
  if (should_increment_message_count(messageType)) {
    ++(*res.first).second.m_count;
  }

  size_t count = (*res.first).second.m_count;

  if (count < (*res.first).second.m_cutoff) {
      return CutoffStatus::MSG_DISPLAY;
  } else if (count == (*res.first).second.m_cutoff) {
      return CutoffStatus::MSG_CUTOFF;
  } else {
      return CutoffStatus::MSG_CUTOFF_EXCEEDED;
  }
}

void
register_message_type(
  unsigned              message_type,
  unsigned              max_count,
  const char *          name)
{
  MessageTypeInfo &message_info = get_message_type_info(message_type);

  message_info.m_maxCount = max_count;
  message_info.m_name = name;
}

  
unsigned
get_message_count(
  unsigned              message_type)
{
  return get_message_type_info(message_type).m_count;
}

unsigned get_message_printed_count(unsigned messageType)
{
  STK_ThrowRequireMsg(messageType==MSG_WARNING,"Only count printed warning messages.");
  unsigned totalPrinted = 0;
  for(auto& iter : s_messageIdMap) {
    totalPrinted += std::min(iter.second.m_count, iter.second.m_cutoff);
  }
  return totalPrinted;
}

unsigned get_message_printed_count(const MessageCode& messageCode)
{
  MessageIdMap::iterator iter = s_messageIdMap.find(MessageIdMap::key_type(messageCode.m_id, std::string("")));
  if (iter != s_messageIdMap.end()) {
    return std::min(iter->second.m_count, iter->second.m_cutoff);
  }
  return 0;
}


void increment_message_count(unsigned message_type)
{
  if (should_increment_message_count(message_type)) {
    ++get_message_type_info(message_type).m_count;
  }
}

void
reset_message_count(
  unsigned              message_type)
{
  get_message_type_info(message_type).m_count = 0;
  for(auto& iter : s_messageIdMap) {
    iter.second.m_count = 0;
  }
}


const std::string &
get_message_name(
  unsigned              message_type)
{
  return get_message_type_info(message_type).m_name;
}


void
set_max_message_count(
  unsigned              message_type,
  unsigned              max_count)
{
  get_message_type_info(message_type).m_maxCount = max_count;
}


unsigned
get_max_message_count(
  unsigned              message_type)
{
  return get_message_type_info(message_type).m_maxCount;
}


void
report_message(
  const char *          message,
  unsigned              message_type, 
  const MessageCode &   message_code)
{
  if (message_type & MSG_DEFERRED) {
    report(message, message_type);
  } else { 
    increment_message_count(message_type);
    const unsigned count = get_message_count(message_type);
    const unsigned max_count = get_max_message_count(message_type); 
  
    if (count == max_count) {
      report(message, message_type);

      std::ostringstream s;
      s << "Maximum " << get_message_name(message_type) << " count has been exceeded and will no longer be displayed";
      report(s.str().c_str(), MSG_WARNING | MSG_SYMMETRIC);
    }

    else if (count < max_count) {
      CutoffStatus cutoff = count_message(message_type, message_code);
    
      if (cutoff == CutoffStatus::MSG_CUTOFF) {
        report(message, message_type);

        std::ostringstream s;
        s << "Maximum count for this " << get_message_name(message_type) << " (previous message)"
          << " has been exceeded and will no longer be displayed";
        report(s.str().c_str(), MSG_WARNING | MSG_SYMMETRIC);
      }
    
      else if (cutoff == CutoffStatus::MSG_DISPLAY)
      {
        report(message, message_type);
      }
    }
  }
}

void
add_deferred_message(
  int                   message_type,
  MessageId             message_id,
  size_t                throttle_cutoff,
  int                   throttle_group,
  const char *          header,
  const char *          aggegrate)
{
  std::ostringstream s;
  s << header << " " << aggegrate;
  
  report(s.str().c_str(), message_type | MSG_DEFERRED);

  std::pair<MessageIdMap::iterator, bool> res = s_deferredMessageIdMap.insert(MessageIdMap::value_type(MessageIdMap::key_type(message_id, header), Throttle(throttle_cutoff, throttle_group)));
  size_t count = ++(*res.first).second.m_count;

  if (count <= throttle_cutoff) {
    s_deferredMessageVector.push_back(DeferredMessage(message_type, message_id, throttle_cutoff, throttle_group, header, aggegrate));
  }
}


/// @todo REFACTOR Should the message counts be broadcast to the other processors?

void
report_deferred_messages(
  MPI_Comm       comm)
{
#ifdef STK_HAS_MPI
  const int p_root = 0 ;
  int p_size = stk::parallel_machine_size(comm);
  int p_rank = stk::parallel_machine_rank(comm);

  for (auto it = s_deferredMessageVector.begin(); it != s_deferredMessageVector.end(); ++it)
  {
    (*it).m_rank = p_rank;
  }
  
  Marshal mout;
  mout << s_deferredMessageVector;

  DeferredMessageVector deferred_message_vector;

  // Gather the send counts on root processor
  std::string send_string(mout.stream.str());
  int send_count = send_string.size();
  std::vector<int> recv_count(p_size, 0);
  int * const recv_count_ptr = recv_count.data() ;

  int result = MPI_Gather(&send_count, 1, MPI_INT, recv_count_ptr, 1, MPI_INT, p_root, comm);

  if (MPI_SUCCESS != result) {
    std::ostringstream message ;
    message << "stk::report_deferred_messages FAILED: MPI_Gather = " << result ;
    throw std::runtime_error(message.str());
  }

  // Receive counts are only non-zero on the root processor:
  std::vector<int> recv_displ(p_size + 1, 0);

  for (int i = 0 ; i < p_size ; ++i) {
    recv_displ[i + 1] = recv_displ[i] + recv_count[i] ;
  }

  const int recv_size = recv_displ[p_size] ;
 
  std::vector<char> buffer(recv_size);

  {
    const char * const send_ptr = send_string.data();
    char * const recv_ptr = recv_size ? buffer.data() : nullptr ;
    int * const recv_displ_ptr = recv_displ.data() ;

    result = MPI_Gatherv(const_cast<char*>(send_ptr), send_count, MPI_CHAR,
                         recv_ptr, recv_count_ptr, recv_displ_ptr, MPI_CHAR,
                         p_root, comm);
    if (MPI_SUCCESS != result) {
      std::ostringstream message ;
      message << "stk::report_deferred_messages FAILED: MPI_Gatherv = " << result ;
      throw std::runtime_error(message.str());
    }

    if (p_rank == p_root) {
      for (int i = 0; i < p_size; ++i) {
        Marshal min(std::string(recv_ptr + recv_displ[i], recv_ptr + recv_displ[i + 1]));
        min >> deferred_message_vector;
      }

      std::stable_sort(deferred_message_vector.begin(), deferred_message_vector.end(), DeferredMessageLess());      

      DeferredMessageVector::const_iterator current_message_it = deferred_message_vector.begin();
      while (current_message_it != deferred_message_vector.end()) {
        const DeferredMessage &current_message = (*current_message_it);
        
        DeferredMessageVector::const_iterator end = current_message_it + 1;
        while (end != deferred_message_vector.end()
               && current_message.m_messageId == (*end).m_messageId
               && current_message.m_header == (*end).m_header)
          ++end;
        
        std::ostringstream s;
        
        s << current_message.m_header << current_message.m_aggregate;

        for (DeferredMessageVector::const_iterator it1 = current_message_it + 1; it1 != end; ++it1) {  
          bool print = true;
          for (DeferredMessageVector::const_iterator it2 = current_message_it; it2 != it1; ++it2)
            if ((*it1).m_aggregate == (*it2).m_aggregate) {
              print = false;
              break;
            }
          if (print) {
            if (!(*it1).m_aggregate.find('\n'))
              s << ", ";
            s << (*it1).m_aggregate;
          }
        }

        report_message(s.str().c_str(), current_message.m_type | stk::MSG_SYMMETRIC, MessageCode(current_message.m_messageId, current_message.m_throttleCutoff, current_message.m_throttleGroup));
        
        current_message_it = end;
      }
    }
  }

  s_deferredMessageIdMap.clear();
  s_deferredMessageVector.clear();
#endif
}

std::ostream &
operator<<(
  std::ostream &        os,
  const MessageType &   message_type) 
{
  os << get_message_type_info(message_type).m_name;
  return os;
}

} // namespace stk

