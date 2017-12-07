// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
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

#ifndef stk_util_util_OctTreeOps_hpp
#define stk_util_util_OctTreeOps_hpp

#include <utility>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <iostream>
#include <stdexcept>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <stk_util/environment/ReportHandler.hpp>

namespace stk {
namespace search {

template <class DomainBoundingBox, class RangeBoundingBox>
void communicate(
  stk::ParallelMachine arg_comm ,
  const std::set< std::pair< typename DomainBoundingBox::second_type,  typename RangeBoundingBox::second_type > > & send_relation ,
        std::set< std::pair< typename DomainBoundingBox::second_type,  typename RangeBoundingBox::second_type > > & recv_relation ,
        bool communicateRangeBoxInfo = true )
{
  typedef typename DomainBoundingBox::second_type DomainKey;
  typedef typename RangeBoundingBox::second_type RangeKey;
  typedef std::pair<DomainKey, RangeKey> ValueType ;

  CommSparse commSparse( arg_comm );

  const int p_rank = commSparse.parallel_rank();
  const int p_size = commSparse.parallel_size();

  typename std::set< ValueType >::const_iterator i ;

  for ( i = send_relation.begin() ; i != send_relation.end() ; ++i ) {
    const ValueType & val = *i ;
    if ( static_cast<int>(val.first.proc()) == p_rank || ( communicateRangeBoxInfo && static_cast<int>(val.second.proc()) == p_rank) )
    {
      recv_relation.insert( val );
    }
    if ( static_cast<int>(val.first.proc()) != p_rank ) {
      CommBuffer & buf = commSparse.send_buffer( val.first.proc() );
      buf.skip<ValueType>( 1 );
    }
    if ( communicateRangeBoxInfo )
    {
        if ( static_cast<int>(val.second.proc()) != p_rank && val.second.proc() != val.first.proc() ) {
          CommBuffer & buf = commSparse.send_buffer( val.second.proc() );
          buf.skip<ValueType>( 1 );
        }
    }
  }

  commSparse.allocate_buffers();

  for ( i = send_relation.begin() ; i != send_relation.end() ; ++i ) {
    const ValueType & val = *i ;
    if ( static_cast<int>(val.first.proc()) != p_rank ) {
      CommBuffer & buf = commSparse.send_buffer( val.first.proc() );
      buf.pack<ValueType>( val );
    }
    if ( communicateRangeBoxInfo )
    {
        if ( static_cast<int>(val.second.proc()) != p_rank && val.second.proc() != val.first.proc() ) {
          CommBuffer & buf = commSparse.send_buffer( val.second.proc() );
          buf.pack<ValueType>( val );
        }
    }
  }

  commSparse.communicate();

  for ( int p = 0 ; p < p_size ; ++p ) {
    CommBuffer & buf = commSparse.recv_buffer( p );
    while ( buf.remaining() ) {
      ValueType val ;
      buf.unpack<ValueType>( val );
      recv_relation.insert( val );
    }
  }
}


template <typename DomainKey, typename RangeKey>
void communicateVector(
  stk::ParallelMachine arg_comm ,
  const std::vector< std::pair< DomainKey, RangeKey> > & send_relation ,
        std::vector< std::pair< DomainKey, RangeKey> > & recv_relation ,
        bool communicateRangeBoxInfo = true )
{
  typedef std::pair<DomainKey, RangeKey> ValueType ;

  CommSparse commSparse( arg_comm );

  const int p_rank = commSparse.parallel_rank();
  const int p_size = commSparse.parallel_size();

  typename std::vector< ValueType >::const_iterator i ;

  for ( i = send_relation.begin() ; i != send_relation.end() ; ++i ) {
    const ValueType & val = *i ;
    if ( static_cast<int>(val.first.proc()) == p_rank || ( communicateRangeBoxInfo && static_cast<int>(val.second.proc()) == p_rank) )
    {
      recv_relation.push_back( val );
    }
    if ( static_cast<int>(val.first.proc()) != p_rank ) {
      CommBuffer & buf = commSparse.send_buffer( val.first.proc() );
      buf.skip<ValueType>( 1 );
    }
    if ( communicateRangeBoxInfo )
    {
        if ( static_cast<int>(val.second.proc()) != p_rank && val.second.proc() != val.first.proc() ) {
          CommBuffer & buf = commSparse.send_buffer( val.second.proc() );
          buf.skip<ValueType>( 1 );
        }
    }
  }

  commSparse.allocate_buffers();

  for ( i = send_relation.begin() ; i != send_relation.end() ; ++i ) {
    const ValueType & val = *i ;
    if ( static_cast<int>(val.first.proc()) != p_rank ) {
      CommBuffer & buf = commSparse.send_buffer( val.first.proc() );
      buf.pack<ValueType>( val );
    }
    if ( communicateRangeBoxInfo )
    {
        if ( static_cast<int>(val.second.proc()) != p_rank && val.second.proc() != val.first.proc() ) {
          CommBuffer & buf = commSparse.send_buffer( val.second.proc() );
          buf.pack<ValueType>( val );
        }
    }
  }

  commSparse.communicate();

  for ( int p = 0 ; p < p_size ; ++p ) {
    CommBuffer & buf = commSparse.recv_buffer( p );
    while ( buf.remaining() ) {
      ValueType val ;
      buf.unpack<ValueType>( val );
      recv_relation.push_back( val );
    }
  }
}

} // namespace search
} // namespace stk

#endif

