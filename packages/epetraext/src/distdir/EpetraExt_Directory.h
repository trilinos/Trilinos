//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
#ifndef EPETRAEXT_DIRECTORY_H
#define EPETRAEXT_DIRECTORY_H

#if defined(EpetraExt_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The EpetraExt package is deprecated"
#endif
#endif

// ----------- Includes ----------

#include <map>
#include <vector>

#include <cmath>

#include <Teuchos_RCP.hpp>

#include <EpetraExt_Functors.h>

namespace EpetraExt {

///
/** Distributed Directory Tool
 */
template <typename KT, typename DT, class DH, class AC, class MG>
class Directory
{

public:

  typedef typename std::map< KT, Teuchos::RCP<DT> > DataMap;
  typedef typename DataMap::iterator                DataMapIter;
  typedef typename DataMap::const_iterator          DataMapCIter;

  typedef typename std::multimap< KT, Teuchos::RCP<DT> > DataRecvMap;
  typedef typename DataRecvMap::iterator                 DataRecvMapIter;
  typedef typename DataRecvMap::const_iterator           DataRecvMapCIter;

  typedef typename std::vector<KT>          KeyList;
  typedef typename KeyList::iterator        KeyListIter;
  typedef typename KeyList::const_iterator  KeyListCIter;

  typedef typename std::vector<int>   ProcList;
  typedef typename ProcList::iterator ProcListIter;

  typedef typename std::pair<int,KT> ProcKeyPair;
  typedef typename std::vector<ProcKeyPair> ProcKeyList;
  typedef typename ProcKeyList::iterator ProcKeyListIter;

  typedef typename AC::iterator       ContainerIter;
  typedef typename AC::const_iterator ContainerCIter;

  // Constructors
  Directory( MG migrate,
             DH distHash )
  : migrate_(migrate),
    distHash_(distHash)
  {}

  // Destructor
  ~Directory() {}

private:
  // No public copy construction, assignment, or equality operators
  Directory( const Directory & );

  Directory & operator=( const Directory & );

  bool operator==( const Directory & ) const;
  bool operator!=( const Directory & ) const;

public:

  // Add objects from directory.
  void addEntries( DataMap const & entries );

  // Remove objects from directory.
  void deleteEntries( KeyList & keys );

  // Get the items in the directory.
  void getEntries( KeyList & keys,
                   DataMap & entries );

  AC & container() { return container_; }
  ContainerIter & begin() { return container_.begin(); }
  ContainerIter & end() { return container_.end(); }

protected:

#ifdef EPETRA_MPI
  void pushKeys_( KeyList &, KeyList &, ProcList & );
  void pushData_( DataMap const &, DataRecvMap &, ProcList & );
#endif

  MG migrate_;
  DH distHash_;
  AC container_;

};

///
/*** Hash function for processor assignment
 */
template <typename T>
class Hash
{
  int operator()( const T & in ) { assert(0); return 0; }
};

template <>
class Hash<std::string>
{
  float size_;

 public:

  Hash( int size )
  : size_( static_cast<double>(size) )
  {}

  int operator()( const std::string & in )
  {
    int slen = in.length();
    int sum = 0;
    for( int i = 0; i < slen; ++i )
      sum += static_cast<int>( in[i] );

    return static_cast<int>( fmod( static_cast<float>( sum ), size_ ) );
  }
};

///
/** Sorts a given container: deal with a problem with some STL impl.'s
 */
template < typename T, typename U >
void SortContainer2( T & firstContainer, U & secondContainer )
{
  typedef typename std::multimap< typename T::value_type, typename U::value_type> UTMultiMap;

  UTMultiMap SortMap;

  typename T::iterator iterT = firstContainer.begin();
  typename T::iterator endT = firstContainer.end();
  typename U::iterator iterU = secondContainer.begin();
  typename U::iterator endU = secondContainer.end();

  for( ; (iterT!=endT)||(iterU!=endU) ; ++iterT, ++iterU )
    SortMap.insert( typename UTMultiMap::value_type( *iterT, *iterU ) );

  firstContainer.clear();
  secondContainer.clear();

  typename UTMultiMap::iterator iterUTM = SortMap.begin();
  typename UTMultiMap::iterator endUTM = SortMap.end();

  for( ; iterUTM != endUTM; ++iterUTM )
  {
    firstContainer.push_back( iterUTM->first );
    secondContainer.push_back( iterUTM->second );
  }
}

///
/** Checks if data in a container is sorted
 */
template < typename T >
bool IsSorted( T & container )
{
  if( container.size() < 2 ) return true;

  typename T::iterator iterT = container.begin();
  typename T::iterator endT = container.end();
  typename T::iterator iterTPlus = iterT;
  iterTPlus++;

  for( ; iterTPlus != endT; ++iterT, ++iterTPlus )
    if( !(*iterT<*iterTPlus) ) return false;

  return true;
}

template <typename KT, typename DT, class DH, class AC, class MG>
void
Directory<KT,DT,DH,AC,MG>::
addEntries( DataMap const & entries )
{
#ifdef EPETRA_MPI

  DataRecvMap newEntries;
  ProcList procs;
  pushData_( entries, newEntries, procs );

  DataRecvMapCIter citDM = newEntries.begin();
  DataRecvMapCIter cendDM = newEntries.end();

#else

  DataMapCIter citDM  = entries.begin();
  DataMapCIter cendDM = entries.end();

#endif

  for( ; citDM != cendDM; ++citDM )
      container_.insert( *citDM );
}

template <typename KT, typename DT, class DH, class AC, class MG>
void
Directory<KT,DT,DH,AC,MG>::
deleteEntries( KeyList & keys )
{
#ifdef EPETRA_MPI

  KeyList newKeys;
  ProcList procs;
  pushKeys_( keys, newKeys, procs );

  KeyListCIter citKL = newKeys.begin();
  KeyListCIter cendKL = newKeys.end();

#else

  KeyListCIter citKL  = keys.begin();
  KeyListCIter cendKL = keys.end();

#endif

  for( ; citKL != cendKL; ++citKL )
    container_.erase( *citKL );
}

template <typename KT, typename DT, class DH, class AC, class MG>
void
Directory<KT,DT,DH,AC,MG>::
getEntries( KeyList & keys,
            DataMap & entries )
{
#ifdef EPETRA_MPI

  //Push Keys to owning processors
  KeyList newKeys;
  ProcList procs;
  pushKeys_( keys, newKeys, procs );

  KeyListCIter citKL  = newKeys.begin();
  KeyListCIter cendKL = newKeys.end();

  //Rvs migrate to move data from directory back to requesting procs
  DataMap newEntries;
  for( ; citKL != cendKL; ++citKL )
  {
    if( !container_.count( *citKL ) )
      throw "Data not in directory: " + *citKL + "\n";

    newEntries[*citKL] = (container_.lower_bound( *citKL ))->second;
  }

  migrate_.rvs( procs, newKeys, newEntries, entries );

#else

  KeyListCIter citKL  = keys.begin();
  KeyListCIter cendKL = keys.end();
  for( ; citKL != cendKL; ++citKL )
  {
    if( !container_.count( *citKL ) )
      throw "Data not in directory: " + *citKL + "\n";

    entries[*citKL] = (container_.lower_bound( *citKL ))->second;
  }

#endif
}

#ifdef EPETRA_MPI

template <typename KT, typename DT, class DH, class AC, class MG>
void
Directory<KT,DT,DH,AC,MG>::
pushKeys_( KeyList & sKeys,
           KeyList & rKeys,
           ProcList & procs )
{
  KeyListCIter itKL  = sKeys.begin();
  KeyListCIter endKL = sKeys.end();

  procs.clear();
  for( ; itKL != endKL; ++itKL )
    procs.push_back( distHash_(*itKL) );

  if( !IsSorted( procs ) ) SortContainer2( procs, sKeys );

  migrate_( procs, sKeys, rKeys );
}

template <typename KT, typename DT, class DH, class AC, class MG>
void
Directory<KT,DT,DH,AC,MG>::
pushData_( DataMap const & sData,
           DataRecvMap & rData,
           ProcList & procs )
{
  DataMapCIter itDM  = sData.begin();
  DataMapCIter endDM = sData.end();

  procs.clear();
  for( ; itDM != endDM; ++itDM )
    procs.push_back( distHash_(itDM->first) );

  migrate_( procs, sData, rData );
}

#endif

} //namespace EpetraExt

#endif
