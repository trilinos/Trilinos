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
#ifndef EPETRAEXT_MIGRATE_H
#define EPETRAEXT_MIGRATE_H

#if defined(EpetraExt_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The EpetraExt package is deprecated"
#endif
#endif

// ---------- Includes ----------

#include <string>
#include <vector>
#include <map>

#include <Teuchos_RCP.hpp>

#include <EpetraExt_PackTraits.h>

#ifdef EPETRA_MPI
 #include <Epetra_MpiComm.h>
 #include <Epetra_MpiDistributor.h>
#endif

namespace EpetraExt {

///
/** Data Migration Utility used by EpetraExt::Directory
 */
template <typename KT, typename DT>
class Migrate
{
 public:

  typedef typename std::map< KT, Teuchos::RCP<DT> >         DataMap;
  typedef typename DataMap::iterator        DataMapIter;
  typedef typename DataMap::const_iterator  DataMapCIter;

  typedef typename DataMap::value_type        DataPair;

  typedef typename std::vector<KT>          KeyList;
  typedef typename KeyList::iterator        KeyListIter;
  typedef typename KeyList::const_iterator  KeyListCIter;

  typedef typename std::vector<int>   ProcList;
  typedef typename ProcList::iterator ProcListIter;

  typedef typename std::vector<char> Buffer;

  // Constructor
  Migrate( Epetra_Comm & comm )
  : comm_(comm),
    imports_(0),
    importSize_(0)
  {}

  // Destructor
  ~Migrate()
  { if( importSize_ ) delete [] imports_; }

 private:

  // No public copy construction, assignment, or equality operators.
  Migrate();

  bool operator==( Migrate const & right ) const;
  bool operator!=( Migrate const & right ) const;

 public:

  void operator()( std::vector<int> const & pList,
                   std::vector<KT> const & iKeys,
                   std::vector<KT> & oKeys );

  void operator()( std::vector<int> const & pList,
                   std::map< KT, Teuchos::RCP<DT> > const & iData,
                   std::multimap< KT, Teuchos::RCP<DT> > & oData );

  void rvs( std::vector<int> const & pList,
            std::vector<KT> const & keys,
            std::map< KT, Teuchos::RCP<DT> > & iData,
            std::map< KT, Teuchos::RCP<DT> > & oData );

 protected:

  Epetra_Comm & comm_;

  char * imports_;
  int    importSize_;

  Buffer exports_;

};

template <typename DT>
class Migrate1
{
 public:

  typedef typename Teuchos::RCP<DT> DataPtr;
  typedef typename std::vector<DataPtr>   DataContainer;
  typedef typename DataContainer::iterator        DataContainerIter;
  typedef typename DataContainer::const_iterator  DataContainerCIter;

  typedef typename std::vector<int>   ProcList;
  typedef typename ProcList::iterator ProcListIter;

  typedef typename std::vector<char>  Buffer;

  // Constructor
  Migrate1( Epetra_Comm & comm )
  : comm_(comm),
    imports_(0),
    importSize_(0)
  {}

  // Destructor
  ~Migrate1()
  { if( importSize_ ) delete [] imports_; }

 private:

  // No public copy construction, assignment, or equality operators.
  Migrate1();

  bool operator==( Migrate1 const & right ) const;
  bool operator!=( Migrate1 const & right ) const;

 public:

  void operator()( std::vector<int> const & pList,
                   std::vector< Teuchos::RCP<DT> > const & iData,
                   std::vector< Teuchos::RCP<DT> > & oData );

  void rvs( std::vector<int> const & pList,
            std::vector< Teuchos::RCP<DT> > const & iData,
            std::vector< Teuchos::RCP<DT> > & oData );

 protected:

  Epetra_Comm & comm_;

  char * imports_;
  int    importSize_;

  Buffer exports_;

};

template <typename KT, typename DT>
void
Migrate<KT,DT>::
operator()( std::vector<int> const & pList,
            std::vector<KT> const & iKeys,
            std::vector<KT> & oKeys )
{
#ifdef EPETRA_MPI
  Epetra_MpiDistributor distributor( dynamic_cast<Epetra_MpiComm&>(comm_) );

  int exportCnt = pList.size();

  int max_size = 0;
  KeyListCIter citKL = iKeys.begin();
  KeyListCIter cendKL = iKeys.end();
  for( ; citKL != cendKL; ++citKL )
    max_size = std::max( max_size, static_cast<int>(PackTraits<KT>::size( *citKL )) );

  int importCnt;
  distributor.CreateFromSends( exportCnt, &(pList[0]), true, importCnt );

  int max_all;
  comm_.MaxAll( &max_size, &max_all, 1 );

  exports_.resize( max_all * exportCnt );

  if( importSize_ < (max_all*importCnt) )
  {
    if( importSize_ ) delete [] imports_;
    importSize_ = (max_all*importCnt);
    imports_ = new char[importSize_];
  }

  int pos = 0;
  citKL = iKeys.begin();
  for( int i = 0; citKL != cendKL; ++citKL, ++i )
  {
    pos = max_all * i;
    PackTraits<KT>::pack( *citKL, &(exports_[0]), (max_all*exportCnt ), pos );
  }

  distributor.Do( &(exports_[0]), max_all, importSize_, imports_ );

  oKeys.resize( importCnt );
  for( int i = 0; i < importCnt; ++i )
  {
    pos = max_all * i;
    PackTraits<KT>::unpack( oKeys[i], &(imports_[0]), (max_all*importCnt), pos );
  }
#else
  //Just Copy Data
  oKeys = iKeys;
#endif
}

template <typename KT, typename DT>
void
Migrate<KT,DT>::
operator()( std::vector<int> const & pList,
            std::map< KT, Teuchos::RCP<DT> > const & iData,
            std::multimap< KT, Teuchos::RCP<DT> > & oData )
{
#ifdef EPETRA_MPI
  Epetra_MpiDistributor distributor( dynamic_cast<Epetra_MpiComm&>(comm_) );

  int exportCnt = pList.size();

  int max_size = 0;
  DataMapCIter citDM  = iData.begin();
  DataMapCIter cendDM = iData.end();
  for( ; citDM != cendDM; ++citDM )
    max_size = std::max( max_size, static_cast<int>(PackTraits<KT>::size( citDM->first ))
               + static_cast<int>(PackTraits<DT>::size( *(citDM->second) )) );

  int importCnt;
  distributor.CreateFromSends( exportCnt, &(pList[0]), true, importCnt );

  int max_all;
  comm_.MaxAll( &max_size, &max_all, 1 );

  exports_.resize( max_all * exportCnt );

  if( importSize_ < (max_all*importCnt) )
  {
    if( importSize_ ) delete [] imports_;
    importSize_ = (max_all*importCnt);
    imports_ = new char[importSize_];
  }

  int pos = 0;
  citDM = iData.begin();
  for( int i = 0; citDM != cendDM; ++citDM, ++i )
  {
    pos = max_all * i;
    PackTraits<KT>::pack( citDM->first, &(exports_[0]), (max_all*exportCnt ), pos );
    PackTraits<DT>::pack( *(citDM->second), &(exports_[0]), (max_all*exportCnt ), pos );
  }

  distributor.Do( &(exports_[0]), max_all, importSize_, imports_ );

  oData.clear();
  KT key;
  for( int i = 0; i < importCnt; ++i )
  {
    pos = max_all * i;
    PackTraits<KT>::unpack( key, &(imports_[0]), (max_all*importCnt), pos );
    Teuchos::RCP<DT> data = rcp( new DT );
    PackTraits<DT>::unpack( *data, &(imports_[0]), (max_all*importCnt), pos );
    oData.insert( DataPair( key, data ) );
  }
#else
  //Just Copy Data
  DataMapCIter citDM  = iData.begin();
  DataMapCIter cendDM = iData.end();
  for( ; citDM != cendDM; ++citDM )
    oData.insert( *citDM );
#endif
}

template <typename KT, typename DT>
void
Migrate<KT,DT>::
rvs( std::vector<int> const & pList,
     std::vector<KT> const & keys,
     std::map< KT, Teuchos::RCP<DT> > & iData,
     std::map< KT, Teuchos::RCP<DT> > & oData )
{
#ifdef EPETRA_MPI
  Epetra_MpiDistributor distributor( dynamic_cast<Epetra_MpiComm&>(comm_) );

  int importCnt = pList.size();
  int exportCnt;

  distributor.CreateFromSends( importCnt, &(pList[0]), true, exportCnt );

  if( exportCnt != keys.size() )
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
     "Xyce::Parallel::Migrate::rvs Failed Size Match!\n" );

  int max_size = 0;
  KeyListCIter citKL  = keys.begin();
  KeyListCIter cendKL = keys.end();
  for( ; citKL != cendKL; ++citKL )
    max_size = std::max( max_size, static_cast<int>(PackTraits<KT>::size( *citKL ) )
                                 + static_cast<int>(PackTraits<DT>::size( *(iData[*citKL]) ) ) );

  int max_all;
  comm_.MaxAll( &max_size, &max_all, 1 );

  exports_.resize( max_all * exportCnt );

  if( importSize_ < (max_all*importCnt) )
  {
    if( importSize_ ) delete [] imports_;
    importSize_ = (max_all*importCnt);
    imports_ = new char[importSize_];
  }

  int pos = 0;
  int i = 0;
  citKL  = keys.begin();
  for( ; citKL != cendKL; ++citKL, ++i )
  {
    pos = max_all * i;
    PackTraits<KT>::pack( *citKL, &(exports_[0]), (max_all*exportCnt ), pos );
    PackTraits<DT>::pack( *(iData[*citKL]), &(exports_[0]), (max_all*exportCnt ), pos );
  }

  distributor.DoReverse( &(exports_[0]), max_all, importSize_, imports_ );

  oData.clear();
  KT key;
  for( int i = 0; i < importCnt; ++i )
  {
    pos = max_all * i;
    PackTraits<KT>::unpack( key, &(imports_[0]), (max_all*importCnt), pos );
    Teuchos::RCP<DT> data = rcp( new DT );
    PackTraits<DT>::unpack( *data, &(imports_[0]), (max_all*importCnt), pos );
    oData[key] = data;
  }
#else
  oData = iData;
#endif
}

template <typename DT>
void
Migrate1<DT>::
operator()( std::vector<int> const & pList,
            std::vector< Teuchos::RCP<DT> > const & iData,
            std::vector< Teuchos::RCP<DT> > & oData )
{
#ifdef EPETRA_MPI
  Epetra_MpiDistributor distributor( dynamic_cast<Epetra_MpiComm&>(comm_) );

  int exportCnt = pList.size();

  int max_size = 0;
  DataContainerCIter citDC = iData.begin();
  DataContainerCIter cendDC = iData.end();
  for( ; citDC != cendDC; ++citDC )
    max_size = std::max( max_size, static_cast<int>(PackTraits<DT>::size( **citDC )) );

  int importCnt;
  distributor.CreateFromSends( exportCnt, &(pList[0]), true, importCnt );

  int max_all;
  comm_.MaxAll( &max_size, &max_all, 1 );

  exports_.resize( max_all * exportCnt );

  if( importSize_ < (max_all*importCnt) )
  {
    if( importSize_ ) delete [] imports_;
    importSize_ = (max_all*importCnt);
    imports_ = new char[importSize_];
  }

  int pos = 0;
  citDC = iData.begin();
  for( int i = 0; citDC != cendDC; ++citDC, ++i )
  {
    pos = max_all * i;
    PackTraits<DT>::pack( **citKL, &(exports_[0]), (max_all*exportCnt ), pos );
  }

  distributor.Do( &(exports_[0]), max_all, importSize_, imports_ );

  oData.clear();
  for( int i = 0; i < importCnt; ++i )
  {
    pos = max_all * i;
    Teuchos::RCP<DT> data = rcp( new DT );
    PackTraits<DT>::unpack( *data, &(imports_[0]), (max_all*importCnt), pos );
    oData.push_back( data );
  }
#else
  //Just Copy Data
  oData = iData;
#endif
}

template <typename DT>
void
Migrate1<DT>::
rvs( std::vector<int> const & pList,
     std::vector< Teuchos::RCP<DT> > const & iData,
     std::vector< Teuchos::RCP<DT> > & oData )
{
#ifdef EPETRA_MPI
  Epetra_MpiDistributor distributor( dynamic_cast<Epetra_MpiComm&>(comm_) );

  int importCnt = pList.size();
  int exportCnt;

  distributor.CreateFromSends( importCnt, &(pList[0]), true, exportCnt );

  if( exportCnt != keys.size() )
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
     "Xyce::Parallel::Migrate::rvs Failed Size Match!\n" );

  int max_size = 0;
  DataContainerCIter citDC  = iData.begin();
  DataContainerCIter cendDC = iData.end();
  for( ; citDC != cendDC; ++citDC )
    max_size = std::max( max_size, PackTraits<DT>::size( **citDC ) );

  int max_all;
  comm_.MaxAll( &max_size, &max_all, 1 );

  exports_.resize( max_all * exportCnt );

  if( importSize_ < (max_all*importCnt) )
  {
    if( importSize_ ) delete [] imports_;
    importSize_ = (max_all*importCnt);
    imports_ = new char[importSize_];
  }

  int i = 0;
  int pos = 0;
  citDC  = iData.begin();
  for( ; citDC != cendDC; ++citDC, ++i )
  {
    pos = max_all * i;
    PackTraits<DT>::pack( **citDC, &(exports_[0]), (max_all*exportCnt ), pos );
  }

  distributor.DoReverse( &(exports_[0]), max_all, importSize_, imports_ );

  oData.clear();
  for( int i = 0; i < importCnt; ++i )
  {
    pos = max_all * i;
    Teuchos::RCP<DT> data = rcp( new DT );
    PackTraits<DT>::unpack( *data, &(imports_[0]), (max_all*importCnt), pos );
    oData.push_back( data );
  }
#else
  oData = iData;
#endif
}

} //namespace EpetraExt

#endif
