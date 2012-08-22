/*
//@HEADER
//************************************************************************
//
//              Isorropia: Partitioning and Load Balancing Package
//                Copyright (2006) Sandia Corporation
//
//Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
//license for use of this work by or on behalf of the U.S. Government.
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
//************************************************************************
//@HEADER
*/

#ifndef EPETRA_LEVELSOLVER_HPP_
#define EPETRA_LEVELSOLVER_HPP_

#include <Epetra_Operator.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Import.h>
#include <Epetra_SerialComm.h>
#include <Epetra_LocalMap.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_DefaultSerialComm.hpp>

#include <Isorropia_LevelScheduler.hpp> 
#include <Isorropia_EpetraLevelScheduler.hpp> 

#include <Tpetra_MultiVector.hpp>




#include <Kokkos_MultiVector.hpp>
#include <Kokkos_CrsMatrix.hpp>
#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultArithmetic.hpp"
#include "Kokkos_DefaultSparseMultiply.hpp"
#include "Kokkos_Version.hpp"
#include "Kokkos_SerialNode.hpp"
#ifdef HAVE_KOKKOSCLASSIC_TBB
#include "Kokkos_TBBNode.hpp"
#endif
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
#include "Kokkos_TPINode.hpp"
#endif
#ifdef HAVE_KOKKOSCLASSIC_THRUST
#include "Kokkos_ThrustGPUNode.hpp"
#endif

//typedef 


// TODO: only one parallel compute buffer; triangular solve/apply should be in situ

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// class declaration for Epetra_LevelSolver
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class Node> 
class Epetra_LevelSolver : public virtual Epetra_Operator, public virtual Epetra_Object {
  public:

  typedef Tpetra::MultiVector<double,int,int,Node> TMV;
  typedef Tpetra::Map<int, int, Node> tpetraMap;


  //Epetra_LevelSolver(const Epetra_Map &Map, Node &node);
  Epetra_LevelSolver(const Epetra_Map &Map, Teuchos::RCP<Kokkos::TPINode>& node);

  virtual ~Epetra_LevelSolver();
    int Analyze(const Epetra_CrsGraph &G);
    int SetLevelInfo(int numLevels, const int *lsizes, const int *P);
    virtual int Setup(const Epetra_CrsMatrix &L);
    // Methods from Epetra_Operator
    int SetUseTranspose (bool UseTranspose);
    int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const {std::cout << "Epetra_LevelSolver::Apply not implemented" << std::endl; return 1;}          // affect the solve
    int Apply(const TMV &X, TMV &Y) const;          // affect the solve

    //int Epetra_LevelSolver<Node>::Apply(const Tpetra::MultiVector<double,int,int,Node> &Y, 
    //                                Epetra_MultiVector<double,int,int,Node> &X) const 



    int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;   // affect the multiplication
    double NormInf() const;
    const char *Label() const;
    bool UseTranspose() const;
    bool HasNormInf() const;
    const Epetra_Comm &Comm() const;
    const Epetra_Map & OperatorDomainMap() const;
    const Epetra_Map & OperatorRangeMap() const;
    void Print(ostream &os, int verbosity=1) const;
    void setUnitDiag(bool ud);
    bool getUnitDiag() const;
    void setIgnorePerm(bool ip);
    bool getIgnorePerm() const;

    const Teuchos::RCP<const tpetraMap > &getTpetraMap() const {return tmap_RCP_;}

    Teuchos::RCP<Kokkos::TPINode> & getNode() const {return node_;}

  protected:
    //typedef typename Node::template buffer<double>::buffer_t  dbuf;
    //    typedef typename Node::template buffer<int>::buffer_t  obuf;
    //    typedef typename Node::template buffer<Kokkos::size_type>::buffer_t stbuf;

    // constant elements of identity
    Teuchos::RCP<Kokkos::TPINode> &node_;
      //Node &node_;
    const Epetra_Map &map_;
    Teuchos::RCP< const tpetraMap > tmap_RCP_;
    Teuchos::RCP< const tpetraMap > tmapPerm_RCP_;


    // accounting
    bool setupCalled_, unitDiag_, ignorePerm_;
    int numRows_, numLevels_;
    //mutable Kokkos::size_type buf_alloc_len_;

    // computational attributes
    Teuchos::ArrayRCP<int> pinds_;
    Teuchos::Array<int> lsizes_;
    Teuchos::Array<double> Dblock_;
    Teuchos::Array<double> DInvblock_;
    Teuchos::ArrayRCP<double> dBuffOrig;
    Teuchos::ArrayRCP<double> dInvBuffOrig;

    Teuchos::Array<Teuchos::RCP<Kokkos::DefaultSparseMultiply<double,int,Node> > > BblockOps_;
    Teuchos::RCP<Epetra_Import> importer_;
    Teuchos::RCP<Tpetra::Import< int, int, Node > >Timporter_;

    // workspace
    mutable Teuchos::RCP<TMV> importMV_;
    double meanLsize_, stddevLsize_, meanLnnz_, stddevLnnz_;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// implementation for Epetra_LevelSolver
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class Node>
//Epetra_LevelSolver<Node>::Epetra_LevelSolver(const Epetra_Map &Map, Node &node)
Epetra_LevelSolver<Node>::Epetra_LevelSolver(const Epetra_Map &Map, Teuchos::RCP<Kokkos::TPINode> &node)
: node_(node), map_(Map),setupCalled_(false), unitDiag_(false), ignorePerm_(false), numRows_(Map.NumMyPoints()), numLevels_(0), //buf_alloc_len_(0),
  dBuffOrig(Teuchos::null), dInvBuffOrig(Teuchos::null), meanLsize_(-1), stddevLsize_(-1), meanLnnz_(-1), stddevLnnz_(-1)
{
  assert(map_.IndexBase() == 0);
}

template <class Node>
Epetra_LevelSolver<Node>::~Epetra_LevelSolver() {}

template <class Node> const char *Epetra_LevelSolver<Node>::Label() const { return(Epetra_Object::Label()); }
template <class Node> bool Epetra_LevelSolver<Node>::HasNormInf() const { return false; }
template <class Node> double Epetra_LevelSolver<Node>::NormInf() const { return 0.0; }
template <class Node> void Epetra_LevelSolver<Node>::setUnitDiag(bool ud) {unitDiag_ = ud;}
template <class Node> bool Epetra_LevelSolver<Node>::getUnitDiag() const  {return unitDiag_;}
template <class Node> void Epetra_LevelSolver<Node>::setIgnorePerm(bool ip) {ignorePerm_ = ip;}
template <class Node> bool Epetra_LevelSolver<Node>::getIgnorePerm() const  {return ignorePerm_;}

// no support for transpose operations right now, in Apply() or ApplyInverse()
template <class Node> int Epetra_LevelSolver<Node>::SetUseTranspose (bool UseTranspose) { return -1; }
template <class Node> bool Epetra_LevelSolver<Node>::UseTranspose() const { return false; }

template <class Node> const Epetra_Comm & Epetra_LevelSolver<Node>::Comm() const { return map_.Comm(); }
template <class Node> const Epetra_Map & Epetra_LevelSolver<Node>::OperatorDomainMap() const { return map_; }
template <class Node> const Epetra_Map & Epetra_LevelSolver<Node>::OperatorRangeMap() const { return map_; }

template <class Node>
int Epetra_LevelSolver<Node>::SetLevelInfo(int numLevels, const int *lsizes, const int *P)
{
  if (setupCalled_ || pinds_ != Teuchos::null) EPETRA_CHK_ERR(-1);
  if (numLevels < 0 || numLevels > numRows_) EPETRA_CHK_ERR(-2);
  numLevels_ = numLevels;
  if (numRows_) pinds_ = Teuchos::arcp<int>(numRows_);
  for (int i=0; i<numRows_; ++i) 
  {
    if (P[i] < map_.MinLID() || P[i] > map_.MaxLID()) EPETRA_CHK_ERR(-3);
    pinds_[i] = P[i];
  }
  if (numLevels > 0) lsizes_.resize(numLevels);
  int lsizsum = 0;
  for (int i=0; i<numLevels; ++i) 
  {
    if (lsizes[i] < 0) EPETRA_CHK_ERR(-3);
    lsizes_[i] = lsizes[i];
    lsizsum   += lsizes[i];
  }
  if (lsizsum != numRows_) EPETRA_CHK_ERR(-4);
  meanLsize_ = 0.0;
  stddevLsize_ = 0.0;
  for (int i=0; i<numLevels_; ++i) 
  {
    meanLsize_ += lsizes_[i];
  }
  meanLsize_ = meanLsize_ / (double)(numLevels_);
  for (int i=0; i<numLevels_; ++i) 
  {
    double tmp = lsizes_[i] - meanLsize_;
    stddevLsize_ += tmp*tmp;
  }
  stddevLsize_ = sqrt(stddevLsize_ / (double)(numLevels_));
  return 0;
}

template <class Node>
int Epetra_LevelSolver<Node>::Analyze(const Epetra_CrsGraph &g)
{
   if (g.NumMyRows() < 1)
   {
      return SetLevelInfo(0, NULL, NULL); 
   } 

   Teuchos::RCP<const Epetra_CrsGraph> graph = Teuchos::rcp(&g); 

   Isorropia::Epetra::LevelScheduler level(graph); 

   graph.release(); 

   int nlevels = level.numLevels(); 

   int *levelSize = new int [nlevels]; 
   int **rows = new int * [nlevels]; 
   int *permutation = new int [g.NumMyRows()];

   if (!levelSize || !rows || !permutation)
   {
      return 1;
   }

  for (int i=0; i < nlevels; i++)
  {
    levelSize[i] = level.numElemsWithLevel(i);
    rows[i] = new int [levelSize[i]];
    if (!rows[i]){
      return 1;
    }
    level.elemsWithLevel(i, rows[i], levelSize[i]);
  }

  // Create map from row r to row rNew, the row which will replace it
  ///
  int nextRow = 0;
  for (int i=0; i < nlevels; i++)
  {
    for (int j=0; j < levelSize[i]; j++)
    {
      int rowID = rows[i][j];  // a row in level i
      permutation[nextRow++] = rowID;
    }
  }

  int failed = SetLevelInfo(nlevels, levelSize, permutation);

  for (int i=0; i < nlevels; i++)
  {
    if(rows[i]!=0)
    {
      delete [] rows[i];
    }
  }
  if(rows!=0)
  {
    delete [] rows;
  }
  if(levelSize!=0)
  {
    delete [] levelSize;
  }
  if(permutation!=0)
  {
    delete [] permutation;
  }

  return failed;
}

template <class Node>
int Epetra_LevelSolver<Node>::Setup(const Epetra_CrsMatrix &L) 
{
  // need MyCols() and MyRows() to be the same; will assume they are in the same order as well, for now
  if (setupCalled_ == true)                             EPETRA_CHK_ERR(-1);
  if (L.RowMap().SameAs(L.ColMap()) != true)            EPETRA_CHK_ERR(-1);
  if (L.RowMap().SameAs(L.OperatorDomainMap()) != true) EPETRA_CHK_ERR(-1);
  if (L.RowMap().SameAs(L.OperatorRangeMap()) != true)  EPETRA_CHK_ERR(-1);
  // assume this for now
  if (L.RowMap().LinearMap() == false)                  EPETRA_CHK_ERR(-2);
  // necessary, if the matrix is invertible
  if (L.NumMyDiagonals() != L.NumMyRows())              EPETRA_CHK_ERR(-3);
  if (L.LowerTriangular() == false)                     EPETRA_CHK_ERR(-4);

  // permutation indices were given as LIDs; for creating pmap, they will need to 
  // be changed into GIDs

  // setup pmap: this is a permutation of our domain/range map. 
  // it will be our row map and our domain map, used to construct the import and export objects.
  for (int i=0; i<numRows_; ++i) 
  {
    pinds_[i] = map_.GID(pinds_[i]);
    if (pinds_[i] == -1)                                EPETRA_CHK_ERR(-5);
  }
  Epetra_Map pmap(map_.NumGlobalPoints(),numRows_,&*pinds_,0,map_.Comm());

  // setup Import and Export objects
  importer_ = Teuchos::rcp(new Epetra_Import(pmap,map_)); //pmap is targ, map_ is src


  // debugging: these objects should only be doing permutations
  if (importer_->NumRemoteIDs() != 0 || importer_->NumExportIDs() != 0 ) EPETRA_CHK_ERR(-6);
  // will need pmap later, to create import and export multivectors; can get it from the importer_

  if (numRows_ > 0) 
  {
    // get diagonal, then permute it using importer_
    Epetra_Vector PDvec(pmap,false),
                   DVec(map_,false);
    L.ExtractDiagonalCopy(DVec);
    EPETRA_CHK_ERR( PDvec.Import(DVec,*importer_,Insert) );
    // copy to compute buffer
    Dblock_.resize(numRows_); 
    std::copy(&PDvec[0],&PDvec[0]+numRows_,Dblock_.begin());
  }
  // Calculate and store reciprocal of diagonal, perhaps this is a bad optimization due to extra space
  DInvblock_.resize(numRows_);
  std::copy(Dblock_.begin(),Dblock_.end(),DInvblock_.begin());
  for(int i=0;i<numRows_;i++)
  {
    DInvblock_[i] = 1.0 / DInvblock_[i];
  }

  dBuffOrig = Teuchos::arcp(const_cast<double *>(&Dblock_[0]),0,numRows_,false); 
  dInvBuffOrig = Teuchos::arcp(const_cast<double *>(&DInvblock_[0]),0,numRows_,false); 

  // determine number of non-zeros for each (permuted) row, neglecting the diagonal
  Teuchos::Array<size_t> NEPR(numRows_);
  for (int r=0; r<numRows_; ++r) 
  {
    NEPR[r] = L.NumGlobalEntries(pinds_[r]) - 1;
  }
  //meanLnnz_ = std::accumulate(NEPR.begin()+lsizes_[0],NEPR.end(),0) / (double)(numLevels_-1);

  // create array of pointers to hold the Vector and TPICrsMatrix objects we will create below
  if (numLevels_ > 1) BblockOps_.resize(numLevels_-1);

  typedef Kokkos::MultiVector<double,Node> MV;
  typedef Kokkos::CrsMatrix<double,Node> MAT;
  typedef Kokkos::CrsGraph<int,Node> GRPH;

  //typedef Kokkos::DefaultSparseMultiply<MAT,MV> MATVEC;
typedef Kokkos::DefaultSparseMultiply<double,int> MATVEC;

  Epetra_SerialComm scomm;  // no communication necessary, no communication allowed
  // make this simple: create a single column map for all serial TPICrsMatrix objects
  Epetra_LocalMap CMap(L.NumMyRows(),0,scomm);
  int loffset = lsizes_[0];
  // we need to modify column indices, so create some storage to hold them
  const int MaxNumEntries = L.MaxNumEntries();
  Teuchos::Array<double> tmpvals(MaxNumEntries);
  Teuchos::Array<   int> tmpinds(MaxNumEntries);
  stddevLnnz_ = 0.0;

  int numRowsInLevel;
  int nnzInLevel;
  for (int i=1; i<numLevels_; ++i) 
  {
    numRowsInLevel = lsizes_[i];

    // setup B block
    MAT B(numRowsInLevel,node_);
    GRPH G(numRowsInLevel,node_);

    nnzInLevel=0;  // perhaps more efficient way to calculate
    for(int r=0;r<numRowsInLevel;r++)
    {
      nnzInLevel += NEPR[loffset+r];
    }

   if(nnzInLevel!=0)
   {
    Teuchos::ArrayRCP<size_t> offsets = node_->template allocBuffer<size_t> (numRowsInLevel+1);
    Teuchos::ArrayRCP<int>   inds = node_->template allocBuffer<int>(nnzInLevel);
    Teuchos::ArrayRCP<double>    vals = node_->template allocBuffer<double>(nnzInLevel);

    // hosts view of data
    Teuchos::ArrayRCP<size_t>  offsets_h = node_->template viewBufferNonConst<size_t> (Kokkos::WriteOnly,numRowsInLevel+1,offsets);
    Teuchos::ArrayRCP<int>    inds_h = node_->template viewBufferNonConst<int>(Kokkos::WriteOnly,nnzInLevel,inds);
    Teuchos::ArrayRCP<double>     vals_h = node_->template viewBufferNonConst<double>(Kokkos::WriteOnly,nnzInLevel,vals);

    //B.initializeProfile(lsizes_[i],&NEPR[loffset]);
    //    double tmp = B.getNumEntries() - meanLnnz_;
    //    stddevLnnz_ += tmp*tmp;
    // fill it
    int startIndx=0;

    offsets_h[startIndx]=0;
    for (int r=0; r<lsizes_[i]; ++r) 
    {
      // also check that the last element on the row is the diagonal
      int numentries;
      int GID = pinds_[loffset+r];
      EPETRA_CHK_ERR( L.ExtractGlobalRowCopy(GID,MaxNumEntries,numentries,&tmpvals[0],&tmpinds[0]) );
      if (tmpinds[numentries-1] != GID) EPETRA_CHK_ERR(-8);
      if (numentries > 1) 
      {
        // permute column indices and transform to local
        for (int c=0; c<numentries-1; ++c) 
        {
          inds_h[startIndx+c] = pmap.LID(tmpinds[c]);
          vals_h[startIndx+c] = tmpvals[c];
        }
        //B.insertEntries(r,numentries-1,&tmpinds[0],&tmpvals[0]);

        startIndx += numentries-1;
      }
      offsets_h[r+1]=startIndx;      
    }

    offsets_h = Teuchos::null;
    inds_h = Teuchos::null;
    vals_h = Teuchos::null;

    G.setPackedStructure(offsets,inds);
    B.setPackedValues(vals);

    // pass matrix to matvec object
    BblockOps_[i-1] = Teuchos::rcp(new MATVEC(node_));
    BblockOps_[i-1]->initializeStructure(G,Teuchos::View);
    BblockOps_[i-1]->initializeValues(B,Teuchos::View);
    loffset += lsizes_[i];

   }
   else
   {
    BblockOps_[i-1] = Teuchos::null;
    loffset += lsizes_[i];
   }

  }
  //  stddevLnnz_ = sqrt(stddevLnnz_ / (double)(numLevels_-1));

  ///////////////////////////////////////////////////////////////////
  // Build Tpetra map from Epetra Map (source map for Importer)
  ///////////////////////////////////////////////////////////////////
  typedef Tpetra::Map<int, int, Node> tpetraMap;

  Teuchos::RCP< const Teuchos::Comm< int > > TCommRCP = Teuchos::rcp(new Teuchos::SerialComm<int>());

  int *tmpElementList = new int[importer_->SourceMap().NumGlobalElements()];
  importer_->SourceMap().MyGlobalElements(tmpElementList);
  Teuchos::ArrayView<int> aview(tmpElementList,importer_->SourceMap().NumGlobalElements());

  tmap_RCP_ = Teuchos::rcp(new tpetraMap(importer_->SourceMap().NumGlobalElements(), aview, 0, TCommRCP,
                                                getNode()));

  if(tmpElementList!=0)
  {
    delete [] tmpElementList; tmpElementList=0;  // not sure if this is good
  }
  ///////////////////////////////////////////////////////////////////                                                                          

  ///////////////////////////////////////////////////////////////////
  // Build permuted Tpetra map from Epetra Map (target map for importer)
  ///////////////////////////////////////////////////////////////////
  tmpElementList = new int[importer_->TargetMap().NumGlobalElements()];

  importer_->TargetMap().MyGlobalElements(tmpElementList);

  Teuchos::ArrayView<int> aview2(tmpElementList,importer_->TargetMap().NumGlobalElements());

  tmapPerm_RCP_ = Teuchos::rcp(new tpetraMap(importer_->TargetMap().NumGlobalElements(), aview2, 0, TCommRCP, 
 								 node_ )); 

  if(tmpElementList!=0)
  {
    delete [] tmpElementList; tmpElementList=0;  // not sure if this is good
  }
  ///////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  // Build tpetra importer from epetra importer
  //////////////////////////////////////////////////////////////////
  Timporter_ = Teuchos::rcp(new Tpetra::Import< int, int, Node >(tmap_RCP_,tmapPerm_RCP_));
  ///////////////////////////////////////////////////////////////////

  // don't need pinds_ anymore; stored in pmap (which is stored in importer_)
  pinds_ = Teuchos::null;
  // delete local storage
  setupCalled_ = true;
  return 0;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <class Node>
int Epetra_LevelSolver<Node>::Apply(const Tpetra::MultiVector<double,int,int,Node> &Y, 
                                    Tpetra::MultiVector<double,int,int,Node> &X) const 
{
  typedef Kokkos::MultiVector<double,Node>   KMV;
  typedef Kokkos::DefaultArithmetic<KMV> DMVA;
  if (numLevels_ == 0) return 0;
  // solve L*X = Y into X
  using std::endl;
  using std::cout;
  using Teuchos::ArrayRCP;

  const size_t NumVectors = X.getNumVectors();
  TEUCHOS_TEST_FOR_EXCEPT(NumVectors != 1); // FINISH: sparse mat-vec doesn't support this yet
  if (NumVectors != Y.getNumVectors()) EPETRA_CHK_ERR(-1);
  //
  // if we require permutation, then we must have importMV
  // this is because permutation cannot be done in situ, and we are only allowed to write to X
  //
  // if we don't require permutation, we may still require importMV (for temp storage)
  // this is because the Kokkos multivector must be strided.
  // therefore, if X isn't strided, we must use some other writeable strided storage.
  // since we can't write to Y, this means temp storage.
  // 

   const Kokkos::MultiVector<double,Node> *lclX = &(X.getLocalMV());
   const Kokkos::MultiVector<double,Node> *lclY = &(Y.getLocalMV());
   if(lclX==lclY)
   {
     std::cout << "Not handling case where X=Y Yet!!!" << std::endl;
     return 0;
   }

   bool Xstrided = X.isConstantStride();
   bool needTmp  = (!Xstrided) || (!ignorePerm_);

   KMV *Xbptr = 0;

   int MVstride;
   if (needTmp) 
   {                        
     if (importMV_ != Teuchos::null && (importMV_->getNumVectors() != NumVectors)) // delete it if it is the wrong size 
     { 
       importMV_ = Teuchos::null; 
     } 
     if (importMV_ == Teuchos::null) // allocate it if it doesn't exist 
     {   
        importMV_ = Teuchos::rcp(new TMV(tmapPerm_RCP_,NumVectors,false)); 
     } 
     if (!ignorePerm_) 
     { 
       importMV_->doImport(Y,*Timporter_,Tpetra::INSERT);     // we needed importMV to permute 
     } 
     else 
     { 
       (*importMV_) = Y;                                           // we needed importMV for strided workspace 
     } 

     MVstride = importMV_->getStride(); 
     Xbptr = &(importMV_->getLocalMVNonConst()); // not sure if necessary
   } 
   else 
   { 
     X = Y;         // Need to put Y data into X 

     MVstride = X.getStride(); 
     Xbptr = &(X.getLocalMVNonConst());
   } 

   KMV & Xb = *Xbptr;

  /* Example with four levels:
     [D1             ] [X1] = [Y1]
     [B2  D2         ] [X2]   [Y2]
     [   B3   D3     ] [X3]   [Y3]
     [       B4    D4] [X4]   [Y4]

     X1 = D1 \ Y1
     X2 = D2 \ (Y2 - B2*[X1]      )
     X3 = D3 \ (Y3 - B3*[X1;X2]   )
     X4 = D4 \ (Y4 - B4*[X1;X2;X3])

     The B blocks are CrsMatrix objects; the D blocks are Vector objects.
     There is one less B block than there are levels; index into this array is
     therefore decremented by one.

     The Xi,Yi are created using the pointers in Xptrs,Yptrs, which are updated at every step.
     The [X1:Xi-1] uses the pointers from importMV_
     The maps for creating all of these are local maps with serial communicators, created when we
     created the B and D blocks.
     */

   KMV Xd(node_), D(node_),Dinv(node_);

   ///////////////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////
   Teuchos::ArrayRCP<double> offXbufOrig = Xb.getValuesNonConst();

   int Dloffset =0;
   Teuchos::ArrayRCP<double> offXbuf = offXbufOrig.persistingView(Dloffset,lsizes_[0]);
   Teuchos::ArrayRCP<double> offDbuf = dBuffOrig.persistingView(Dloffset,lsizes_[0]);
   Teuchos::ArrayRCP<double> offDInvbuf = dInvBuffOrig.persistingView(Dloffset,lsizes_[0]);
   ///////////////////////////////////////////////////////////////////////

   ///////////////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////
   for (int i=0; i<numLevels_; ++i)  
   { 
     // point the vector/multivectors to the current diagonal blocks 
     D.initializeValues(lsizes_[i],1,offDbuf,lsizes_[i]); 
     Xd.initializeValues(lsizes_[i],NumVectors,offXbuf,MVstride); 

     if (i != 0 && BblockOps_[i-1] != Teuchos::null)  
     { 
       BblockOps_[i-1]->multiply(Teuchos::NO_TRANS,-1.0,Xb,1.0,Xd);      // Xd -= B*Xb 
     } 

     // scale Xd by diagonal block 
     if (unitDiag_ == false) 
     { 
       //DMVA::Divide(Xd,(const KMV&)D);         /// Assumes ones on diagonal for now
       // Perhaps adding a specific function to Kokkos would make this more efficient
       //KMV recipD(D);
       Dinv.initializeValues(lsizes_[i],1,offDInvbuf,lsizes_[i]); 
       //DMVA::Recip(Dinv,D);                          // recipD = 1 ./ D  -- Should really do this once in setup
       DMVA::ElemMult(Xd,0,1.0,Xd,Dinv);             // Xd(i,j) = 0 * Xd(i,j) + 1.0 * Xd(i,j) * recipD(i,1) = Xd(i,j) * recipD(i,1)
     } 

     if(i!=numLevels_-1)
     {
       // increment the pointers for next Xd and D 
       Dloffset += lsizes_[i];
       offXbuf = offXbufOrig.persistingView(Dloffset,lsizes_[i+1]);
       offDbuf = dBuffOrig.persistingView(Dloffset,lsizes_[i+1]);
       offDInvbuf = dInvBuffOrig.persistingView(Dloffset,lsizes_[i+1]);
     }

   }
   ///////////////////////////////////////////////////////////////////////

  if (needTmp)
  {
    if (!ignorePerm_)
    {
      X.doExport(*importMV_,*Timporter_,Tpetra::INSERT);
    }
    else {
      X = (*importMV_);
    }
  }

  return 0;
}
////////////////////////////////////////////////////////////////////////////////



/* template <class Node> */
/* int Epetra_LevelSolver<Node>::Apply(const Epetra_MultiVector &Y, Epetra_MultiVector &X) const  */
/* { */
/*   typedef Kokkos::MultiVector<double,Node>   MV; */
/*   typedef Kokkos::DefaultArithmetic<MV>        DMVA; */
/*   if (numLevels_ == 0) return 0; */
/*   // solve L*X = Y into X */
/*   using std::endl; */
/*   using std::cout; */
/*   using Teuchos::ArrayRCP; */

/* /\*   using Kokkos::MultiVector; *\/ */
/* /\*   using Kokkos::CrsMatrix; *\/ */
/* /\*   using Kokkos::CrsGraph; *\/ */
/* /\*   using Kokkos::DefaultArithmetic; *\/ */
/* /\*   using Kokkos::DefaultSparseMultiply; *\/ */
/* /\*   using Kokkos::SerialNode; *\/ */
/* /\*   using Teuchos::RCP; *\/ */
/* /\*   using Teuchos::rcp; *\/ */
/* /\*   using Teuchos::null; *\/ */


/*   const int NumVectors = X.NumVectors(); */
/*   TEUCHOS_TEST_FOR_EXCEPT(NumVectors != 1); // FINISH: sparse mat-vec doesn't support this yet */
/*   if (NumVectors != Y.NumVectors()) EPETRA_CHK_ERR(-1); */
/*   // */
/*   // if we require permutation, then we must have importMV */
/*   // this is because permutation cannot be done in situ, and we are only allowed to write to X */
/*   // */
/*   // if we don't require permutation, we may still require importMV (for temp storage) */
/*   // this is because the Kokkos multivector must be strided. */
/*   // therefore, if X isn't strided, we must use some other writeable strided storage. */
/*   // since we can't write to Y, this means temp storage. */
/*   //  */
/*   bool Xstrided = X.ConstantStride(); */
/*   bool XeqY     = (X.Values() == Y.Values()) && X.ConstantStride() && Y.ConstantStride(); */
/*   bool needTmp  = (!Xstrided) || (!ignorePerm_); */


/*   //dbuf MVbuffer; */
/*   ArrayRCP<double> MVbuffer; */
/*   // */

/*   //  MVbuffer = node_->template allocBuffer<double>(Y.GlobalLength() * Y.NumVectors()); */

/*    int MVstride; */
/*    if (needTmp) // delete it if it is the wrong size  */
/*    {                         */
/*      if (importMV_ != Teuchos::null && (importMV_->NumVectors() != NumVectors))  */
/*      {  */
/*        importMV_ = Teuchos::null;  */
/*      }  */
/* /\*     if (importMV_ == Teuchos::null) {   // allocate it if it doesn't exist *\/ */
/* /\*       importMV_ = Teuchos::rcp(new Epetra_MultiVector(importer_->TargetMap(),NumVectors,false)); *\/ */
/* /\*     } *\/ */
/* /\*     if (!ignorePerm_) { *\/ */
/* /\*       EPETRA_CHK_ERR(importMV_->Import(Y,*importer_,Insert));     // we needed importMV to permute *\/ */
/* /\*     } *\/ */
/* /\*     else { *\/ */
/* /\*       (*importMV_) = Y;                                           // we needed importMV for strided workspace *\/ */
/* /\*     } *\/ */
/* /\*     MVbuffer = importMV_->Values(); *\/ */
/* /\*     MVstride = importMV_->Stride(); *\/ */
/*    }  */
/*    else  */
/*    {  */
/*      if (!XeqY) X = Y;         // Need to put Y data into X  */
/* /\*     MVbuffer = X.Values();    // X provides our strided workspace.  *\/ */
/*      MVstride = X.Stride();  */
/*    }  */
/*   /\* Example with four levels: */
/*      [D1             ] [X1] = [Y1] */
/*      [B2  D2         ] [X2]   [Y2] */
/*      [   B3   D3     ] [X3]   [Y3] */
/*      [       B4    D4] [X4]   [Y4] */

/*      X1 = D1 \ Y1 */
/*      X2 = D2 \ (Y2 - B2*[X1]      ) */
/*      X3 = D3 \ (Y3 - B3*[X1;X2]   ) */
/*      X4 = D4 \ (Y4 - B4*[X1;X2;X3]) */

/*      The B blocks are CrsMatrix objects; the D blocks are Vector objects. */
/*      There is one less B block than there are levels; index into this array is */
/*      therefore decremented by one. */

/*      The Xi,Yi are created using the pointers in Xptrs,Yptrs, which are updated at every step. */
/*      The [X1:Xi-1] uses the pointers from importMV_ */
/*      The maps for creating all of these are local maps with serial communicators, created when we */
/*      created the B and D blocks. */
/*      *\/ */

/* /\*   MV Xb(node_), Xd(node_), D(node_); *\/ */
/* /\*   Xb.initializeValues(numRows_,NumVectors,MVbuffer,MVstride); *\/ */




/* /\*   dbuf offXbuf = MVbuffer, *\/ */
/* /\*        offDbuf = const_cast<double *>(&Dblock_[0]); *\/ */
/* /\*   for (int i=0; i<numLevels_; ++i)  *\/ */
/* /\*   { *\/ */

/* /\*     // point the vector/multivectors to the current diagonal blocks *\/ */
/* /\*     D.initializeValues(lsizes_[i],1,offDbuf,lsizes_[i]); *\/ */
/* /\*     Xd.initializeValues(lsizes_[i],NumVectors,offXbuf,MVstride); *\/ */

/* /\*     if (i != 0)  *\/ */
/* /\*     { *\/ */
/* /\*       BblockOps_[i-1]->Apply(false,-1.0,Xb,1.0,Xd);             // Xd -= B*Xb *\/ */
/* /\*     } *\/ */

/* /\*     // scale Xd by diagonal block *\/ */
/* /\*     if (unitDiag_ == false) { *\/ */
/* /\*       DMVA::Divide(Xd,(const MV&)D); *\/ */
/* /\*     } *\/ */
/* /\*     // increment the pointers for next Xd and D *\/ */
/* /\*     offXbuf = offXbuf + lsizes_[i]; *\/ */
/* /\*     offDbuf = offDbuf + lsizes_[i]; *\/ */
/* /\*   } *\/ */


/* /\*   if (needTmp)  *\/ */
/* /\*   { *\/ */
/* /\*     if (!ignorePerm_)  *\/ */
/* /\*     { *\/ */
/* /\*       EPETRA_CHK_ERR(X.Export(*importMV_,*importer_,Insert)); *\/ */
/* /\*     } *\/ */
/* /\*     else { *\/ */
/* /\*       X = (*importMV_); *\/ */
/* /\*     } *\/ */
/* /\*   } *\/ */
/*   return 0; */
/* } */

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <class Node> 
int Epetra_LevelSolver<Node>::ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const  
{
  std::cout << "Epetra_LevelSolver::ApplyInverse currently not implemented" << std::endl;
 
/*   typedef Kokkos::MultiVector<double,Node>   MV; */
/*   typedef Kokkos::DefaultArithmetic<MV>        DMVA; */
/*   if (numLevels_ == 0) return 0; */
/*   // apply L*X into Y */
/*   using std::endl; */
/*   using std::cout; */
/*   const int NumVectors = X.NumVectors(); */
/*   if (NumVectors != Y.NumVectors()) EPETRA_CHK_ERR(-1); */
/*   TEUCHOS_TEST_FOR_EXCEPT(NumVectors != 1); // FINISH: sparse mat-vec doesn't support this yet */
/*   // */
/*   // if we require permutation, then we must have importMV */
/*   // this is because permutation cannot be done in situ, and we are only allowed to write to Y */
/*   // */
/*   // if we don't require permutation, we may still require importMV (for temp storage) */
/*   // this is because the Kokkos multivector must be strided. */
/*   // therefore, if Y isn't strided, we must use some other writeable strided storage. */
/*   // since we can't write to X, this means temp storage. */
/*   //  */
/*   bool Ystrided = Y.ConstantStride(); */
/*   bool XeqY     = (X.Values() == Y.Values()) && X.ConstantStride() && Y.ConstantStride(); */
/*   bool needTmp  = (!Ystrided) || (!ignorePerm_); */
/*   dbuf MVbuffer; */
/*   int MVstride; */
/*   if (needTmp) {                        // delete it if it is the wrong size */
/*     if (importMV_ != Teuchos::null && (importMV_->NumVectors() != NumVectors)) { */
/*       importMV_ == Teuchos::null; */
/*     } */
/*     if (importMV_ == Teuchos::null) {   // allocate it if it doesn't exist */
/*       importMV_ = Teuchos::rcp(new Epetra_MultiVector(importer_->TargetMap(),NumVectors,false)); */
/*     } */
/*     if (!ignorePerm_) { */
/*       EPETRA_CHK_ERR(importMV_->Import(X,*importer_,Insert));     // we needed importMV to permute */
/*     } */
/*     else { */
/*       (*importMV_) = X;                                           // we needed importMV for strided workspace */
/*     } */
/*     MVbuffer = importMV_->Values(); */
/*     MVstride = importMV_->Stride(); */
/*   } */
/*   else { */
/*     if (!XeqY) Y = X;         // Need to put X data into Y */
/*     MVbuffer = Y.Values();    // X provides our strided workspace.  */
/*     MVstride = Y.Stride(); */
/*   } */
/*   /\* Example with four levels: */
/*          [D1             ] [X1] = [Y1] */
/*          [B2  D2         ] [X2]   [Y2] */
/*          [   B3   D3     ] [X3]   [Y3] */
/*          [       B4    D4] [X4]   [Y4] */
      
/*       Y1 = D1*X1 */
/*       Y2 = D2*X2 + B2*[X1] */
/*       Y3 = D3*X3 + B3*[X1;X2] */
/*       Y4 = D4*X4 + B4*[X1;X2;X3] */

/*       To do this in situ, we must iterate backwards, computing Y4 before Y3 ... before Y1. */
/*       Otherwise, the computation of Y1 would overwrite X1, which is needed for Y2:Y4 */
/*    *\/ */
/*   MV Xb(node_), Yd(node_), D(node_); */
/*   Xb.initializeValues(numRows_,NumVectors,MVbuffer,MVstride); */
/*   dbuf offYbuf = MVbuffer+numRows_,     // point to just after last element in first column */
/*        offDbuf = const_cast<double *>(&Dblock_[0]+numRows_); */
/*   for (int i=numLevels_-1; i>=0; --i) { */
/*     // decrement the pointers for previous Xd and Yd */
/*     offYbuf = offYbuf - lsizes_[i]; */
/*     offDbuf = offDbuf - lsizes_[i]; */
/*     // point the vector/multivectors to the current diagonal blocks */
/*     D.initializeValues( lsizes_[i],1         ,offDbuf,lsizes_[i]); */
/*     Yd.initializeValues(lsizes_[i],NumVectors,offYbuf,MVstride); */
/*     if (unitDiag_ == false) { */
/*       DMVA::Multiply(Yd,(const MV&)D);                 // Yd = D*Yd = D*Xd */
/*     } */
/*     if (i != 0) { */
/*       BblockOps_[i-1]->Apply(false,1.0,Xb,1.0,Yd);    // Yd = Yd + B*Xb = D*Xd + B*Xb */
/*     } */
/*   } */
/*   if (needTmp) { */
/*     if (!ignorePerm_) { */
/*       EPETRA_CHK_ERR(Y.Export(*importMV_,*importer_,Insert)); */
/*     } */
/*     else { */
/*       Y = (*importMV_); */
/*     } */
/*   } */
   return 0; 
}
////////////////////////////////////////////////////////////////////////////////


template <class Node>
void Epetra_LevelSolver<Node>::Print(ostream &os, int verbosity) const 
{
  // verbosity:
  // 0  prints nothing
  // 1  prints level statistics
  // 2  prints detailed level info
  using std::endl;
  const Epetra_Comm &comm = map_.Comm();
  for (int p=0; p<comm.NumProc(); ++p) {
    if (p == comm.MyPID()) {
      if (verbosity == 0) return;
      else if (verbosity == 1) {
        os << "PID: " << p << endl;
        os << "Num Levels: " << numLevels_ << endl;
        if (meanLsize_ != -1.0) {
          os << "Mean level size: " << meanLsize_ << endl;
          os << "Std dev level size: " << stddevLsize_ << endl;
        }
	//        if (meanLnnz_ != -1.0) {
        //  os << "Mean level nnz: " << meanLnnz_ << endl;
        //  os << "Std dev level nnz: " << stddevLnnz_ << endl;
        //}
      }
      else if (verbosity == 2) {
        //if (p == 0) {
        //  os << "LocalLevel    Rows     Nonzeros" << endl;
        //}
        for (int i=0; i<numLevels_; ++i) {
          os << i;
          if (lsizes_.size()) os << "\t" << lsizes_[i];
          // if (BblockOps_.size()) {
          //   if (i > 0) os << "\t" << BblockOps_[i-1]->getNumEntries();
          //   else       os << "\t0";
          // }
          os << endl;
        }
      }
    }
  }
  comm.Barrier();
  comm.Barrier();
  comm.Barrier();
}
#endif
