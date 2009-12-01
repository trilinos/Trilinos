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

#include <Kokkos_MultiVector.hpp>
#include <Kokkos_CrsMatrix.hpp>
#include <Kokkos_DefaultSparseMultiply.hpp>
#include <Kokkos_DefaultArithmetic.hpp>
#include <Isorropia_LevelScheduler.hpp> 
#include <Isorropia_EpetraLevelScheduler.hpp> 

// TODO: only one parallel compute buffer; triangular solve/apply should be in situ

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// class declaration for Epetra_LevelSolver
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class Node> 
class Epetra_LevelSolver : public virtual Epetra_Operator, public virtual Epetra_Object {
  public:
    Epetra_LevelSolver(const Epetra_Map &Map, Node &node);
    virtual ~Epetra_LevelSolver();
    int Analyze(const Epetra_CrsGraph &G);
    int SetLevelInfo(int numLevels, const int *lsizes, const int *P);
    virtual int Setup(const Epetra_CrsMatrix &L);
    // Methods from Epetra_Operator
    int SetUseTranspose (bool UseTranspose);
    int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;          // affect the solve
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
  protected:
    typedef typename Node::template buffer<double           >::buffer_t  dbuf;
    typedef typename Node::template buffer<int              >::buffer_t  obuf;
    typedef typename Node::template buffer<Kokkos::size_type>::buffer_t stbuf;

    // constant elements of identity
    Node &node_;
    const Epetra_Map &map_;

    // accounting
    bool setupCalled_, unitDiag_, ignorePerm_;
    int numRows_, numLevels_;
    mutable Kokkos::size_type buf_alloc_len_;

    // computational attributes
    Teuchos::ArrayRCP<int> pinds_;
    Teuchos::Array<int> lsizes_;
    Teuchos::Array<double> Dblock_;
    Teuchos::Array< Teuchos::RCP< 
        Kokkos::DefaultSparseMultiply<Kokkos::CrsMatrix<double,int,Node>, Kokkos::MultiVector<double,int,Node> >
    > > BblockOps_;
    Teuchos::RCP<Epetra_Import> importer_;
    // workspace
    mutable Teuchos::RCP<Epetra_MultiVector> importMV_;
    double meanLsize_, stddevLsize_, meanLnnz_, stddevLnnz_;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// implementation for Epetra_LevelSolver
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class Node>
Epetra_LevelSolver<Node>::Epetra_LevelSolver(const Epetra_Map &Map, Node &node)
  : node_(node), map_(Map), setupCalled_(false), unitDiag_(false), ignorePerm_(false), numRows_(Map.NumMyPoints()), numLevels_(0), buf_alloc_len_(0),
    meanLsize_(-1), stddevLsize_(-1), meanLnnz_(-1), stddevLnnz_(-1)
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
  for (int i=0; i<numRows_; ++i) {
    if (P[i] < map_.MinLID() || P[i] > map_.MaxLID()) EPETRA_CHK_ERR(-3);
    pinds_[i] = P[i];
  }
  if (numLevels > 0) lsizes_.resize(numLevels);
  int lsizsum = 0;
  for (int i=0; i<numLevels; ++i) {
    if (lsizes[i] < 0) EPETRA_CHK_ERR(-3);
    lsizes_[i] = lsizes[i];
    lsizsum   += lsizes[i];
  }
  if (lsizsum != numRows_) EPETRA_CHK_ERR(-4);
  meanLsize_ = 0.0;
  stddevLsize_ = 0.0;
  for (int i=0; i<numLevels_; ++i) {
    meanLsize_ += lsizes_[i];
  }
  meanLsize_ = meanLsize_ / (double)(numLevels_);
  for (int i=0; i<numLevels_; ++i) {
    double tmp = lsizes_[i] - meanLsize_;
    stddevLsize_ += tmp*tmp;
  }
  stddevLsize_ = sqrt(stddevLsize_ / (double)(numLevels_));
  return 0;
}

template <class Node>
int Epetra_LevelSolver<Node>::Analyze(const Epetra_CrsGraph &g)
{
  if (g.NumMyRows() < 1){
    return SetLevelInfo(0, NULL, NULL);
  }

  Teuchos::RCP<const Epetra_CrsGraph> graph = Teuchos::rcp(&g);

  Isorropia::Epetra::LevelScheduler level(graph);

  graph.release();

  int nlevels = level.numLevels();

  int *levelSize = new int [nlevels];
  int **rows = new int * [nlevels];
  int *permutation = new int [g.NumMyRows()];

  if (!levelSize || !rows || !permutation){
    return 1;
  }

  for (int i=0; i < nlevels; i++){
    levelSize[i] = level.numElemsWithLevel(i);
    rows[i] = new int [levelSize[i]];
    if (!rows[i]){
      return 1;
    }
    level.elemsWithLevel(i, rows[i], levelSize[i]);
  }

  // Create map from row r to row rNew, the row which will replace it

  int nextRow = 0;
  for (int i=0; i < nlevels; i++){
    for (int j=0; j < levelSize[i]; j++){
      int rowID = rows[i][j];  // a row in level i
      permutation[nextRow++] = rowID;
    }
  }

  int failed = SetLevelInfo(nlevels, levelSize, permutation);

  for (int i=0; i < nlevels; i++){
    delete [] rows[i];
  }
  delete [] rows;
  delete [] levelSize;
  delete [] permutation;

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
  for (int i=0; i<numRows_; ++i) {
    pinds_[i] = map_.GID(pinds_[i]);
    if (pinds_[i] == -1)                                EPETRA_CHK_ERR(-5);
  }
  Epetra_Map pmap(map_.NumGlobalPoints(),numRows_,&*pinds_,0,map_.Comm());

  // setup Import and Export objects
  importer_ = Teuchos::rcp(new Epetra_Import(pmap,map_));
  // debugging: these objects should only be doing permutations
  if (importer_->NumRemoteIDs() != 0 || importer_->NumExportIDs() != 0 ) EPETRA_CHK_ERR(-6);
  // will need pmap later, to create import and export multivectors; can get it from the importer_

  if (numRows_ > 0) {
    // get diagonal, then permute it using importer_
    Epetra_Vector PDvec(pmap,false),
                   DVec(map_,false);
    L.ExtractDiagonalCopy(DVec);
    EPETRA_CHK_ERR( PDvec.Import(DVec,*importer_,Insert) );
    // copy to compute buffer
    Dblock_.resize(numRows_); 
    std::copy(&PDvec[0],&PDvec[0]+numRows_,Dblock_.begin());
  }

  // determine number of non-zeros for each (permuted) row, neglecting the diagonal
  Teuchos::Array<Kokkos::size_type> NEPR(numRows_);
  for (int r=0; r<numRows_; ++r) {
    NEPR[r] = L.NumGlobalEntries(pinds_[r]) - 1;
  }
  meanLnnz_ = std::accumulate(NEPR.begin()+lsizes_[0],NEPR.end(),0) / (double)(numLevels_-1);

  // create array of pointers to hold the Vector and TPICrsMatrix objects we will create below
  if (numLevels_ > 1) BblockOps_.resize(numLevels_-1);

  typedef Kokkos::MultiVector<double,int,Node> MV;
  typedef Kokkos::CrsMatrix<double,int,Node> MAT;
  typedef Kokkos::DefaultSparseMultiply<MAT,MV> MATVEC;

  Epetra_SerialComm scomm;  // no communication necessary, no communication allowed
  // make this simple: create a single column map for all serial TPICrsMatrix objects
  Epetra_LocalMap CMap(L.NumMyRows(),0,scomm);
  int loffset = lsizes_[0];
  // we need to modify column indices, so create some storage to hold them
  const int MaxNumEntries = L.MaxNumEntries();
  Teuchos::Array<double> vals(MaxNumEntries);
  Teuchos::Array<   int> inds(MaxNumEntries);
  stddevLnnz_ = 0.0;
  for (int i=1; i<numLevels_; ++i) {
    // setup B block
    MAT B(node_);
    B.initializeProfile(lsizes_[i],&NEPR[loffset]);
    double tmp = B.getNumEntries() - meanLnnz_;
    stddevLnnz_ += tmp*tmp;
    // fill it
    for (int r=0; r<lsizes_[i]; ++r) {
      // also check that the last element on the row is the diagonal
      int numentries;
      int GID = pinds_[loffset+r];
      EPETRA_CHK_ERR( L.ExtractGlobalRowCopy(GID,MaxNumEntries,numentries,&vals[0],&inds[0]) );
      if (inds[numentries-1] != GID) EPETRA_CHK_ERR(-8);
      if (numentries > 1) {
        // permute column indices and transform to local
        for (int c=0; c<numentries-1; ++c) {
          inds[c] = pmap.LID(inds[c]);
        }
        B.insertEntries(r,numentries-1,&inds[0],&vals[0]);
      }
    }
    // pass matrix to matvec object
    BblockOps_[i-1] = Teuchos::rcp(new MATVEC(node_));
    BblockOps_[i-1]->initializeStructure(B,false);
    BblockOps_[i-1]->initializeValues(B,false);
    loffset += lsizes_[i];
  }
  stddevLnnz_ = sqrt(stddevLnnz_ / (double)(numLevels_-1));
  // don't need pinds_ anymore; stored in pmap (which is stored in importer_)
  pinds_ = Teuchos::null;
  // delete local storage
  setupCalled_ = true;
  return 0;
}


template <class Node>
int Epetra_LevelSolver<Node>::Apply(const Epetra_MultiVector &Y, Epetra_MultiVector &X) const 
{
  typedef Kokkos::MultiVector<double,int,Node>   MV;
  typedef Kokkos::DefaultArithmetic<MV>        DMVA;
  if (numLevels_ == 0) return 0;
  // solve L*X = Y into X
  using std::endl;
  using std::cout;
  const int NumVectors = X.NumVectors();
  TEST_FOR_EXCEPT(NumVectors != 1); // FINISH: sparse mat-vec doesn't support this yet
  if (NumVectors != Y.NumVectors()) EPETRA_CHK_ERR(-1);
  //
  // if we require permutation, then we must have importMV
  // this is because permutation cannot be done in situ, and we are only allowed to write to X
  //
  // if we don't require permutation, we may still require importMV (for temp storage)
  // this is because the Kokkos multivector must be strided.
  // therefore, if X isn't strided, we must use some other writeable strided storage.
  // since we can't write to Y, this means temp storage.
  // 
  bool Xstrided = X.ConstantStride();
  bool XeqY     = (X.Values() == Y.Values()) && X.ConstantStride() && Y.ConstantStride();
  bool needTmp  = (!Xstrided) || (!ignorePerm_);
  dbuf MVbuffer;
  int MVstride;
  if (needTmp) {                        // delete it if it is the wrong size
    if (importMV_ != Teuchos::null && (importMV_->NumVectors() != NumVectors)) {
      importMV_ == Teuchos::null;
    }
    if (importMV_ == Teuchos::null) {   // allocate it if it doesn't exist
      importMV_ = Teuchos::rcp(new Epetra_MultiVector(importer_->TargetMap(),NumVectors,false));
    }
    if (!ignorePerm_) {
      EPETRA_CHK_ERR(importMV_->Import(Y,*importer_,Insert));     // we needed importMV to permute
    }
    else {
      (*importMV_) = Y;                                           // we needed importMV for strided workspace
    }
    MVbuffer = importMV_->Values();
    MVstride = importMV_->Stride();
  }
  else {
    if (!XeqY) X = Y;         // Need to put Y data into X
    MVbuffer = X.Values();    // X provides our strided workspace. 
    MVstride = X.Stride();
  }
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
  MV Xb(node_), Xd(node_), D(node_);
  Xb.initializeValues(numRows_,NumVectors,MVbuffer,MVstride);
  dbuf offXbuf = MVbuffer,
       offDbuf = const_cast<double *>(&Dblock_[0]);
  for (int i=0; i<numLevels_; ++i) {
    // point the vector/multivectors to the current diagonal blocks
    D.initializeValues(lsizes_[i],1,offDbuf,lsizes_[i]);
    Xd.initializeValues(lsizes_[i],NumVectors,offXbuf,MVstride);
    if (i != 0) {
      BblockOps_[i-1]->Apply(false,-1.0,Xb,1.0,Xd);             // Xd -= B*Xb
    }
    // scale Xd by diagonal block
    if (unitDiag_ == false) {
      DMVA::Divide(Xd,(const MV&)D);
    }
    // increment the pointers for next Xd and D
    offXbuf = offXbuf + lsizes_[i];
    offDbuf = offDbuf + lsizes_[i];
  }
  if (needTmp) {
    if (!ignorePerm_) {
      EPETRA_CHK_ERR(X.Export(*importMV_,*importer_,Insert));
    }
    else {
      X = (*importMV_);
    }
  }
  return 0;
}


template <class Node>
int Epetra_LevelSolver<Node>::ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const 
{
  typedef Kokkos::MultiVector<double,int,Node>   MV;
  typedef Kokkos::DefaultArithmetic<MV>        DMVA;
  if (numLevels_ == 0) return 0;
  // apply L*X into Y
  using std::endl;
  using std::cout;
  const int NumVectors = X.NumVectors();
  if (NumVectors != Y.NumVectors()) EPETRA_CHK_ERR(-1);
  TEST_FOR_EXCEPT(NumVectors != 1); // FINISH: sparse mat-vec doesn't support this yet
  //
  // if we require permutation, then we must have importMV
  // this is because permutation cannot be done in situ, and we are only allowed to write to Y
  //
  // if we don't require permutation, we may still require importMV (for temp storage)
  // this is because the Kokkos multivector must be strided.
  // therefore, if Y isn't strided, we must use some other writeable strided storage.
  // since we can't write to X, this means temp storage.
  // 
  bool Ystrided = Y.ConstantStride();
  bool XeqY     = (X.Values() == Y.Values()) && X.ConstantStride() && Y.ConstantStride();
  bool needTmp  = (!Ystrided) || (!ignorePerm_);
  dbuf MVbuffer;
  int MVstride;
  if (needTmp) {                        // delete it if it is the wrong size
    if (importMV_ != Teuchos::null && (importMV_->NumVectors() != NumVectors)) {
      importMV_ == Teuchos::null;
    }
    if (importMV_ == Teuchos::null) {   // allocate it if it doesn't exist
      importMV_ = Teuchos::rcp(new Epetra_MultiVector(importer_->TargetMap(),NumVectors,false));
    }
    if (!ignorePerm_) {
      EPETRA_CHK_ERR(importMV_->Import(X,*importer_,Insert));     // we needed importMV to permute
    }
    else {
      (*importMV_) = X;                                           // we needed importMV for strided workspace
    }
    MVbuffer = importMV_->Values();
    MVstride = importMV_->Stride();
  }
  else {
    if (!XeqY) Y = X;         // Need to put X data into Y
    MVbuffer = Y.Values();    // X provides our strided workspace. 
    MVstride = Y.Stride();
  }
  /* Example with four levels:
         [D1             ] [X1] = [Y1]
         [B2  D2         ] [X2]   [Y2]
         [   B3   D3     ] [X3]   [Y3]
         [       B4    D4] [X4]   [Y4]
      
      Y1 = D1*X1
      Y2 = D2*X2 + B2*[X1]
      Y3 = D3*X3 + B3*[X1;X2]
      Y4 = D4*X4 + B4*[X1;X2;X3]

      To do this in situ, we must iterate backwards, computing Y4 before Y3 ... before Y1.
      Otherwise, the computation of Y1 would overwrite X1, which is needed for Y2:Y4
   */
  MV Xb(node_), Yd(node_), D(node_);
  Xb.initializeValues(numRows_,NumVectors,MVbuffer,MVstride);
  dbuf offYbuf = MVbuffer+numRows_,     // point to just after last element in first column
       offDbuf = const_cast<double *>(&Dblock_[0]+numRows_);
  for (int i=numLevels_-1; i>=0; --i) {
    // decrement the pointers for previous Xd and Yd
    offYbuf = offYbuf - lsizes_[i];
    offDbuf = offDbuf - lsizes_[i];
    // point the vector/multivectors to the current diagonal blocks
    D.initializeValues( lsizes_[i],1         ,offDbuf,lsizes_[i]);
    Yd.initializeValues(lsizes_[i],NumVectors,offYbuf,MVstride);
    if (unitDiag_ == false) {
      DMVA::Multiply(Yd,(const MV&)D);                 // Yd = D*Yd = D*Xd
    }
    if (i != 0) {
      BblockOps_[i-1]->Apply(false,1.0,Xb,1.0,Yd);    // Yd = Yd + B*Xb = D*Xd + B*Xb
    }
  }
  if (needTmp) {
    if (!ignorePerm_) {
      EPETRA_CHK_ERR(Y.Export(*importMV_,*importer_,Insert));
    }
    else {
      Y = (*importMV_);
    }
  }
  return 0;
}


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
        if (meanLnnz_ != -1.0) {
          os << "Mean level nnz: " << meanLnnz_ << endl;
          os << "Std dev level nnz: " << stddevLnnz_ << endl;
        }
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
