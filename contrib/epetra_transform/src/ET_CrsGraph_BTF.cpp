
#include <vector>

#include <ET_CrsGraph_BTF.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>

extern "C" {
extern void mattrans_( int*, int*, int*, int*, int*, int* );
extern void genbtf_( int*, int*, int*, int*, int*, int*, int*, int*, int*,
                     int*, int*, int*, int*, int*, int*, int*, int*, int*,
                     int*, int*, int* );
}

namespace Epetra_Transform {

std::auto_ptr<Epetra_CrsGraph> CrsGraph_BTF::operator()( const Epetra_CrsGraph & original )
{
  if( original.RowMap().DistributedGlobal() )
  { cout << "FAIL for Global!\n"; exit(0); }
  if( original.IndicesAreGlobal() )
  { cout << "FAIL for Global Indices!\n"; exit(0); }

  int n = original.NumMyRows();
  int nnz = original.NumMyNonzeros();

  //create std CRS format
  vector<int> ia(n+1,0);
  vector<int> ja(nnz);
  int cnt;
  for( int i = 0; i < n; ++i )
  {
    int * tmpP = &ja[ia[i]];
    original.ExtractMyRowCopy( i, nnz-ia[i], cnt, tmpP );
    ia[i+1] = ia[i] + cnt;
  }

  //convert to Fortran indexing
  for( int i = 0; i < n+1; ++i ) ++ia[i];
  for( int i = 0; i < nnz; ++i ) ++ja[i];

#ifdef BTF_VERBOSE
  cout << "-----------------------------------------\n";
  cout << "CRS Format Graph\n";
  cout << "-----------------------------------------\n";
  for( int i = 0; i < n; ++i )
  {
    cout << i << ": " << ia[i+1] << ": ";
    for( int j = ia[i]-1; j<ia[i+1]-1; ++j )
      cout << " " << ja[j];
    cout << endl;
  }
  cout << "-----------------------------------------\n";
#endif

  vector<int> iat(n+1);
  vector<int> jat(nnz);
  int * jaf = &ja[0];
  int * iaf = &ia[0];
  int * jatf = &jat[0];
  int * iatf = &iat[0];
  mattrans_( &n, &n, jaf, iaf, jatf, iatf );
    
#ifdef BTF_VERBOSE
  cout << "-----------------------------------------\n";
  cout << "CCS Format Graph\n";
  cout << "-----------------------------------------\n";
  for( int i = 0; i < n; ++i )
  {
    cout << i << ": " << iat[i+1] << ": ";
    for( int j = iat[i]-1; j<iat[i+1]-1; ++j )
      cout << " " << jat[j];
    cout << endl;
  }
  cout << "-----------------------------------------\n";
#endif

  vector<int> w(10*n);

  vector<int> rowperm(n);
  vector<int> colperm(n);

  //horizontal block
  int nhrows, nhcols, hrzcmp;
  //square block
  int nsrows, sqcmpn;
  //vertial block
  int nvrows, nvcols, vrtcmp;

  vector<int> rcmstr(n+1);
  vector<int> ccmstr(n+1);

  int msglvl = 0;
  int output = 6;

  genbtf_( &n, &n, &iat[0], &jat[0], &ia[0], &ja[0], &w[0],
          &rowperm[0], &colperm[0], &nhrows, &nhcols,
          &hrzcmp, &nsrows, &sqcmpn, &nvrows, &nvcols, &vrtcmp,
          &rcmstr[0], &ccmstr[0], &msglvl, &output );

  //convert back to C indexing
  for( int i = 0; i < n; ++i )
  {
    --rowperm[i];
    --colperm[i];
  }
  for( int i = 0; (i<n+1) && (rcmstr[i]!=n+1); ++i )
  {
    --rcmstr[i];
    --ccmstr[i];
  }

#ifdef BTF_VERBOSE
  cout << "-----------------------------------------\n";
  cout << "BTF Output\n";
  cout << "-----------------------------------------\n";
  cout << "RowPerm and ColPerm\n";
  for( int i = 0; i<n; ++i )
    cout << rowperm[i] << "\t" << colperm[i] << endl;
  if( hrzcmp )
  {
    cout << "Num Horizontal: Rows, Cols, Comps\n";
    cout << nhrows << "\t" << nhcols << "\t" << hrzcmp << endl;
  }
  cout << "Num Square: Rows, Comps\n";
  cout << nsrows << "\t" << sqcmpn << endl;
  if( vrtcmp )
  {
    cout << "Num Vertical: Rows, Cols, Comps\n";
    cout << nvrows << "\t" << nvcols << "\t" << vrtcmp << endl;
  }
  cout << "Row, Col of upper left pt in blocks\n";
  for( int i = 0; (i<n+1) && (rcmstr[i]!=n+1); ++i )
    cout << i << " " << rcmstr[i] << " " << ccmstr[i] << endl;
  cout << "-----------------------------------------\n";
#endif

  if( hrzcmp || vrtcmp )
  { cout << "FAILED! hrz cmp's:" << hrzcmp << " vrtcmp: " << vrtcmp << endl;
    exit(0); }

  //convert rowperm to OLD->NEW
  //reverse ordering of permutation to get upper triangular
  vector<int> rowperm_t( n );
  vector<int> colperm_t( n );
  for( int i = 0; i < n; ++i )
  {
//    rowperm_t[ rowperm[i] ] = n-i;
//    colperm[i] = n-colperm[i];
    rowperm_t[i] = rowperm[(n-1)-i];
    colperm_t[i] = colperm[(n-1)-i];
  }

  //Generate New Domain and Range Maps
  //for now, assume they start out as identical
  const Epetra_BlockMap & OldMap = original.RowMap();
  vector<int> myElements( n );
  OldMap.MyGlobalElements( &myElements[0] );

  vector<int> newDomainElements( n );
  vector<int> newRangeElements( n );
  for( int i = 0; i < n; ++i )
  {
    newDomainElements[ i ] = myElements[ rowperm_t[i] ];
    newRangeElements[ i ] = myElements[ colperm_t[i] ];
cout << i << "\t" << rowperm_t[i] << "\t" << colperm[i] << "\t" << myElements[i] << endl;
  }

  Epetra_Map * DomainMap = new Epetra_Map( n, n, &newDomainElements[0], OldMap.IndexBase(), OldMap.Comm() );
  Epetra_Map * RangeMap = new Epetra_Map( n, n, &newRangeElements[0], OldMap.IndexBase(), OldMap.Comm() );

#ifdef BTF_VERBOSE
  cout << "New Domain Map\n";
  cout << *DomainMap << endl;
  cout << "New Range Map\n";
  cout << *RangeMap << endl;
#endif

  //Generate New Graph
  std::auto_ptr<Epetra_CrsGraph> NewGraph( new Epetra_CrsGraph( Copy, *DomainMap, 0 ) );
  Epetra_Import Importer( *DomainMap, OldMap );
  NewGraph->Import( original, Importer, Insert );
  NewGraph->TransformToLocal( DomainMap, RangeMap );

#ifdef BTF_VERBOSE
  cout << "New CrsGraph\n";
  cout << *NewGraph << endl;
#endif

  return NewGraph;
}

} //namespace Epetra_Transform
