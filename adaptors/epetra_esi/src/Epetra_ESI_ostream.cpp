#ifndef _Epetra_ESI_ostream_cpp_
#define _Epetra_ESI_ostream_cpp_

//
//A convenience: an ostream << operator for epetra_esi::CrsMatrix.
//
ostream& operator<<(ostream& os, epetra_esi::CrsMatrix<double,int>& mat)
{
  //This function allows an instance of a "epetra_esi::CrsMatrix" to
  //be given directly to an ostream for output.
  //
  esi::IndexSpace<int>* indexspace = NULL, *dummy = NULL;
  int err = mat.getIndexSpaces(indexspace, dummy);
  if (err != 0) {
    cerr << "Epetra_ESI ostrm ERROR in mat.getIndexSpaces"<<endl;
    return(os);
  }

  int localSize = 0, localOffset = 0;
  err = indexspace->getLocalSize(localSize);
  err += indexspace->getLocalPartitionOffset(localOffset);

  if (err != 0) {
    cerr <<"Epetra_ESI ostrm ERROR in indexspace->getLocalSize or offset."<<endl;
    return(os);
  }


  Epetra_Array<double> coefs(50);
  Epetra_Array<int> colInds(50);

  for(int i=0; i<localSize; i++) {
    int numnz = 0;
    err = mat.getRowNonzeros(i+localOffset, numnz);
    if (err != 0) {
      cerr << "Epetra_ESI ostrm ERROR in mat.getRowNonzeros" << endl;
      return(os);
    }

    coefs.resize(numnz);
    colInds.resize(numnz);
    err = mat.copyOutRow(i+localOffset, coefs.dataPtr(), colInds.dataPtr(),
                         colInds.length(), numnz);
    if (err != 0) {
      cerr << "Epetra_ESI ostrm ERROR in mat.copyOutRow " << i+localOffset<<endl;
      return(os);
    }

    for(int j=0; j<numnz; j++) {
      os << " " << i+localOffset << " " << colInds[j] << " " << coefs[j]<<endl;
    }
  }

  return(os);
};

#endif

