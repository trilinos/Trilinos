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

#include "Isorropia_EpetraMatcher.hpp"
#include "Isorropia_EpetraRedistributor.hpp"

#ifdef HAVE_EPETRAEXT
#include "EpetraExt_Reindex_CrsMatrix.h"
#include "EpetraExt_CrsMatrixIn.h"
#endif


using namespace std;

int main(int argc, char** argv) {

#ifdef ISORROPIA_HAVE_OMP
        if(argc>2)
        {
                int rc=0;
#ifdef HAVE_EPETRAEXT
                  const Epetra_SerialComm Comm;

                  Epetra_CrsMatrix *matrixPtr;
                  rc = EpetraExt::MatrixMarketFileToCrsMatrix(argv[1], Comm, matrixPtr);

                  if (rc < 0){
            cout << "error reading input file" << std::endl;
            return 1;
                  }

                Teuchos::ParameterList paramlist;
                paramlist.set("Matching Algorithm",argv[2]);
        Isorropia::Epetra::Matcher pm(matrixPtr,paramlist);

        //Teuchos::RCP<const Epetra_CrsMatrix> r(Teuchos::RCP<const
        //Epetra_CrsMatrix>(matrixPtr,true));
        //Isorropia::Epetra::Matcher pm(r,paramlist);

        pm.match();

        cout << endl << "Original Matrix:" << endl;
        std::cout<<*matrixPtr<<std::endl;

        Teuchos::RCP<Epetra_CrsMatrix> perm_matrix =
                                        pm.applyColumnPermutation();
        cout << endl << "After Column permutation:" << endl;
        cout << *perm_matrix << endl;

        perm_matrix = pm.applyRowPermutation();
        cout << endl << "After Row permutation:" << endl;
        cout << *perm_matrix << endl;

#else
                 fail = 0;
         cout << "Matching test requires EpetraExt" << std::endl;
         return 1;
#endif
        }
        else
    {
                cout<<endl<<" Usage: ./Isorropia_parallel_matching.exe <mtx file>" <<
             " <Algorithm>" << endl;
        cout << "\t Algorithm: PHK, PHKDW, PDFS,PPF" << endl;
        cout << "Requires Isorropia to be compiled with OpenMP, OMP_NUM_THREADS"
              << "set to at least one and the test requires EpetraExt." <<
               endl << endl;
    }
#else
    cout << "Matching in Isorropia requires OpenMP." << endl;
    cout << "Please recompile with OpenMP enabled and try again" << endl;
#endif

        return 0;
}
