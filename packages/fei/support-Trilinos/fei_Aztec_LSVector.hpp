/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/

#ifndef _fei_Aztec_LSVector_hpp_
#define _fei_Aztec_LSVector_hpp_


#include <fei_SharedPtr.hpp>

//
// This class provides a vector that can be used with the AztecDMSR_matrix
// and the AztecDVBR_Matrix.
//
// An important restriction to note:
//
// * An Aztec_LSVector can not be constructed until AFTER the AztecDMSR_Matrix
//   (or AztecDVBR_Matrix) that it is to be used with has been completely
//   initialized and filled (e.g., A.loadComplete() has been called (which
//   means, most importantly, that AZ_transform has been called)). This is
//   because the local data array for an aztec vector must be allocated
//   with enough extra space to hold 'boundary elements' that are exchanged
//   with other processors during the calculation of a parallel matrix-vector
//   product, and we don't know how much memory that requires until after
//   AZ_transform has been called.
//
// * Also, the calling code is responsible for keeping track of any 
//   re-ordering that AZ_transform has done. i.e., Aztec_LSVector is just
//   like a raw array with respect to indexing of entries. If v is an
//   instantiation of an Aztec_LSVector, then v[9] literally returns the
//   entry at position 9 (the 10th entry, since indexing is 0-based).
//
namespace fei_trilinos {

class Aztec_Map;

/**==========================================================================**/
class Aztec_LSVector {
  public:
    // Constructor.
    Aztec_LSVector(fei::SharedPtr<Aztec_Map> map, int* data_org);

    Aztec_LSVector(const Aztec_LSVector& source);  // copy constructor

    virtual ~Aztec_LSVector ();

    Aztec_LSVector* newVector() const;

    // Mathematical functions.
    double dotProd (const Aztec_LSVector& y) const;
    void scale (double s);
    void addVec (double s, const Aztec_LSVector& c);
    double norm () const;
    double norm1 () const;
 
    // operator=
    Aztec_LSVector& operator = (const Aztec_LSVector& rhs);
    
    // Access functions.
    double& operator [] (int index);
    const double& operator [] (int index) const;
    
    void put (double scalar);

    const double* startPointer() const {return localCoeffs_;};

    //Special function
    bool readFromFile(const char *fileName);
    bool writeToFile(const char *fileName) const;
    
  protected:
    virtual void assign(const Aztec_LSVector& rhs);
    
  private:
    void checkInput();
    int inUpdate(int globalIndex, int& localIndex) const;

    fei::SharedPtr<Aztec_Map> amap_;
    double *localCoeffs_;        // local vector coefficients
    int length_;
};

}//namespace fei_trilinos

#endif
