// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
// @HEADER

/*! \file Galeri_core_Workspace
 *
 * \author Marzio Sala, ETHZ
 *
 * \date Last modified on Aug-06
 */

#ifndef GALERI_CORE_WORKSPACE_H
#define GALERI_CORE_WORKSPACE_H

#include "Teuchos_TestForException.hpp"

class Epetra_Comm;
class Epetra_RowMatrix;
class Epetra_MultiVector;

#define GALERI_MAX(x,y) (( (x) > (y) ) ? x : y)
#define GALERI_MIN(x,y) (( (x) < (y) ) ? x : y) 
#define GALERI_SGN(x) (((x) < 0.0) ? -1.0 : 1.0) 

namespace Galeri {
namespace core {

/*!
 * \class Workspace
 *
 * \brief Function class containing a few static methods and constants, to be
 * used as workspace tools.
 */
class Workspace 
{
  public:
    //! Default constructor.
    Workspace() {}

    //! Default destructor.
    ~Workspace() {}

    //! Sets the number of dimension (1, 2, or 3).
    static void setNumDimensions(const int numDimensions)
    {
      numDimensions_ = numDimensions;
    }

    //! Returns the number of dimensions used (1, 2, or 3).
    static int getNumDimensions()
    {
      return(numDimensions_);
    }

    //! Solves a serial linear system using LAPACK (WARNING: ONLY SMALL MATRICES)
    static void solve_LAPACK(Epetra_RowMatrix& matrix,
                             Epetra_MultiVector& LHS,
                             Epetra_MultiVector& RHS);

    /*! \brief Creates a multivector that can hold a component of the specified 
     * multivector.
     *
     * From the input Epetra_MultiVector, defined on an Epetra_BlockMap,
     * defines a new Epetra_Map that can host one component of the
     * multivector, and allocates an Epetra_MultiVector based on it. This is
     * useful in vector problems, when the input vector contains several PDE
     * componets, and one of them has to be extracted.
     */
    static 
    Epetra_MultiVector* createMultiVectorComponent(const Epetra_MultiVector& input);

    /*! \brief Extracts the component of the specified \c equation from the \c
     * input Epetra_MultiVector, and stores it in \c output.
     */
    static 
    void extractMultiVectorComponent(const Epetra_MultiVector& input,
                                     const int equation,
                                     Epetra_MultiVector& output);

    //! Input default value for "min".
    static const int MIN;
    //! Input default value for "max".
    static const int MAX;
    //! Default value for uninitialized objects.
    static const int UNINITIALIZED;
    //! Default value for initialized objects (should be > UNINITIALIZED).
    static const int INITIALIZED;
    //! Value for status with freezed connectivity (should be > INITIALIZED).
    static const int CONNECTIVITY_FREEZED;
    //! Value for status with freezed coordinates (should be > INITIALIZED).
    static const int COORDINATES_FREEZED;

  private:
    //! Number of dimensions in the computations.
    static int numDimensions_;
}; // class Workspace

} // namespace core
} // namespace Galeri

#endif
