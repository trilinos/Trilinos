// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

/*!
 * \file Amesos_MC64.h
 *
 * \class Amesos_MC64
 *
 * \brief Interface to MC64, reordering and scaling algorithm.
 *
 * \author Marzio Sala, ETHZ.
 *
 * \date Last updated on Feb-06.
 */

#ifndef AMESOS_MC64_H
#define AMESOS_MC64_H

#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif

#include "Amesos_ConfigDefs.h"
#if defined(HAVE_AMESOS_MC64)
#include "Amesos_Scaling.h"
#include "Amesos_Reordering.h"

class Epetra_RowMatrix;

//! Amesos_MC64: Interface to the reordering and scaling algorithm HSL's MC64.
/*! Amesos_MC64 is an interface to HSL's MC64. Amesos must be configured with
 * option --enable-amesos-mc64 to take advantage of this class. The location
 * of the MC64 algorithm is specified by using --with-ldflags and
 * --with-libs.
 *
 *  \author Marzio Sala, ETHZ.
 *
 *  \date Last updated on Feb-06.
 */
class Amesos_MC64 : public Amesos_Scaling, public Amesos_Reordering
{
  public:
    //! Constructor, requires matrix to reorder and scale.
    Amesos_MC64(const Epetra_RowMatrix& A, int JOB, 
                const bool StoreTranspose = false,
                const bool analyze = false);

    //! Destructor.
    ~Amesos_MC64() {}

    //! Returns the specified position (in FORTRAN style) of the INFO array.
    int GetINFO(const int pos) const
    {
      if (pos <= 0 || pos > 10)
        throw(-1);

      return(INFO_[pos - 1]);
    }

    //! Returns the row permutation vector, or 0 if not computed.
    int* GetRowPerm() 
    {
      return(0);
    }

    //! Returns the column permutation vector, or 0 if not computed.
    int* GetColPerm()
    {
      return((int*)&CPERM_[0]);
    }

    //! Returns the row scaling vector, or 0 if not computed.
    double* GetRowScaling()
    {
      return((double*)&DW_[0]);
    }

    //! Returns the column scaling vector, or 0 if not computed.
    double* GetColScaling();

    //! Returns a pointer to the internally stored CPERM int array.
    int* GetCPERM()
    {
      return((int*)&CPERM_[0]);
    }

    //! Returns a pointer to the internally stored DW double array.
    double* GetDW()
    {
      return((double*)&DW_[0]);
    }

  private:
    //
    //! Computes the column reordering CPERM and scaling DW.
    int Compute(int JOB, const bool StoreTranspose, const bool analyze);

    const Epetra_RowMatrix& A_;
    int ICNTL_[10];
    int INFO_[10];
    std::vector<int> CPERM_;
    std::vector<double> DW_;
};

#endif
#endif
