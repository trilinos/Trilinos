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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//
/*!
 * \file Amesos_Reordering.h
 *
 * \class Amesos_Reordering
 *
 * \brief Base class for reordering procedures.
 *
 * \author Marzio Sala, ETHZ.
 *
 * \date Last updated on Feb-06.
 */

#ifndef AMESOS_REORDERING_H
#define AMESOS_REORDERING_H

//! Amesos_Reordering: base class for reordering procedures.
class Amesos_Reordering
{
  public:

    ~Amesos_Reordering() {}

    //! Returns the row permutation vector, or 0 if not computed.
    virtual int* GetRowPerm() = 0;

    //! Returns the column permutation vector, or 0 if not computed.
    virtual int* GetColPerm() = 0;
};

#endif
