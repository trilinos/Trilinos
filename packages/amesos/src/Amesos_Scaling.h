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
//
/*!
 * \file Amesos_Scaling.h
 *
 * \class Amesos_Scaling
 *
 * \brief Base class for scaling procedures.
 *
 * \author Marzio Sala, ETHZ.
 *
 * \date Last updated on Feb-06.
 */

#ifndef AMESOS_SCALING_H
#define AMESOS_SCALING_H

#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif

//! Amesos_Scaling: base class for scaling procedures.
class Amesos_Scaling
{
  public:

    ~Amesos_Scaling() {}

    //! Returns the row scaling vector, or 0 if not computed.
    virtual double* GetRowScaling() = 0;

    //! Returns the column scaling vector, or 0 if not computed.
    virtual double* GetColScaling() = 0;
};

#endif
