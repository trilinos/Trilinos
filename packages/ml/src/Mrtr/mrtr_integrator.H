/*
#@HEADER
# ************************************************************************
#
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifdef TRILINOS_PACKAGE

#ifndef MRTR_INTEGRATOR_H
#define MRTR_INTEGRATOR_H

#include <ctime>
#include <iostream>
#include <iomanip>
#include <vector>

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_RefCountPtr.hpp"

#include "mrtr_overlap.H"

using namespace std;
using namespace Teuchos;

// ----------   User Defined Includes   ----------

namespace MRTR
{

// forward declarations
class MRTR::Interface;
class MRTR::Segment;
class MRTR::Node;

class Integrator 
{
public:
  
  // ctor
  explicit Integrator(int ngp, bool oneD);
  
  // dtor
  virtual ~Integrator();
  
  // get the order of integration (== # gaussian points)
  inline int Ngp() { return ngp_; }
  
  // return the coordinate of a given gaussian point gp
  inline double Coordinate(int gp) { return coords_[gp]; }
  inline double* Coordinate(int* gp) { return &coords_[*gp*2]; }
  
  // return the weight of a given gaussian point gp
  inline double Weight(int gp) { return weights_[gp]; }


  // integrate a 1D overlap between a slave and a master segment
  Epetra_SerialDenseMatrix* Integrate(MRTR::Segment& sseg, double sxia, double sxib,
                                      MRTR::Segment& mseg, double mxia, double mxib);
  // assemble the result -Mdense from the integration above into M
  bool Assemble(MRTR::Interface& inter, MRTR::Segment& sseg, MRTR::Segment& mseg, 
                Epetra_CrsMatrix& M, Epetra_SerialDenseMatrix& Mdense);


  // integrate a 1D slave segment
  Epetra_SerialDenseMatrix* Integrate(MRTR::Segment& sseg, double sxia, double sxib);
  // assemble the result Ddense from the integration above into D
  bool Assemble(MRTR::Interface& inter, MRTR::Segment& sseg,  
                Epetra_CrsMatrix& D, Epetra_SerialDenseMatrix& Ddense);


  // integrate modification of M in 1D
  Epetra_SerialDenseMatrix* Integrate_2D_Mmod(MRTR::Segment& sseg, double sxia, double sxib,
                                              MRTR::Segment& mseg, double mxia, double mxib);
  // assemble the result -Mdense from the integration above into M
  bool Assemble_2D_Mod(MRTR::Interface& inter, MRTR::Segment& sseg, MRTR::Segment& mseg, 
                Epetra_CrsMatrix& M, Epetra_SerialDenseMatrix& Mmod);


  // integrate a 2D triangle overlap segment, master part M and slave part D
  bool Integrate(RefCountPtr<MRTR::Segment> actseg,
                 MRTR::Segment& sseg, 
                 MRTR::Segment& mseg,
                 Epetra_SerialDenseMatrix** Ddense, 
                 Epetra_SerialDenseMatrix** Mdense, 
                 MRTR::Overlap& overlap);

  // Assemble D into nodes
  bool Assemble(MRTR::Interface& inter,MRTR::Segment& sseg,Epetra_SerialDenseMatrix& Ddense);

  // Assemble M into nodes
  bool Assemble(MRTR::Interface& inter,MRTR::Segment& sseg,MRTR::Segment& mseg,
                Epetra_SerialDenseMatrix& Mdense);

private:  
  // don't want = operator
  Integrator operator = (const Integrator& old);
  // don't want copy-ctor
  Integrator(MRTR::Integrator& old);

private:

  bool           oneD_;    // dimension of the integration, true if 1-dimensional
  int            ngp_;     // order of the integration
  vector<double> coords_;  // vector holding gaussian point coordinates
  vector<double> weights_; // vector holding gaussian point weights


};

} // namespace MRTR

#endif // MRTR_INTEGRATOR_H
#endif // TRILINOS_PACKAGE
