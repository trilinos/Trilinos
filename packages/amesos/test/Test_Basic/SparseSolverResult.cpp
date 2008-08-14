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
// @HEADER

#include "SparseSolverResult.h"
#include <iomanip>

void SparseSolverResult::Print(std::ostream & os) const {
  if ( error != UnUsedDbl ) { os << "error=" << error << " " ; }
  if ( residual != UnUsedDbl ) { os << "residual=" << residual << std::endl ; }
  if ( Anorm != UnUsedDbl ) { os << "Anorm=" << Anorm << " " ; }
  if ( Bnorm != UnUsedDbl ) { os << "Bnorm=" << Bnorm << std::endl ; }
  if ( Xnorm != UnUsedDbl ) { os << "Xnorm=" << Xnorm << std::endl ; }

}

void SparseSolverResult::PrintSummary(std::ostream & os) const {
  os << std::setw(9) << std::setprecision(2) << Anorm ; 
  if ( error != UnUsedDbl && Xnorm != UnUsedDbl ) 
    { os << std::setw(9) << std::setprecision(2) << error/Xnorm << " " ; } else os << " error is unknown " ; 
  if ( residual != UnUsedDbl && Bnorm != UnUsedDbl ) 
    { os << std::setw(9) << std::setprecision(2) << residual/Bnorm << " " ; } else os << " residual is unknown " ; 
  os << std::setw(10) << std::setprecision(4) << total_time  ; 
#if 0
  double TotalWallClock = 0.0 ; 
  if ( RedistribTime_.WallTime() != UnUsedDbl )
    { TotalWallClock += RedistribTime_.WallTime() ; } ; 
  if ( SymbolicTime_.WallTime() != UnUsedDbl )
    { TotalWallClock += SymbolicTime_.WallTime() ; } ; 
  if ( FactorTime_.WallTime() != UnUsedDbl )
    { TotalWallClock += FactorTime_.WallTime() ; } ; 
  if ( SolveTime_.WallTime() != UnUsedDbl )
    { TotalWallClock += SolveTime_.WallTime() ; } ; 
  os << std::setw(10) << std::setprecision(4) << TotalWallClock  ; 
#endif
  if ( first_time != UnUsedDbl ) 
    os << std::setw(10) << std::setprecision(4) << first_time  ; 
  else
    os << "        na " ; 
  if ( middle_time != UnUsedDbl ) 
    os << std::setw(10) << std::setprecision(4) << last_time - middle_time ; 
  else
    os << "        na " ; 
  if ( last_time != UnUsedDbl ) 
    os << std::setw(10) << std::setprecision(4) << last_time - first_time ; 
  else
    os << "        na " ; 

}
