#include "SparseSolverResult.h"
#include <iomanip>

void SparseSolverResult::Print(ostream & os) const {
  if ( error != UnUsedDbl ) { os << "error=" << error << " " ; }
  if ( residual != UnUsedDbl ) { os << "residual=" << residual << endl ; }
  if ( Anorm != UnUsedDbl ) { os << "Anorm=" << Anorm << " " ; }
  if ( Bnorm != UnUsedDbl ) { os << "Bnorm=" << Bnorm << endl ; }
  if ( Xnorm != UnUsedDbl ) { os << "Xnorm=" << Xnorm << endl ; }

  { os << "Data redistributions times: " << setw(9) << setprecision(3) ; } 
  { os << "Wall( " ;
  if (RedistribTime_.WallTime() == UnUsedDbl ) 
    os << "Unknown"; else os  << RedistribTime_.WallTime();
  os << " ) " ; }
  { os << "System( " ;
  if (RedistribTime_.SystemTime() == UnUsedDbl ) 
    os << "Unknown"; else os  << RedistribTime_.SystemTime();
  os << " ) " ; }
  { os << "User( " ;
  if (RedistribTime_.UserTime() == UnUsedDbl ) 
    os << "Unknown"; else os  << RedistribTime_.UserTime();
  os << " )" << endl ; }

  { os << "Symbolic factorization times:  Wall( " ;
  if (SymbolicTime_.WallTime() == UnUsedDbl ) 
    os << "Unknown"; else os  << SymbolicTime_.WallTime();
  os << " ) " ; }
  { os << "System( " ;
  if (SymbolicTime_.SystemTime() == UnUsedDbl ) 
    os << "Unknown"; else os  << SymbolicTime_.SystemTime();
  os << " ) " ; }
  { os << "User( " ;
  if (SymbolicTime_.UserTime() == UnUsedDbl ) 
    os << "Unknown"; else os  << SymbolicTime_.UserTime();
  os << " )" << endl ; }

  { os << "Numeric factorization times:   Wall( " ;
  if (FactorTime_.WallTime() == UnUsedDbl ) 
    os << "Unknown"; else os  << FactorTime_.WallTime();
  os << " ) " ; }
  { os << "System( " ;
  if (FactorTime_.SystemTime() == UnUsedDbl ) 
    os << "Unknown"; else os  << FactorTime_.SystemTime();
  os << " ) " ; }
  { os << "User( " ;
  if (FactorTime_.UserTime() == UnUsedDbl ) 
    os << "Unknown"; else os  << FactorTime_.UserTime();
  os << " )" << endl ; }

  { os << "Backsolve times:               Wall( " ;
  if (SolveTime_.WallTime() == UnUsedDbl ) 
    os << "Unknown"; else os  << SolveTime_.WallTime();
  os << " ) " ; }
  { os << "System( " ;
  if (SolveTime_.SystemTime() == UnUsedDbl ) 
    os << "Unknown"; else os  << SolveTime_.SystemTime();
  os << " ) " ; }
  { os << "User( " ;
  if (SolveTime_.UserTime() == UnUsedDbl ) 
    os << "Unknown"; else os  << SolveTime_.UserTime();
  os << " )" << endl ; }


}

void SparseSolverResult::PrintSummary(ostream & os) const {
  os << setw(9) << setprecision(2) << Anorm ; 
  if ( error != UnUsedDbl && Xnorm != UnUsedDbl ) 
    { os << setw(9) << setprecision(2) << error/Xnorm << " " ; } else os << " error is unknown " ; 
  if ( residual != UnUsedDbl && Bnorm != UnUsedDbl ) 
    { os << setw(9) << setprecision(2) << residual/Bnorm << " " ; } else os << " residual is unknown " ; 
  double TotalWallClock = 0.0 ; 
  if ( RedistribTime_.WallTime() != UnUsedDbl )
    { TotalWallClock += RedistribTime_.WallTime() ; } ; 
  if ( SymbolicTime_.WallTime() != UnUsedDbl )
    { TotalWallClock += SymbolicTime_.WallTime() ; } ; 
  if ( FactorTime_.WallTime() != UnUsedDbl )
    { TotalWallClock += FactorTime_.WallTime() ; } ; 
  if ( SolveTime_.WallTime() != UnUsedDbl )
    { TotalWallClock += SolveTime_.WallTime() ; } ; 
  os << setw(10) << setprecision(4) << total_time  ; 
  os << setw(10) << setprecision(4) << TotalWallClock  ; 
  if ( true || first_time != 0.0 ) {
    os << setw(10) << setprecision(4) << first_time  ; 
    os << setw(10) << setprecision(4) << last_time - middle_time ; 
    os << setw(10) << setprecision(4) << last_time - first_time ; 
    
  }

}
