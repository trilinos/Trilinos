#include "SparseSolverResult.h"
#include <iomanip>

void SparseSolverResult::Print(ostream & os) const {
  if ( error != UnUsedDbl ) { os << "error=" << error << " " ; }
  if ( residual != UnUsedDbl ) { os << "residual=" << residual << endl ; }
  if ( Anorm != UnUsedDbl ) { os << "Anorm=" << Anorm << " " ; }
  if ( Bnorm != UnUsedDbl ) { os << "Bnorm=" << Bnorm << endl ; }
  if ( Xnorm != UnUsedDbl ) { os << "Xnorm=" << Xnorm << endl ; }

}

void SparseSolverResult::PrintSummary(ostream & os) const {
  os << setw(9) << setprecision(2) << Anorm ; 
  if ( error != UnUsedDbl && Xnorm != UnUsedDbl ) 
    { os << setw(9) << setprecision(2) << error/Xnorm << " " ; } else os << " error is unknown " ; 
  if ( residual != UnUsedDbl && Bnorm != UnUsedDbl ) 
    { os << setw(9) << setprecision(2) << residual/Bnorm << " " ; } else os << " residual is unknown " ; 
  os << setw(10) << setprecision(4) << total_time  ; 
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
  os << setw(10) << setprecision(4) << TotalWallClock  ; 
#endif
  os << setw(10) << setprecision(4) << first_time  ; 
  os << setw(10) << setprecision(4) << last_time - middle_time ; 
  os << setw(10) << setprecision(4) << last_time - first_time ; 

}
