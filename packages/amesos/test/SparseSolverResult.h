#ifndef SPARSESOLVERRESULT
#define SPARSESOLVERRESULT
#include "Epetra_Object.h"
const double UnUsedDbl = 13e130;

class SparseSolverResult : public Epetra_Object { 
  
 public:
  SparseSolverResult() :  
    first_time( UnUsedDbl ), middle_time( UnUsedDbl ), 
    last_time( UnUsedDbl ), total_time( UnUsedDbl ),
    error( UnUsedDbl), residual(UnUsedDbl), 
    Anorm( UnUsedDbl ), Xnorm( UnUsedDbl ), Bnorm( UnUsedDbl ) 
  { 
  ; } ; 
  ~SparseSolverResult(){};

  void Set_First_Time ( double First_Time_in ) { first_time = First_Time_in; } ; 
  inline double Get_First_Time ( ) { return first_time ; } ; 

  void Set_Middle_Time ( double Middle_Time_in ) { middle_time = Middle_Time_in; } ; 
  inline double Get_Middle_Time ( ) { return middle_time ; } ; 

  void Set_Last_Time ( double Last_Time_in ) { last_time = Last_Time_in; } ; 
  inline double Get_Last_Time ( ) { return last_time ; } ; 

  void Set_Total_Time ( double Total_Time_in ) { total_time = Total_Time_in; } ; 
  inline double Get_Total_Time ( ) { return total_time ; } ; 

  void Set_Bnorm ( double bnorm_in ) { Bnorm = bnorm_in; } ; 
  inline double Get_Bnorm ( ) { return Bnorm ; } ; 
  void Set_Xnorm ( double xnorm_in ) { Xnorm = xnorm_in; } ; 
  inline double Get_Xnorm ( ) { return Xnorm ; } ; 
  void Set_Residual ( double residual_in ) { residual = residual_in; } ; 
  inline double Get_Residual ( ) { return residual ; } ; 
  void Set_Error ( double error_in ) { error = error_in; } ; 
  inline double Get_Error ( ) { return error; } ; 
  void Set_Anorm ( double anorm_in ) { Anorm = anorm_in; } ; 
  inline double Get_Anorm ( ) { return Anorm; } ; 

  virtual void Print(ostream & os) const;
  virtual void PrintSummary(ostream & os) const;


 private:
  double first_time ;
  double middle_time ;
  double last_time ;
  double total_time ;
  double error ;
  double residual ;
  double Xnorm ;
  double Bnorm ;
  double Anorm ;

} ;

#endif
