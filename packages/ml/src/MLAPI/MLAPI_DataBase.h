#ifndef MLAPI_DATABASE_H
#define MLAPI_DATABASE_H

#include "MLAPI_BaseObject.h"
#include "Teuchos_ParameterList.hpp"

namespace MLAPI {

const int ALL_LEVELS = -1;

class SmootherDataBase : public BaseObject {

public:

  // @{ Constructors and destructors

  //! Initializes to default values.
  SmootherDataBase(const int MaxLevels = 20) 
  {
    MaxLevels_ = 20;
    Type_      = "symmetric Gauss-Seidel";
    Damping_   = 0.67;
    Sweeps_    = 1;
    ILUFill_   = 0;
    ILUTFill_  = 1.0;
    ICTFill_   = 1.0;
    CurrentLevel_ = 0;
  }
  
  SmootherDataBase(Teuchos::ParameterList& List)
  {
    MaxLevels_ = List.get("max levels", 10);
    Type_      = List.get("smoother: type", "symmetric Gauss-Seidel");
    Damping_   = List.get("smoother: damping factor", 0.67);
    Sweeps_    = List.get("smoother: sweeps", 1);
    ILUFill_   = 0;
    ILUTFill_  = 1.0;
    ICTFill_   = 1.0;
    CurrentLevel_ = 0;
  }

  //! Dtor
  ~SmootherDataBase()
  {}

  // @}
  // @{ Set and get methods

  int GetMaxLevels() const
  {
    return(MaxLevels_);
  }

  void SetCurrentLevel(const int CurrentLevel)
  {
    CurrentLevel_ = CurrentLevel;
  }

  int GetCurrentLevel() const
  {
    return(CurrentLevel_);
  }

  string GetType() const
  {
    return(Type_);
  }
  
  void SetType(const string Type, const int level = ALL_LEVELS)
  {
    Type_ = Type;
  }

  void SetSweeps(const int Sweeps, const int level = ALL_LEVELS)
  {
    Sweeps_ = Sweeps;
  }

  int GetSweeps() const
  {
    return(Sweeps_);
  }

  void SetDamping(const int Damping, const int level = ALL_LEVELS)
  {
    Damping_ = Damping;
  }

  double GetDamping() const
  {
    return(Damping_);
  }

  void SetILUFill(const int ILUFill, const int level = ALL_LEVELS)
  {
    ILUFill_ = ILUFill;
  }

  int GetILUFill() const
  {
    return(ILUFill_);
  }

  void SetILUTFill(const double ILUTFill, const int level = ALL_LEVELS)
  {
    ILUTFill_ = ILUTFill;
  }

  double GetILUTFill() const
  {
    return(ILUTFill_);
  }

  void SetICTFill(const double ICTFill, const int level = ALL_LEVELS)
  {
    ICTFill_ = ICTFill;
  }

  double GetICTFill() const
  {
    return(ICTFill_);
  }

  //! Prints basic information about \c this object.
  ostream& Print(ostream& os, bool verbose) const 
  {
    os << std::endl;
    if (MyPID() == 0) {
      os << "*** MLAPI::SmootherDataBase ***" << endl;
      os << "Label                 = " << GetLabel() << endl;
      os << endl;
    }
    return(os);
  }

private:  

  //! Maximum number of levels
  int    MaxLevels_;
  int    CurrentLevel_;
  string Type_;
  int    Sweeps_;
  double Damping_;
  int    ILUFill_;
  double ILUTFill_;
  double ICTFill_;

};

class CoarseSolverDataBase : public BaseObject {

public:

  // @{ Constructors and destructors

  //! Initializes to default values.
  CoarseSolverDataBase(const int MaxLevels = 20) 
  {
    Type_ = "Uncoupled";
  }
  
  CoarseSolverDataBase(Teuchos::ParameterList& List)
  {
    Type_     = List.get("coarse: type", "Amesos-KLU");
  }

  //! Dtor
  ~CoarseSolverDataBase()
  {}

  // @}
  // @{ Set and get methods

  //! Sets the coarsening scheme.
  void SetType(const string Type, const int level = ALL_LEVELS)
  {
     Type_ = Type;
  }

  //! Gets the coarsening string.
  string GetType() const
  {
    return(Type_);
  }

  //! Prints basic information about \c this object.
  ostream& Print(ostream& os, bool verbose) const 
  {
    os << std::endl;
    if (MyPID() == 0) {
      os << "*** MLAPI::CoarseSolverDataBase ***" << endl;
      os << endl;
    }
    return(os);
  }

private:  

  //! Coarsening scheme.
  string Type_;

};

class AggregationDataBase : public BaseObject {

public:

  // @{ Constructors and destructors

  //! Initializes to default values.
  AggregationDataBase(const int MaxLevels = 20) 
  {
    MaxLevels_     = MaxLevels_;
    Type_          = "Uncoupled";
    NPA_           = 64;
    EigenAnalysis_ = "Anorm";
    Damping_       = 1.333;
    CurrentLevel_  = 0;
  }

  AggregationDataBase(Teuchos::ParameterList& List)
  {
    MaxLevels_     = List.get("max levels", 10);
    Type_          = List.get("aggregation: type", "Uncoupled");
    Damping_       = List.get("aggregation: damping factor", 1.333);
    Threshold_     = List.get("aggregation: threshold", 0.0);
    NPA_           = List.get("aggregation: nodes per aggregate", 64);
    EigenAnalysis_ = List.get("eigen-analysis: type", "Anorm");
    CurrentLevel_  = 0;
  }

  //! Dtor
  ~AggregationDataBase()
  {}

  // @}
  // @{ Set and get methods

  //! Sets the coarsening scheme.
  void SetType(const string Type, const int level = ALL_LEVELS)
  {
    Type_ = Type;
  }

  //! Gets the coarsening string.
  string GetType() const
  {
    return(Type_);
  }

  //! Sets the threshold to be used during the aggregate construction.
  void SetThreshold(const double thresh, const int level = ALL_LEVELS)
  {
     Threshold_ = thresh;
  }

  //! Gets the threshold to be used during the aggregate construction.
  double GetThreshold() const
  {
    return(Threshold_);
  }

  //! Sets the ideal number of nodes on each aggregate (for METIS and ParMETIS only).
  void SetNodesPerAggregate(const int NPA, const int level = ALL_LEVELS)
  {
     NPA_ = NPA;
  }

  //! Gets the ideal number of nodes on each aggregate (for METIS and ParMETIS only).
  int GetNodesPerAggregate() const
  {
    return(NPA_);
  }

  int GetMaxLevels() const
  {
    return(MaxLevels_);
  }

  double GetDamping() const
  {
    return(Damping_);
  }
  
  void SetDamping(const double Damping, const int level = ALL_LEVELS)
  {
    Damping_ = Damping;
  }
  
  string GetEigenAnalysisType() const
  {
    return(EigenAnalysis_);
  }
  
  void SetEigenAnalysisType(const string EigenAnalysis, const int level = ALL_LEVELS)
  {
    EigenAnalysis_ = EigenAnalysis;
  }
  
  void SetCurrentLevel(const int CurrentLevel)
  {
    CurrentLevel_ = CurrentLevel;
  }

  int GetCurrentLevel() const
  {
    return(CurrentLevel_);
  }

  //! Prints basic information about \c this object.
  ostream& Print(ostream& os, bool verbose) const 
  {
    os << std::endl;
    if (MyPID() == 0) {
      os << "*** MLAPI::AggregationDataBase ***" << endl;
      os << "Label                 = " << GetLabel() << endl;
      os << endl;
    }
    return(os);
  }

private:  

  //! Maximum number of levels
  int MaxLevels_;
  //! Threshold for aggregate construction.
  double Threshold_;
  //! Number of nodes per aggregate (METIS and ParMETIS only)
  int    NPA_;
  string Type_;
  //! Damping factor
  double Damping_;
  //! Eigen-analysis type
  string EigenAnalysis_;
  int CurrentLevel_;

};
} // namespace MLAPI

#endif // MLAPI_DATABASE_H
