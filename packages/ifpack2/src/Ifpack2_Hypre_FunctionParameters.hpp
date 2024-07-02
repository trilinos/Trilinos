// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef IFPACK2_HYPRE_FUNCTIONPARAMETERS_HPP
#define IFPACK2_HYPRE_FUNCTIONPARAMETERS_HPP

#include "Ifpack2_ConfigDefs.hpp"
#if defined(HAVE_IFPACK2_HYPRE) && defined(HAVE_IFPACK2_MPI)

#include <sstream>
#include "HYPRE_utilities.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "krylov.h"
#include "_hypre_parcsr_mv.h"
#include "_hypre_IJ_mv.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE.h"

// Hypre forward declarations (to avoid downstream header pollution)
struct hypre_IJMatrix_struct;
typedef struct hypre_IJMatrix_struct *HYPRE_IJMatrix;
struct hypre_IJVector_struct;
typedef struct hypre_IJVector_struct *HYPRE_IJVector;
struct hypre_ParCSRMatrix_struct;
typedef struct hypre_ParCSRMatrix_struct* HYPRE_ParCSRMatrix;
struct hypre_ParVector_struct;
typedef struct hypre_ParVector_struct * HYPRE_ParVector;
struct hypre_Solver_struct;
typedef struct hypre_Solver_struct *HYPRE_Solver;
struct hypre_ParVector_struct;
typedef struct hypre_ParVector_struct hypre_ParVector;
//struct hypre_Vector;

#ifndef HYPRE_ENUMS
#define HYPRE_ENUMS
//! This enumerated type defines the allowed solvers and preconditioners in Hypre. Some can be used as both solver and preconditioner.
  enum Hypre_Solver{
    BoomerAMG,
    ParaSails,
    Euclid,
    AMS,
    Hybrid,
    PCG,
    GMRES,
    FlexGMRES,
    LGMRES,
    BiCGSTAB
  };

  //! This enumerated type defines the two options for applying inverse, either solve or apply the preconditioner.
  enum Hypre_Chooser{
    Hypre_Is_Solver,
    Hypre_Is_Preconditioner
  };
#endif //HYPRE_ENUMS

// The Python script that generates the ParameterMap needs to be after these typedefs
typedef HYPRE_Int (*int_func)(HYPRE_Solver, HYPRE_Int);
typedef HYPRE_Int (*double_func)(HYPRE_Solver, HYPRE_Real);
typedef HYPRE_Int (*double_int_func)(HYPRE_Solver, HYPRE_Real, HYPRE_Int);
typedef HYPRE_Int (*int_double_func)(HYPRE_Solver, HYPRE_Int, HYPRE_Real);
typedef HYPRE_Int (*int_int_func)(HYPRE_Solver, HYPRE_Int, HYPRE_Int);
typedef HYPRE_Int (*int_star_func)(HYPRE_Solver, HYPRE_Int*);
typedef HYPRE_Int (*int_star_star_func)(HYPRE_Solver, HYPRE_Int**);
typedef HYPRE_Int (*double_star_func)(HYPRE_Solver, HYPRE_Real*);
typedef HYPRE_Int (*int_int_double_double_func)(HYPRE_Solver, HYPRE_Int, HYPRE_Int, HYPRE_Real, HYPRE_Real);
typedef HYPRE_Int (*int_int_int_double_int_int_func)(HYPRE_Solver, HYPRE_Int, HYPRE_Int, HYPRE_Int, HYPRE_Real, HYPRE_Int, HYPRE_Int);
typedef HYPRE_Int (*char_star_func)(HYPRE_Solver, char*);


namespace Ifpack2 {

  void IFPACK2_CHK_ERRV(int code);

  void IFPACK2_CHK_ERR(int code);

  //! This class is used to help with passing parameters in the SetParameter() function. Use this class to call Hypre's internal parameters.
  class FunctionParameter {
  public:
    //! Single int constructor.
    FunctionParameter(Hypre_Chooser chooser, int_func funct, HYPRE_Int param1) :
      chooser_(chooser),
      option_(0),
      int_func_(funct),
      int_param1_(param1) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, HYPRE_Int param1) :
      chooser_(chooser),
      option_(0),
      int_func_(hypreMapIntFunc_.at(funct_name)),
      int_param1_(param1) {}

    //! Single HYPRE_Real constructor.
    FunctionParameter(Hypre_Chooser chooser, double_func funct, HYPRE_Real param1):
      chooser_(chooser),
      option_(1),
      double_func_(funct),
      double_param1_(param1) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, HYPRE_Real param1):
      chooser_(chooser),
      option_(1),
      double_func_(hypreMapDoubleFunc_.at(funct_name)),
      double_param1_(param1) {}

    //! Single double, single int constructor.
    FunctionParameter(Hypre_Chooser chooser, double_int_func funct, HYPRE_Real param1, HYPRE_Int param2):
      chooser_(chooser),
      option_(2),
      double_int_func_(funct),
      int_param1_(param2),
      double_param1_(param1) {}

    //! Single int, single HYPRE_Real constructor.
    FunctionParameter(Hypre_Chooser chooser, int_double_func funct, HYPRE_Int param1, HYPRE_Real param2):
      chooser_(chooser),
      option_(10),
      int_double_func_(funct),
      int_param1_(param1),
      double_param1_(param2) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, HYPRE_Real param1, HYPRE_Int param2):
      chooser_(chooser),
      option_(2),
      double_int_func_(hypreMapDoubleIntFunc_.at(funct_name)),
      int_param1_(param2),
      double_param1_(param1) {}

    //! Two ints constructor.
    FunctionParameter(Hypre_Chooser chooser, int_int_func funct, HYPRE_Int param1, HYPRE_Int param2):
      chooser_(chooser),
      option_(3),
      int_int_func_(funct),
      int_param1_(param1),
      int_param2_(param2) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, HYPRE_Int param1, HYPRE_Int param2):
      chooser_(chooser),
      option_(3),
      int_int_func_(hypreMapIntIntFunc_.at(funct_name)),
      int_param1_(param1),
      int_param2_(param2) {}

    //! Int pointer constructor.
    FunctionParameter(Hypre_Chooser chooser, int_star_func funct, HYPRE_Int *param1):
      chooser_(chooser),
      option_(4),
      int_star_func_(funct),
      int_star_param_(param1) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, HYPRE_Int *param1):
      chooser_(chooser),
      option_(4),
      int_star_func_(hypreMapIntStarFunc_.at(funct_name)),
      int_star_param_(param1) {}

    //! HYPRE_Real pointer constructor.
    FunctionParameter(Hypre_Chooser chooser, double_star_func funct, double* param1):
      chooser_(chooser),
      option_(5),
      double_star_func_(funct),
      double_star_param_(param1) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, double* param1):
      chooser_(chooser),
      option_(5),
      double_star_func_(hypreMapDoubleStarFunc_.at(funct_name)),
      double_star_param_(param1) {}

    //! Two ints, two doubles constructor.
    FunctionParameter(Hypre_Chooser chooser, int_int_double_double_func funct, HYPRE_Int param1, HYPRE_Int param2, HYPRE_Real param3, HYPRE_Real param4):
      chooser_(chooser),
      option_(6),
      int_int_double_double_func_(funct),
      int_param1_(param1),
      int_param2_(param2),
      double_param1_(param3),
      double_param2_(param4) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, HYPRE_Int param1, HYPRE_Int param2, HYPRE_Real param3, HYPRE_Real param4):
      chooser_(chooser),
      option_(6),
      int_int_double_double_func_(hypreMapIntIntDoubleDoubleFunc_.at(funct_name)),
      int_param1_(param1),
      int_param2_(param2),
      double_param1_(param3),
      double_param2_(param4) {}

    //! Integer pointer to list of integer pointers
    FunctionParameter(Hypre_Chooser chooser, int_star_star_func funct, HYPRE_Int ** param1):
      chooser_(chooser),
      option_(7),
      int_star_star_func_(funct),
      int_star_star_param_(param1) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, HYPRE_Int** param1):
      chooser_(chooser),
      option_(7),
      int_star_star_func_(hypreMapIntStarStarFunc_.at(funct_name)),
      int_star_star_param_(param1) {}

    //! Five ints, one HYPRE_Real constructor.
    FunctionParameter(Hypre_Chooser chooser, int_int_int_double_int_int_func funct, HYPRE_Int param1, HYPRE_Int param2, HYPRE_Int param3, HYPRE_Real param4, HYPRE_Int param5, HYPRE_Int param6):
      chooser_(chooser),
      option_(8),
      int_int_int_double_int_int_func_(funct),
      int_param1_(param1),
      int_param2_(param2),
      int_param3_(param3),
      int_param4_(param5),
      int_param5_(param6),
      double_param1_(param4) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, HYPRE_Int param1, HYPRE_Int param2, HYPRE_Int param3, HYPRE_Real param4, HYPRE_Int param5, HYPRE_Int param6):
      chooser_(chooser),
      option_(8),
      int_int_int_double_int_int_func_(hypreMapIntIntIntDoubleIntIntFunc_.at(funct_name)),
      int_param1_(param1),
      int_param2_(param2),
      int_param3_(param3),
      int_param4_(param5),
      int_param5_(param6),
      double_param1_(param4) {}

    //! Char pointer constructor.
    FunctionParameter(Hypre_Chooser chooser, char_star_func funct, char *param1):
      chooser_(chooser),
      option_(9),
      char_star_func_(funct),
      char_star_param_(param1) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, char *param1):
      chooser_(chooser),
      option_(9),
      char_star_func_(hypreMapCharStarFunc_.at(funct_name)),
      char_star_param_(param1) {}

    //! Only method of this class. Calls the function pointer with the passed in HYPRE_Solver
    int CallFunction(HYPRE_Solver solver, HYPRE_Solver precond) {
      if(chooser_ == Hypre_Is_Solver){
        if(option_ == 0){
          return int_func_(solver, int_param1_);
        } else if(option_ == 1){
          return double_func_(solver, double_param1_);
        } else if(option_ == 2){
          return double_int_func_(solver, double_param1_, int_param1_);
        } else if (option_ == 3){
          return int_int_func_(solver, int_param1_, int_param2_);
        } else if (option_ == 4){
          return int_star_func_(solver, int_star_param_);
        } else if (option_ == 5){
          return double_star_func_(solver, double_star_param_);
        } else if (option_ == 6) {
          return int_int_double_double_func_(solver, int_param1_, int_param2_, double_param1_, double_param2_);
        } else if (option_ == 7) {
          return int_star_star_func_(solver, int_star_star_param_);
        } else if (option_ == 8) {
          return int_int_int_double_int_int_func_(solver, int_param1_, int_param2_, int_param3_, double_param1_, int_param4_, int_param5_);
        } else if (option_ == 9) {
          return char_star_func_(solver, char_star_param_);
        } else if (option_ == 10) {
          return int_double_func_(solver, int_param1_, double_param1_);
        } else {
          IFPACK2_CHK_ERR(-2);
        }
      } else {
        if(option_ == 0){
          return int_func_(precond, int_param1_);
        } else if(option_ == 1){
          return double_func_(precond, double_param1_);
        } else if(option_ == 2){
          return double_int_func_(precond, double_param1_, int_param1_);
        } else if(option_ == 3) {
          return int_int_func_(precond, int_param1_, int_param2_);
        } else if(option_ == 4) {
          return int_star_func_(precond, int_star_param_);
        } else if(option_ == 5) {
          return double_star_func_(precond, double_star_param_);
        } else if (option_ == 6) {
          return int_int_double_double_func_(precond, int_param1_, int_param2_, double_param1_, double_param2_);
        } else if (option_ == 7) {
          return int_star_star_func_(precond, int_star_star_param_);
        } else if (option_ == 8) {
          return int_int_int_double_int_int_func_(precond, int_param1_, int_param2_, int_param3_, double_param1_, int_param4_, int_param5_);
        } else if (option_ == 9) {
          return char_star_func_(solver, char_star_param_);
        } else if (option_ == 10) {
          return int_double_func_(precond, int_param1_, double_param1_);
        } else {
          IFPACK2_CHK_ERR(-2);
        }
      }
      return 0;
    }

    static bool isFuncIntInt(std::string funct_name) {
      return (hypreMapIntIntFunc_.find(funct_name) != hypreMapIntIntFunc_.end());
    }

    static bool isFuncIntIntDoubleDouble(std::string funct_name) {
      return (hypreMapIntIntDoubleDoubleFunc_.find(funct_name) != hypreMapIntIntDoubleDoubleFunc_.end());
    }

    static bool isFuncIntIntIntDoubleIntInt(std::string funct_name) {
      return (hypreMapIntIntIntDoubleIntIntFunc_.find(funct_name) != hypreMapIntIntIntDoubleIntIntFunc_.end());
    }

    static bool isFuncIntStarStar(std::string funct_name) {
      return (hypreMapIntStarStarFunc_.find(funct_name) != hypreMapIntStarStarFunc_.end());
    }

  private:
    Hypre_Chooser chooser_;
    int option_;
    int_func int_func_;
    double_func double_func_;
    double_int_func double_int_func_;
    int_double_func int_double_func_;
    int_int_func int_int_func_;
    int_star_func int_star_func_;
    double_star_func double_star_func_;
    int_int_double_double_func int_int_double_double_func_;
    int_int_int_double_int_int_func int_int_int_double_int_int_func_;
    int_star_star_func int_star_star_func_;
    char_star_func char_star_func_;
    HYPRE_Int int_param1_;
    HYPRE_Int int_param2_;
    HYPRE_Int int_param3_;
    HYPRE_Int int_param4_;
    HYPRE_Int int_param5_;
    HYPRE_Real double_param1_;
    HYPRE_Real double_param2_;
    HYPRE_Int *int_star_param_;
    HYPRE_Int **int_star_star_param_;
    HYPRE_Real *double_star_param_;
    char *char_star_param_;

    static const std::map<std::string, int_func> hypreMapIntFunc_;
    static const std::map<std::string, double_func> hypreMapDoubleFunc_;
    static const std::map<std::string, double_int_func> hypreMapDoubleIntFunc_;
    static const std::map<std::string, int_double_func> hypreMapIntDoubleFunc_;
    static const std::map<std::string, int_int_func> hypreMapIntIntFunc_;
    static const std::map<std::string, int_star_func> hypreMapIntStarFunc_;
    static const std::map<std::string, double_star_func> hypreMapDoubleStarFunc_;
    static const std::map<std::string, int_int_double_double_func> hypreMapIntIntDoubleDoubleFunc_;
    static const std::map<std::string, int_int_int_double_int_int_func> hypreMapIntIntIntDoubleIntIntFunc_;
    static const std::map<std::string, int_star_star_func> hypreMapIntStarStarFunc_;
    static const std::map<std::string, char_star_func> hypreMapCharStarFunc_;

  };

}

#endif // HAVE_IFPACK2_HYPRE && HAVE_IFPACK2_MPI

#endif /* IFPACK2_HYPRE_FUNCTIONPARAMETERS_HPP */
