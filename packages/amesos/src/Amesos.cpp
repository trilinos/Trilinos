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

#include "Amesos_config.h"
#include "Amesos.h"
#include "Amesos_Klu.h"
#ifdef HAVE_AMESOS_LAPACK
#include "Amesos_Lapack.h"
#endif
#if defined(HAVE_AMESOS_MUMPS) && defined(HAVE_MPI)
#include "Amesos_Mumps.h"
#endif
#ifdef HAVE_AMESOS_SCALAPACK
#include "Amesos_Scalapack.h"
#endif
#ifdef HAVE_AMESOS_UMFPACK
#include "Amesos_Umfpack.h"
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
#include "Amesos_Superludist.h"
#endif
#ifdef HAVE_AMESOS_SUPERLU
#include "Amesos_Superlu.h"
#endif
#ifdef HAVE_AMESOS_DSCPACK
#include "Amesos_Dscpack.h"
#endif
#if defined(HAVE_AMESOS_PARDISO) || defined(HAVE_AMESOS_PARDISO_MKL)
#include "Amesos_Pardiso.h"
#endif
#ifdef HAVE_AMESOS_TAUCS
#include "Amesos_Taucs.h"
#endif
#ifdef HAVE_AMESOS_CSS_MKL
#include "Amesos_CssMKL.h"
#endif
#ifdef HAVE_AMESOS_PARAKLETE
#include "Amesos_Paraklete.h"
#endif
#ifdef HAVE_AMESOS_CSPARSE
#include "Amesos_CSparse.h"
#endif
#include "Epetra_Object.h"

static bool verbose = false; 

Amesos_BaseSolver* Amesos::Create(const char* ClassType, 
				  const Epetra_LinearProblem& LinearProblem ) 
{ 
  std::string CT = ClassType; 
  return(Create(CT,LinearProblem));
}

Amesos_BaseSolver* Amesos::Create(const std::string CT,
				  const Epetra_LinearProblem& LinearProblem )
{

  if ((CT == "Amesos_Lapack") || (CT == "Lapack")) { 
#ifdef HAVE_AMESOS_LAPACK
    return new Amesos_Lapack(LinearProblem); 
#else
    if (verbose) std::cerr << "Amesos_Lapack is not implemented" << std::endl ; 
    return(0); 
#endif
  } 
  
  if ((CT == "Amesos_Klu") || (CT == "Klu")) { 
#ifdef HAVE_AMESOS_KLU
    return new Amesos_Klu(LinearProblem); 
#else
    if (verbose) std::cerr << "Amesos_Klu is not implemented" << std::endl ; 
    return(0); 
#endif
  } 
  
  if ((CT == "Amesos_Umfpack") || (CT == "Umfpack")) { 
#ifdef HAVE_AMESOS_UMFPACK
    return new Amesos_Umfpack(LinearProblem); 
#else
    if (verbose) std::cerr << "Amesos_Umfpack is not implemented" << std::endl ; 
    return(0); 
#endif
  } 
  
  if ((CT == "Amesos_Superlu") || (CT == "Superlu")) { 
#ifdef HAVE_AMESOS_SUPERLU
    return new Amesos_Superlu(LinearProblem); 
#else
    if (verbose) std::cerr << "Amesos_Superlu is not implemented" << std::endl ; 
    return(0); 
#endif
  } 
  
  if ((CT == "Amesos_Superludist") || (CT == "Superludist")) { 
#ifdef HAVE_AMESOS_SUPERLUDIST
    return new Amesos_Superludist(LinearProblem); 
#else
    if (verbose) std::cerr << "Amesos_Superludist is not implemented" << std::endl ; 
    return(0); 
#endif
  } 
  
  if ((CT == "Amesos_Mumps") || (CT == "Mumps")) { 
#if defined(HAVE_AMESOS_MUMPS) && defined(HAVE_MPI)
    return new Amesos_Mumps(LinearProblem); 
#else
    if (verbose) std::cerr << "Amesos_Mumps is not implemented" << std::endl ; 
    return(0); 
#endif
  } 
  
  if ((CT == "Amesos_Scalapack") || (CT == "Scalapack")) { 
#ifdef HAVE_AMESOS_SCALAPACK
    return new Amesos_Scalapack(LinearProblem); 
#else
    if (verbose) std::cerr << "Amesos_Scalapack is not implemented" << std::endl ; 
    return(0); 
#endif
  } 
  
  if ((CT == "Amesos_Dscpack") || (CT == "Dscpack")) { 
#ifdef HAVE_AMESOS_DSCPACK
    return new Amesos_Dscpack(LinearProblem); 
#else
    if (verbose) std::cerr << "Amesos_Dscpack is not implemented" << std::endl ; 
    return(0); 
#endif
  } 
  
  if ((CT == "Amesos_Pardiso") || (CT == "Pardiso")) { 
#if defined(HAVE_AMESOS_PARDISO) || defined(HAVE_AMESOS_PARDISO_MKL)
    return new Amesos_Pardiso(LinearProblem); 
#else
    if (verbose) std::cerr << "Amesos_Pardiso is not implemented" << std::endl ; 
    return(0); 
#endif
  } 
  
  if ((CT == "Amesos_CssMKL") || (CT == "CssMKL")) {
#ifdef HAVE_AMESOS_CSS_MKL
    return new Amesos_CssMKL(LinearProblem);
#else
    if (verbose) std::cerr << "Amesos_CssMKL is not implemented" << std::endl ;
    return(0);
#endif
  }

  if ((CT == "Amesos_Paraklete") || (CT == "Paraklete")) { 
#ifdef HAVE_AMESOS_PARAKLETE
    return new Amesos_Paraklete(LinearProblem); 
#else
    if (verbose) std::cerr << "Amesos_Paraklete is not implemented" << std::endl ; 
    return(0); 
#endif
  }

  if ((CT == "Amesos_Taucs") || (CT == "Taucs")) { 
#ifdef HAVE_AMESOS_TAUCS
    return new Amesos_Taucs(LinearProblem); 
#else
    if (verbose) std::cerr << "Amesos_Taucs is not implemented" << std::endl ; 
    return(0); 
#endif
  } 
  
  if ((CT == "Amesos_CSparse") || (CT == "CSparse")) {
#ifdef HAVE_AMESOS_CSPARSE
    return new Amesos_CSparse(LinearProblem);
#else
    if (verbose) std::cerr << "Amesos_CSparse is not implemented" << std::endl ;
    return(0);
#endif
  }

  if (verbose) std::cerr << "Unknown class type:" << CT << std::endl ; 
  return(0); 
}

// ====================================================================
bool Amesos::Query(const char* ClassType)
{
  std::string CT = ClassType;
  return(Query(CT));
}

// ====================================================================
bool Amesos::Query(const std::string CT) 
{ 

  if ((CT == "Amesos_Lapack") || (CT == "Lapack")) { 
#ifdef HAVE_AMESOS_LAPACK
    return true;
#else
    return false;
#endif
  } 
  
  if ((CT == "Amesos_Klu") || (CT == "Klu")) { 
#ifdef HAVE_AMESOS_KLU
    return true; 
#else
    return false;
#endif
  } 
  
  if ((CT == "Amesos_Umfpack") || (CT == "Umfpack")) { 
#ifdef HAVE_AMESOS_UMFPACK
    return true;
#else
    return false;
#endif
  } 
  
  if ((CT == "Amesos_Superlu") || ( CT == "Superlu")) { 
#ifdef HAVE_AMESOS_SUPERLU
    return true; 
#else
    return false;
#endif
  }

  if ((CT == "Amesos_Superludist") || (CT == "Superludist")) { 
#ifdef HAVE_AMESOS_SUPERLUDIST
    return true;
#else
    return false;
#endif
  } 

  if ((CT == "Amesos_Mumps") || (CT == "Mumps")) { 
#ifdef HAVE_AMESOS_MUMPS
    return true;
#else
    return false;
#endif
  } 

  if ((CT == "Amesos_Scalapack") || (CT == "Scalapack")) { 
#ifdef HAVE_AMESOS_SCALAPACK
    return true;
#else
    return false;
#endif
  } 

  if ((CT == "Amesos_Dscpack") || (CT == "Dscpack")) { 
#ifdef HAVE_AMESOS_DSCPACK
    return true;
#else
    return false;
#endif
  } 
  
  if ((CT == "Amesos_Pardiso") || (CT == "Pardiso")) { 
#if defined(HAVE_AMESOS_PARDISO) || defined(HAVE_AMESOS_PARDISO_MKL)
    return true;
#else
    return false;
#endif
  } 
  
  if ((CT == "Amesos_CssMKL") || (CT == "Css_MKL")) {
#ifdef HAVE_AMESOS_CSS_MKL
    return true;
#else
    return false;
#endif
  }

  if ((CT == "Amesos_Taucs") || (CT == "Taucs")) { 
#ifdef HAVE_AMESOS_TAUCS
    return true;
#else
    return false;
#endif
  }

  if ((CT == "Amesos_Paraklete") || (CT == "Paraklete")) { 
#ifdef HAVE_AMESOS_PARAKLETE
    return true;
#else
    return false;
#endif
  } 

  if ((CT == "Amesos_CSparse") || (CT == "CSparse")) {
#ifdef HAVE_AMESOS_CSPARSE
    return true;
#else
    return false;
#endif
  }

  return(false);

}

Teuchos::ParameterList Amesos::GetValidParameters(){ 
  Teuchos::ParameterList ParamList  ;

  //  Status Parameters - see Amesos_Status.cpp 

  ParamList.set("OutputLevel",1 ) ; 
  ParamList.set("DebugLevel", 0 ) ; 
  ParamList.set("PrintTiming", false ) ; 
  ParamList.set("ComputeVectorNorms", false ) ; 
  ParamList.set("ComputeTrueResidual", false ) ; 


  //  Control Parameters - see Amesos_Control.cpp
  ParamList.set("AddZeroToDiag", false ) ; 
  ParamList.set("AddToDiag", 0.0 ) ; 
  ParamList.set("Refactorize", false ) ; 
  ParamList.set("RcondThreshold", 1e-12 ) ; 
  ParamList.set("MaxProcs", -1 ) ; 
  ParamList.set("MatrixProperty","general" ) ; 
  ParamList.set("ScaleMethod", 0 ) ; 
  ParamList.set("Reindex", false ) ; 


  //  Klu Parameters
  ParamList.set("TrustMe", false ) ;     //  If set, Amesos_Klu trusts that the data can be used in place - see Amesos_Klu.cpp

  //  Superlu Parameters - none 

  //  Dscpack Parameters - none

  //  Superludist Parameters
  ParamList.set("Redistribute", false ) ; 
  Teuchos::ParameterList SuperludistParams;
  { 
    SuperludistParams.set("ReuseSymbolic",false);
    SuperludistParams.set("Fact","SamePattern");
    SuperludistParams.set("Equil",false);
    SuperludistParams.set("ColPerm","NOT SET");
    SuperludistParams.set("RowPerm","NOT SET");
    SuperludistParams.set("perm_c",(int*) 0);
    SuperludistParams.set("perm_r",(int*) 0);
    SuperludistParams.set("IterRefine","NOT SET");
    SuperludistParams.set("ReplaceTinyPivot",true);
    SuperludistParams.set("PrintNonzeros",false);
  }
  ParamList.set("Superludist", SuperludistParams ) ; 
  //  MC64 Parameters - none

  //  Lapack Parameters
  Teuchos::ParameterList LapackParams;
  { 
    LapackParams.set("Equilibrate",true);
  }
  ParamList.set("Lapack", LapackParams ) ; 
  //  Mumps Parameters
  ParamList.set("NoDestroy",false);
  Teuchos::ParameterList MumpsParams;
  { 
    MumpsParams.set("Equilibrate",true);
    // ICNTL0, ICNT1, ..., ICNTL40 
    for (int i = 1 ; i <= 40 ; ++i)
    {
      char what[80];
      sprintf(what, "ICNTL(%d)", i);
      if (MumpsParams.isParameter(what)) 
        MumpsParams.set(what,0);
    }

    // CNTL0, CNTL1, ..., CNTL5
    for (int i = 1 ; i <= 5 ; ++i)
    {
      char what[80];
      sprintf(what, "CNTL(%d)", i);
      if (MumpsParams.isParameter(what)) 
	MumpsParams.set(what,0.0);
    }
    MumpsParams.set("RowScaling",(double *) 0);
    MumpsParams.set("ColScaling",(double *) 0);

  }
  ParamList.set("Mumps", MumpsParams ) ; 

  //  Paraklete Parameters - same as Klu

  //  Pardiso Parameters
  Teuchos::ParameterList PardisoParams;
  { 
    PardisoParams.set("MSGLVL",0);
    PardisoParams.set("IPARM(1)",0);
    PardisoParams.set("IPARM(2)",0);
    PardisoParams.set("IPARM(3)",0);
    PardisoParams.set("IPARM(4)",0);
    PardisoParams.set("IPARM(8)",0);
    PardisoParams.set("IPARM(10)",0);
    PardisoParams.set("IPARM(11)",0);
    PardisoParams.set("IPARM(18)",0);
    PardisoParams.set("IPARM(19)",0);
    PardisoParams.set("IPARM(21)",0);
 
  }
  ParamList.set("Pardiso", PardisoParams ) ; 

  //  Scalapack Parameters
  Teuchos::ParameterList ScalapackParams;
  { 
    ScalapackParams.set("2D distribution",true);
    ScalapackParams.set("grid_nb",32);
  }
  ParamList.set("Scalapack", ScalapackParams ) ; 
  //  Taucs Parameters - none

  //  Umfpack Parameters - none

  return ParamList ; 
} 
