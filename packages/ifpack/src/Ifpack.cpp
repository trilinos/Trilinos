/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/
#include "Ifpack_ConfigDefs.h"
#include "Ifpack.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_BlockRelaxation.h"
#include "Ifpack_IC.h"
#include "Ifpack_ICT.h"
#include "Ifpack_ILU.h"
#include "Ifpack_SILU.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_SPARSKIT.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_DenseContainer.h"
#include "Ifpack_TriDiContainer.h"
#include "Ifpack_SparseContainer.h"
#ifdef HAVE_IFPACK_AMESOS
#include "Ifpack_Amesos.h"
#endif
#ifdef HAVE_IFPACK_HIPS
#include "Ifpack_HIPS.h"
#endif
#ifdef HAVE_IFPACK_SUPERLU
#include "Ifpack_SILU.h"
#endif
#ifdef HAVE_IFPACK_SUPPORTGRAPH
#include "Ifpack_SupportGraph.h"
#endif

#include "Ifpack_Chebyshev.h"
#include "Ifpack_Polynomial.h"
#include "Ifpack_Krylov.h"
#include "Ifpack_IHSS.h"
#include "Ifpack_SORa.h"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StringToIntMap.hpp"
#include "Epetra_CrsMatrix.h"


namespace {

const Teuchos::StringToIntMap
precTypeNameToIntMap(
  "parameter \"Prec Type\"", Ifpack::numPrecTypes, Ifpack::precTypeNames
  );

} // namespace

//==============================================================================
const Ifpack::EPrecType Ifpack::precTypeValues[Ifpack::numPrecTypes] =
{
  POINT_RELAXATION
  ,POINT_RELAXATION_STAND_ALONE
  ,BLOCK_RELAXATION
  ,BLOCK_RELAXATION_STAND_ALONE
  ,BLOCK_RELAXATION_STAND_ALONE_ILU
  ,BLOCK_RELAXATION_STAND_ALONE_ILUT
  ,BLOCK_RELAXATION_STAND_ALONE_IC
#ifdef HAVE_IFPACK_SUPERLU
  ,BLOCK_RELAXATION_STAND_ALONE_SILU
#endif
#ifdef HAVE_IFPACK_AMESOS
  ,BLOCK_RELAXATION_STAND_ALONE_AMESOS
  ,BLOCK_RELAXATION_AMESOS
  ,AMESOS
  ,AMESOS_STAND_ALONE
#endif // HAVE_IFPACK_AMESOS
  ,IC
  ,IC_STAND_ALONE
  ,ICT
  ,ICT_STAND_ALONE
  ,ILU
  ,ILU_STAND_ALONE
  ,ILUT
  ,ILUT_STAND_ALONE
#ifdef HAVE_IFPACK_SPARSKIT
  ,SPARSKIT
#endif // HAVE_IFPACK_SPARSKIT
#ifdef HAVE_IFPACK_HIPS
  ,HIPS
#endif
#ifdef HAVE_HYPRE
  ,HYPRE
#endif
#ifdef HAVE_IFPACK_SUPERLU
  ,SILU
#endif
#if defined (HAVE_IFPACK_SUPPORTGRAPH) && defined (HAVE_IFPACK_AMESOS)
  ,MSF_AMESOS
#endif
#ifdef HAVE_IFPACK_SUPPORTGRAPH
  ,MSF_IC
#endif
  ,CHEBYSHEV
  ,POLYNOMIAL
  ,KRYLOV
  ,IHSS
  ,SORA
  ,TRIDI_RELAXATION
  ,TRIDI_RELAXATION_STAND_ALONE
};

//==============================================================================
const char* Ifpack::precTypeNames[Ifpack::numPrecTypes] =
{
  "point relaxation"
  ,"point relaxation stand-alone"
  ,"block relaxation"
  ,"block relaxation stand-alone"
  ,"block relaxation stand-alone (ILU)"
  ,"block relaxation stand-alone (ILUT)"
  ,"block relaxation stand-alone (IC)"
#ifdef HAVE_IFPACK_SUPERLU
  ,"block relaxation stand-alone (SILU)"
#endif
#ifdef HAVE_IFPACK_AMESOS
  ,"block relaxation stand-alone (Amesos)"
  ,"block relaxation (Amesos)"
  ,"Amesos"
  ,"Amesos stand-alone"
#endif
  ,"IC"
  ,"IC stand-alone"
  ,"ICT"
  ,"ICT stand-alone"
  ,"ILU"
  ,"ILU stand-alone"
  ,"ILUT"
  ,"ILUT stand-alone"
#ifdef HAVE_IFPACK_SPARSKIT
  ,"SPARSKIT"
#endif
#ifdef HAVE_IFPACK_HIPS
  ,"HIPS"
#endif
#ifdef HAVE_HYPRE
  ,"Hypre"
#endif
#ifdef HAVE_IFPACK_SUPERLU
  ,"SILU"
#endif
#if defined (HAVE_IFPACK_SUPPORTGRAPH) && defined (HAVE_IFPACK_AMESOS)
  ,"MSF Amesos"
#endif
#ifdef HAVE_IFPACK_SUPPORTGRAPH
  ,"MSF IC"
#endif
  ,"Chebyshev"
  ,"Polynomial"
  ,"Krylov"
  ,"IHSS"
  ,"SORa"
  ,"tridi relaxation"
  ,"tridi relaxation stand-alone"
};

//==============================================================================
const bool Ifpack::supportsUnsymmetric[Ifpack::numPrecTypes] =
{
  true // point relaxation
  ,true // point relaxation stand-alone
  ,true // block relaxation
  ,true // block relaxation stand-alone
  ,true // block relaxation stand-alone (ILU)
  ,true // block relaxation stand-alone (ILUT)
  ,false // block relaxation stand-alone (IC)
#ifdef HAVE_IFPACK_SUPERLU
  ,true // block relaxation stand-alone (SILU)
#endif
#ifdef HAVE_IFPACK_AMESOS
  ,true // block relaxation stand-alone (Amesos)
  ,true // block relaxation (Amesos)
  ,true // Amesos
  ,true // Amesos stand-alone 
#endif
  ,false // IC
  ,false // IC stand-alone
  ,false // ICT
  ,false // ICT stand-alone
  ,true // ILU
  ,true // ILU stand-alone
  ,true // ILUT
  ,true // ILUT stand-alone
#ifdef HAVE_IFPACK_SPARSKIT
  ,true // SPARSKIT
#endif
#ifdef HAVE_IFPACK_HIPS
  ,true // HIPS
#endif  
#ifdef HAVE_HYPRE
  ,true
#endif
#ifdef HAVE_IFPACK_SUPERLU
  ,true // SuperLU's Supernodal ILUTP
#endif
#if defined (HAVE_IFPACK_SUPPORTGRAPH) && defined (HAVE_IFPACK_AMESOS)
  ,false
#endif
#ifdef HAVE_IFPACK_SUPPORTGRAPH
  ,false
#endif
  ,false // CHEBYSHEV
  ,true  // POLYNOMIAL
  ,true  // KRYLOV
  ,true  // IHSS
  ,true  // SORa
  ,true  // tridi relaxation
  ,true  // tridi relaxation standalone
};

//==============================================================================
Ifpack_Preconditioner* Ifpack::Create(EPrecType PrecType,
                                      Epetra_RowMatrix* Matrix,
                                      const int Overlap,
                                      bool overrideSerialDefault)
{
  const bool serial = (Matrix->Comm().NumProc() == 1);

  switch(PrecType) {
    case POINT_RELAXATION:
      if (serial && !overrideSerialDefault)
        return(new Ifpack_PointRelaxation(Matrix));
      else
        return(new Ifpack_AdditiveSchwarz<Ifpack_PointRelaxation>(Matrix, Overlap));
    case POINT_RELAXATION_STAND_ALONE:
      return(new Ifpack_PointRelaxation(Matrix));
    case BLOCK_RELAXATION:
      if (serial && !overrideSerialDefault)
        return(new Ifpack_BlockRelaxation<Ifpack_DenseContainer>(Matrix));
      else
        return(new Ifpack_AdditiveSchwarz<
             Ifpack_BlockRelaxation<Ifpack_DenseContainer> >(Matrix,Overlap));
    case BLOCK_RELAXATION_STAND_ALONE:
      return(new Ifpack_BlockRelaxation<Ifpack_DenseContainer>(Matrix));
    case BLOCK_RELAXATION_STAND_ALONE_ILU:
      return(new Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_ILU> >(Matrix));
    case BLOCK_RELAXATION_STAND_ALONE_ILUT:
      return(new Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_ILUT> >(Matrix));
    case BLOCK_RELAXATION_STAND_ALONE_IC:
      return(new Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_IC> >(Matrix));
#ifdef HAVE_IFPACK_SUPERLU
    case BLOCK_RELAXATION_STAND_ALONE_SILU:
      return(new Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_SILU> >(Matrix));
#endif
#ifdef HAVE_IFPACK_AMESOS
    case BLOCK_RELAXATION_STAND_ALONE_AMESOS:
      return(new Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_Amesos> >(Matrix));
    case BLOCK_RELAXATION_AMESOS:
      return(new Ifpack_AdditiveSchwarz<
             Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_Amesos> > >(Matrix,Overlap));
    case AMESOS:
      if (serial && !overrideSerialDefault)
        return(new Ifpack_Amesos(Matrix));
      else
        return(new Ifpack_AdditiveSchwarz<Ifpack_Amesos>(Matrix,Overlap));
    case AMESOS_STAND_ALONE:
      return(new Ifpack_Amesos(Matrix));
#endif
    case IC:
      if (serial && !overrideSerialDefault)
        return(new Ifpack_IC(Matrix));
      else
        return(new Ifpack_AdditiveSchwarz<Ifpack_IC>(Matrix,Overlap));
    case IC_STAND_ALONE:
      return(new Ifpack_IC(Matrix));
    case ICT:
      if (serial && !overrideSerialDefault)
        return(new Ifpack_ICT(Matrix));
      else
        return(new Ifpack_AdditiveSchwarz<Ifpack_ICT>(Matrix,Overlap));
    case ICT_STAND_ALONE:
      return(new Ifpack_ICT(Matrix));
    case ILU:
      if (serial && !overrideSerialDefault)
        return(new Ifpack_ILU(Matrix));
      else
        return(new Ifpack_AdditiveSchwarz<Ifpack_ILU>(Matrix,Overlap));
    case ILU_STAND_ALONE:
      return(new Ifpack_ILU(Matrix));
    case ILUT:
      if (serial && !overrideSerialDefault)
        return(new Ifpack_ILUT(Matrix));
      else
        return(new Ifpack_AdditiveSchwarz<Ifpack_ILUT>(Matrix,Overlap));
    case ILUT_STAND_ALONE:
      return(new Ifpack_ILUT(Matrix));
#ifdef HAVE_IFPACK_SPARSKIT
    case SPARSKIT:
      return(new Ifpack_SPARSKIT(Matrix));
#endif
#ifdef HAVE_IFPACK_HIPS
    case HIPS:      
      return(new Ifpack_HIPS(Matrix));
#endif      
#ifdef HAVE_HYPRE
    case HYPRE:
      return(new Ifpack_Hypre(Matrix));
#endif
#ifdef HAVE_IFPACK_SUPERLU
    case SILU:
      return(new Ifpack_SILU(Matrix));
#endif
#if defined (HAVE_IFPACK_SUPPORTGRAPH) && defined (HAVE_IFPACK_AMESOS)
    case MSF_AMESOS:
      if (serial && !overrideSerialDefault)
	return(new Ifpack_SupportGraph<Ifpack_Amesos>(Matrix));
      else
	return(new Ifpack_AdditiveSchwarz<Ifpack_SupportGraph<Ifpack_Amesos> >(Matrix,Overlap));
#endif
#ifdef HAVE_IFPACK_SUPPORTGRAPH
    case MSF_IC:
      if (serial && !overrideSerialDefault)
	return(new Ifpack_SupportGraph<Ifpack_SupportGraph<Ifpack_IC> >(Matrix));
      else
	return(new Ifpack_AdditiveSchwarz<Ifpack_SupportGraph<Ifpack_IC> >(Matrix,Overlap));
#endif
    case CHEBYSHEV:
      return(new Ifpack_Chebyshev(Matrix));
    case POLYNOMIAL:
      return(new Ifpack_Polynomial(Matrix));
    case KRYLOV:
      if (serial && !overrideSerialDefault)
        return(new Ifpack_Krylov(Matrix));
      else
        return(new Ifpack_AdditiveSchwarz<Ifpack_Krylov>(Matrix, Overlap));
#ifdef HAVE_IFPACK_EPETRAEXT
    case IHSS:
      return(new Ifpack_IHSS(Matrix));  
    case SORA:
      return(new Ifpack_SORa(Matrix));  
#endif
    case TRIDI_RELAXATION:
     if (serial && !overrideSerialDefault)
      return(new Ifpack_BlockRelaxation<Ifpack_TriDiContainer>(Matrix));
     else
      return(new Ifpack_AdditiveSchwarz<
             Ifpack_BlockRelaxation<Ifpack_TriDiContainer> >(Matrix,Overlap));
    case TRIDI_RELAXATION_STAND_ALONE:
      return(new Ifpack_BlockRelaxation<Ifpack_TriDiContainer>(Matrix));
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
      // The only way to get here is if some code developer does a cast like
      // (EPrecType)(anyNumber).  You can never get here by passing in a
      // string value for the preconditioner!
  } // end switch
  return 0; // This will never ever be executed!
}

//==============================================================================
Ifpack_Preconditioner* Ifpack::Create(const string PrecType,
                                      Epetra_RowMatrix* Matrix,
                                      const int Overlap,
                                      bool overrideSerialDefault)
{
  try {
    return Ifpack::Create(Teuchos::get<EPrecType>(::precTypeNameToIntMap,PrecType),Matrix,Overlap,overrideSerialDefault);
  }
  catch( const Teuchos::StringToIntMap::DoesNotExist &excpt ) {
    // The old implementation of this function just silently returned a NULL
    // when a preconditiner type name was not recognized.  If you like this
    // behavior then you should use this function.  If you do not like this
    // behavior, then consider using the Ifpack/Thyra adapter
    // Thyra::IfpackPreconditionerFactory or better yet the Stratimikos
    // wrapper class Stratimikos::DefaultLinearSolverBuilder.
  }
  return 0;
}

// ======================================================================
int Ifpack::SetParameters(int argc, char* argv[],
                          Teuchos::ParameterList& List, string& PrecType,
                          int& Overlap)
{
  // THIS FUNCTION IS VERY INCOMPLETE...

  Teuchos::CommandLineProcessor CLP;

  // prec type
  string ifp_prec_type = "ILU";
  CLP.setOption("ifp-prec-type",&ifp_prec_type,"Preconditioner type");
  // overlap among the processors
  int ifp_overlap = 0;
  CLP.setOption("ifp-overlap",&ifp_overlap,"Overlap among processors");
  // relaxation type
  string ifp_relax_type = "Jacobi";
  CLP.setOption("ifp-relax-type",&ifp_relax_type,"Relaxation type");
  // sweeps (for relax only)
  int ifp_relax_sweeps = 1;
  CLP.setOption("ifp-relax-sweeps",
                &ifp_relax_sweeps,"Number of sweeps for relaxation");
  // damping (for relax only)
  double ifp_relax_damping = 1.0;
  CLP.setOption("ifp-relax-damping",
                &ifp_relax_damping,"Damping for relaxation");
  // partitioner type (for block relaxation only)
  string ifp_part_type = "greedy";
  CLP.setOption("ifp-part-type",&ifp_part_type,"Partitioner type");
  // number of local parts (for block relaxation only)
  int ifp_part_local = 1;
  CLP.setOption("ifp-part-local",&ifp_part_local,"number of local partitions");

  // allow users to specify other options for other packages
  CLP.recogniseAllOptions(false);
  CLP.throwExceptions(false);
  CLP.parse(argc,argv);

  // I cannot really set those in the List, I pass them back to the user
  PrecType = ifp_prec_type;
  Overlap = ifp_overlap;

  // set the list here
  List.set("relaxation: type", ifp_relax_type);
  List.set("relaxation: sweeps", ifp_relax_sweeps);
  List.set("relaxation: damping factor", ifp_relax_damping);
  List.set("partitioner: type", ifp_part_type);
  List.set("partitioner: local parts", ifp_part_local);

  return(0);
}
