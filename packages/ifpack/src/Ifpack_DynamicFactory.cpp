/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
//@HEADER
*/

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_DynamicFactory.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_BlockRelaxation.h"
#include "Ifpack_IC.h"
#include "Ifpack_ICT.h"
#include "Ifpack_ILU.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_SPARSKIT.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_DenseContainer.h"
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

#include "Ifpack_Chebyshev.h"
#include "Ifpack_IHSS.h"
#include "Ifpack_SORa.h"

#include "Teuchos_StringToIntMap.hpp"
#include "Epetra_CrsMatrix.h"

// Define builder functions for the default set of preconditioners
namespace {

// POINT_RELAXATION
Ifpack_Preconditioner* buildPointRelaxation(Epetra_RowMatrix* Matrix,
											int Overlap,
											bool serial,
											bool overrideSerialDefault)
{
  if (serial && !overrideSerialDefault)
	return(new Ifpack_PointRelaxation(Matrix));
  else
	return(new Ifpack_AdditiveSchwarz<Ifpack_PointRelaxation>(Matrix, Overlap));
}

// POINT_RELAXATION_STAND_ALONE
Ifpack_Preconditioner* buildPointRelaxationStandAlone(Epetra_RowMatrix* Matrix,
													  int /*Overlap*/,
													  bool /*serial*/,
													  bool /*overrideSerialDefault*/)
{
	return new Ifpack_PointRelaxation(Matrix);
}

// BLOCK_RELAXATION
Ifpack_Preconditioner* buildBlockRelaxation(Epetra_RowMatrix* Matrix,
											int Overlap,
											bool serial,
											bool overrideSerialDefault)
{
  if (serial && !overrideSerialDefault)
	return(new Ifpack_BlockRelaxation<Ifpack_DenseContainer>(Matrix));
  else
	return(new Ifpack_AdditiveSchwarz<
		 Ifpack_BlockRelaxation<Ifpack_DenseContainer> >(Matrix,Overlap));
}

// BLOCK_RELAXATION_STAND_ALONE
Ifpack_Preconditioner* buildBlockRelaxationStandAlone(Epetra_RowMatrix* Matrix,
													  int /*Overlap*/,
													  bool /*serial*/,
													  bool /*overrideSerialDefault*/)
{
	return new Ifpack_BlockRelaxation<Ifpack_DenseContainer>(Matrix);
}

// BLOCK_RELAXATION_STAND_ALONE_ILU
Ifpack_Preconditioner* buildBlockRelaxationStandAloneILU(Epetra_RowMatrix* Matrix,
													  	 int /*Overlap*/,
													  	 bool /*serial*/,
													  	 bool /*overrideSerialDefault*/)
{
	return new Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_ILU> >(Matrix);
}

#ifdef HAVE_IFPACK_AMESOS

// BLOCK_RELAXATION_STAND_ALONE_AMESOS
Ifpack_Preconditioner* buildBlockRelaxationStandAloneAmesos(Epetra_RowMatrix* Matrix,
													  	 	int /*Overlap*/,
													  	 	bool /*serial*/,
													  	 	bool /*overrideSerialDefault*/)
{
	return new Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_Amesos> >(Matrix);
}

// BLOCK_RELAXATION_AMESOS
Ifpack_Preconditioner* buildBlockRelaxationAmesos(Epetra_RowMatrix* Matrix,
												  int Overlap,
												  bool /*serial*/,
												  bool /*overrideSerialDefault*/)
{
	return(new Ifpack_AdditiveSchwarz<
					  Ifpack_BlockRelaxation<
							  Ifpack_SparseContainer<Ifpack_Amesos> > >
				(Matrix,Overlap));
}

// AMESOS
Ifpack_Preconditioner* buildAmesos(Epetra_RowMatrix* Matrix,
								   int Overlap,
								   bool serial,
								   bool overrideSerialDefault)
{
  if (serial && !overrideSerialDefault)
	return(new Ifpack_Amesos(Matrix));
  else
	return(new Ifpack_AdditiveSchwarz<Ifpack_Amesos>(Matrix,Overlap));
}

// AMESOS_STAND_ALONE
Ifpack_Preconditioner* buildAmesosStandAlone(Epetra_RowMatrix* Matrix,
											 int /*Overlap*/,
											 bool /*serial*/,
											 bool /*overrideSerialDefault*/)
{
	return new Ifpack_Amesos(Matrix);
}

#endif // HAVE_IFPACK_AMESOS

// IC
Ifpack_Preconditioner* buildIC(Epetra_RowMatrix* Matrix,
							   int Overlap,
							   bool serial,
							   bool overrideSerialDefault)
{
  if (serial && !overrideSerialDefault)
	return(new Ifpack_IC(Matrix));
  else
	return(new Ifpack_AdditiveSchwarz<Ifpack_IC>(Matrix,Overlap));
}

// IC_STAND_ALONE
Ifpack_Preconditioner* buildICStandAlone(Epetra_RowMatrix* Matrix,
										 int /*Overlap*/,
										 bool /*serial*/,
										 bool /*overrideSerialDefault*/)
{
	return new Ifpack_IC(Matrix);
}

// ICT
Ifpack_Preconditioner* buildICT(Epetra_RowMatrix* Matrix,
							    int Overlap,
							    bool serial,
							    bool overrideSerialDefault)
{
  if (serial && !overrideSerialDefault)
	return(new Ifpack_ICT(Matrix));
  else
	return(new Ifpack_AdditiveSchwarz<Ifpack_ICT>(Matrix,Overlap));
}

// ICT_STAND_ALONE
Ifpack_Preconditioner* buildICTStandAlone(Epetra_RowMatrix* Matrix,
										  int /*Overlap*/,
										  bool /*serial*/,
										  bool /*overrideSerialDefault*/)
{
	return new Ifpack_ICT(Matrix);
}

// ILU
Ifpack_Preconditioner* buildILU(Epetra_RowMatrix* Matrix,
							    int Overlap,
							    bool serial,
							    bool overrideSerialDefault)
{
  if (serial && !overrideSerialDefault)
	return(new Ifpack_ILU(Matrix));
  else
	return(new Ifpack_AdditiveSchwarz<Ifpack_ILU>(Matrix,Overlap));
}

// ILU_STAND_ALONE
Ifpack_Preconditioner* buildILUStandAlone(Epetra_RowMatrix* Matrix,
										  int /*Overlap*/,
										  bool /*serial*/,
										  bool /*overrideSerialDefault*/)
{
	return new Ifpack_ILU(Matrix);
}

// ILUT
Ifpack_Preconditioner* buildILUT(Epetra_RowMatrix* Matrix,
							     int Overlap,
							     bool serial,
							     bool overrideSerialDefault)
{
  if (serial && !overrideSerialDefault)
	return(new Ifpack_ILUT(Matrix));
  else
	return(new Ifpack_AdditiveSchwarz<Ifpack_ILUT>(Matrix,Overlap));
}

// ILUT_STAND_ALONE
Ifpack_Preconditioner* buildILUTStandAlone(Epetra_RowMatrix* Matrix,
										   int /*Overlap*/,
										   bool /*serial*/,
										   bool /*overrideSerialDefault*/)
{
	return new Ifpack_ILUT(Matrix);
}

#ifdef HAVE_IFPACK_SPARSKIT

// SPARSKIT
Ifpack_Preconditioner* buildSPARSKIT(Epetra_RowMatrix* Matrix,
									 int /*Overlap*/,
									 bool /*serial*/,
									 bool /*overrideSerialDefault*/)
{
	return new Ifpack_SPARSKIT(Matrix);
}

#endif // HAVE_IFPACK_SPARSKIT

#ifdef HAVE_IFPACK_HIPS

// HIPS
Ifpack_Preconditioner* buildHIPS(Epetra_RowMatrix* Matrix,
								 int /*Overlap*/,
								 bool /*serial*/,
								 bool /*overrideSerialDefault*/)
{
	return new Ifpack_HIPS(Matrix);
}

#endif // HAVE_IFPACK_HIPS

#ifdef HAVE_HYPRE

// HYPRE
Ifpack_Preconditioner* buildHypre(Epetra_RowMatrix* Matrix,
								  int /*Overlap*/,
								  bool /*serial*/,
								  bool /*overrideSerialDefault*/)
{
	return new Ifpack_Hypre(Matrix);
}

#endif // HAVE_IFPACK_HYPRE

#ifdef HAVE_IFPACK_SUPERLU

// SILU
Ifpack_Preconditioner* buildSILU(Epetra_RowMatrix* Matrix,
								 int /*Overlap*/,
								 bool /*serial*/,
								 bool /*overrideSerialDefault*/)
{
	return new Ifpack_SILU(Matrix);
}

#endif // HAVE_IFPACK_SUPERLU

// CHEBYSHEV
Ifpack_Preconditioner* buildChebyshev(Epetra_RowMatrix* Matrix,
								 	  int /*Overlap*/,
								 	  bool /*serial*/,
								 	  bool /*overrideSerialDefault*/)
{
	return new Ifpack_Chebyshev(Matrix);
}

#ifdef HAVE_IFPACK_EPETRAEXT

// IHSS
Ifpack_Preconditioner* buildIHSS(Epetra_RowMatrix* Matrix,
								 int /*Overlap*/,
								 bool /*serial*/,
								 bool /*overrideSerialDefault*/)
{
	return new Ifpack_IHSS(Matrix);
}

// SORA
Ifpack_Preconditioner* buildSORa(Epetra_RowMatrix* Matrix,
								 int /*Overlap*/,
								 bool /*serial*/,
								 bool /*overrideSerialDefault*/)
{
	return new Ifpack_SORa(Matrix);
}

#endif // HAVE_IFPACK_EPETRAEXT

}

std::map<std::string, Ifpack_DynamicFactory::builderFunction>
	Ifpack_DynamicFactory::PreconditionerMap_;
bool Ifpack_DynamicFactory::Initialized_ = false;
int Ifpack_DynamicFactory::NumPreconditioners_ = 0;

////==============================================================================
//const char* Ifpack_DynamicFactory::precTypeNames[Ifpack::numPrecTypes] =
//{
//#ifdef HAVE_IFPACK_SPARSKIT
//  ,"SPARSKIT"
//#endif
//#ifdef HAVE_IFPACK_HIPS
//  ,"HIPS"
//#endif
//#ifdef HAVE_HYPRE
//  ,"Hypre"
//#endif
//#ifdef HAVE_IFPACK_SUPERLU
//  ,"SILU"
//#endif
//  ,"Chebyshev"
//  ,"IHSS"
//  ,"SORa"
//};

bool Ifpack_DynamicFactory::Initialize()
{
	if (! Initialized_) {
		PreconditionerMap_["point relaxation"] = &buildPointRelaxation;
		PreconditionerMap_["point relaxation stand-alone"]
			  = &buildPointRelaxationStandAlone;
		PreconditionerMap_["block relaxation"] = &buildBlockRelaxation;
		PreconditionerMap_["block relaxation stand-alone"]
			  = &buildBlockRelaxationStandAlone;
		PreconditionerMap_["block relaxation stand-alone (ILU)"]
			  = &buildBlockRelaxationStandAloneILU;

#ifdef HAVE_IFPACK_AMESOS
		PreconditionerMap_["block relaxation stand-alone (Amesos)"]
			  = &buildBlockRelaxationStandAloneAmesos;
		PreconditionerMap_["block relaxation (Amesos)"]
			  = &buildBlockRelaxationAmesos;
		PreconditionerMap_["Amesos"]
			  = &buildAmesos;
		PreconditionerMap_["Amesos stand-alone"]
			  = &buildAmesosStandAlone;
#endif // HAVE_IFPACK_AMESOS

		PreconditionerMap_["IC"] = &buildIC;
		PreconditionerMap_["IC stand-alone"] = &buildICStandAlone;
		PreconditionerMap_["ICT"] = &buildICT;
		PreconditionerMap_["ICT stand-alone"] = &buildICTStandAlone;
		PreconditionerMap_["ILU"] = &buildILU;
		PreconditionerMap_["ILU stand-alone"] = &buildILUStandAlone;
		PreconditionerMap_["ILUT"] = &buildILUT;
		PreconditionerMap_["ILUT stand-alone"] = &buildILUTStandAlone;

#ifdef HAVE_IFPACK_SPARSKIT
		PreconditionerMap_["SPARSKIT"]
			  = &buildSPARSKIT;
#endif

#ifdef HAVE_IFPACK_HIPS
		PreconditionerMap_["HIPS"]
			  = &buildHIPS;
#endif

#ifdef HAVE_HYPRE
		PreconditionerMap_["Hypre"]
			  = &buildHypre;
#endif

#ifdef HAVE_IFPACK_SUPERLU
		PreconditionerMap_["SILU"]
			  = &buildSILU;
#endif

		PreconditionerMap_["Chebyshev"] = &buildChebyshev;

#ifdef HAVE_IFPACK_EPETRAEXT
		PreconditionerMap_["IHSS"] = &buildIHSS;
		PreconditionerMap_["SORa"] = &buildSORa;
#endif

		NumPreconditioners_ =
			    +5
		#ifdef HAVE_IFPACK_AMESOS
			    +4
		#endif
			    +8
		#ifdef HAVE_IFPACK_SPARSKIT
			    +1
		#endif
		#ifdef HAVE_IFPACK_HIPS
			    +1
		#endif
		#ifdef HAVE_HYPRE
			    +1
		#endif
		#ifdef HAVE_IFPACK_SUPERLU
			    +1
		#endif
			    +1
		#ifdef HAVE_IFPACK_EPETRAEXT
			    +2
		#endif
			    ;

		Initialized_ = true;
	}

	return true;
}

int Ifpack_DynamicFactory::RegisterPreconditioner(
		const std::string PrecName,
	    Ifpack_DynamicFactory::builderFunction PrecBuilder)
{
	if (PreconditionerMap_.find(PrecName) == PreconditionerMap_.end()) {
		PreconditionerMap_[PrecName] = PrecBuilder;
		NumPreconditioners_++;
		return 0;
	}
	return 1;
}

void Ifpack_DynamicFactory::Print(std::ostream& os)
{
	os << "Ifpack_DynamicFactory registered preconditioners: " << std::endl;
	for (std::map<std::string, builderFunction>::const_iterator it = PreconditionerMap_.begin();
		 it != PreconditionerMap_.end(); ++it) {
		os << it->first << std::endl;
	}
}

//==============================================================================
Ifpack_Preconditioner* Ifpack_DynamicFactory::Create(const string PrecType,
                                      Epetra_RowMatrix* Matrix,
                                      const int Overlap,
                                      bool overrideSerialDefault)
{
  bool serial = (Matrix->Comm().NumProc() == 1);

  std::map<std::string, builderFunction>::const_iterator it
  		= PreconditionerMap_.find(PrecType);
  bool found = it != PreconditionerMap_.end();
  if (found) {
	builderFunction f = it->second;
	return f(Matrix, Overlap, serial, overrideSerialDefault);
  }

  return 0;
}

// Let's initialize the factory upon compilation
namespace {
  bool init = Ifpack_DynamicFactory::Initialize();
}

