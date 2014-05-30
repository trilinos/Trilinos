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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
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

std::map<std::string, Ifpack_DynamicFactory::builderFunction>
	Ifpack_DynamicFactory::PreconditionerMap_;
bool Ifpack_DynamicFactory::Initialized_ = false;
int Ifpack_DynamicFactory::NumPreconditioners_ = 0;

bool Ifpack_DynamicFactory::Initialize()
{
	if (! Initialized_) {
		PreconditionerMap_["point relaxation"]
			  = &buildPreconditioner<Ifpack_PointRelaxation, false>;
		PreconditionerMap_["point relaxation stand-alone"]
			  = &buildPreconditioner<Ifpack_PointRelaxation, true>;
		PreconditionerMap_["block relaxation"]
              = &buildPreconditioner<Ifpack_BlockRelaxation<Ifpack_DenseContainer>, false>;
		PreconditionerMap_["block relaxation stand-alone"]
			  = &buildPreconditioner<Ifpack_BlockRelaxation<Ifpack_DenseContainer>, true>;
		PreconditionerMap_["block relaxation stand-alone (ILU)"]
              = &buildPreconditioner<Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_ILU> >, true>;

#ifdef HAVE_IFPACK_AMESOS
		PreconditionerMap_["block relaxation stand-alone (Amesos)"]
              = &buildPreconditioner<Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_Amesos> >, true>;
		PreconditionerMap_["block relaxation (Amesos)"]
              = &buildPreconditioner<Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_Amesos> >, false>;
		PreconditionerMap_["Amesos"]
              = &buildPreconditioner<Ifpack_Amesos, false>;
		PreconditionerMap_["Amesos stand-alone"]
              = &buildPreconditioner<Ifpack_Amesos, true>;
#endif // HAVE_IFPACK_AMESOS

		PreconditionerMap_["IC"] = &buildPreconditioner<Ifpack_IC, false>;
		PreconditionerMap_["IC stand-alone"] = &buildPreconditioner<Ifpack_IC, true>;
		PreconditionerMap_["ICT"] = &buildPreconditioner<Ifpack_ICT, false>;
		PreconditionerMap_["ICT stand-alone"] = &buildPreconditioner<Ifpack_ICT, true>;
		PreconditionerMap_["ILU"] = &buildPreconditioner<Ifpack_ILU, false>;
		PreconditionerMap_["ILU stand-alone"] = &buildPreconditioner<Ifpack_ILU, true>;
		PreconditionerMap_["ILUT"] = &buildPreconditioner<Ifpack_ILUT, false>;
		PreconditionerMap_["ILUT stand-alone"] = &buildPreconditioner<Ifpack_ILUT, true>;

#ifdef HAVE_IFPACK_SPARSKIT
		PreconditionerMap_["SPARSKIT"]
			  = &buildPreconditioner<Ifpack_SPARSKIT, true>;
#endif

#ifdef HAVE_IFPACK_HIPS
		PreconditionerMap_["HIPS"]
			  = &buildPreconditioner<Ifpack_HIPS, true>;
#endif

#ifdef HAVE_HYPRE
		PreconditionerMap_["Hypre"]
			  = &buildPreconditioner<Ifpack_Hypre, true>;
#endif

#ifdef HAVE_IFPACK_SUPERLU
		PreconditionerMap_["SILU"]
			  = &buildPreconditioner<Ifpack_SILU, true>;
#endif

		PreconditionerMap_["Chebyshev"]
              = &buildPreconditioner<Ifpack_Chebyshev, true>;

#ifdef HAVE_IFPACK_EPETRAEXT
		PreconditionerMap_["IHSS"]
              = &buildPreconditioner<Ifpack_IHSS, true>;
		PreconditionerMap_["SORa"]
              = &buildPreconditioner<Ifpack_SORa, true>;
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

