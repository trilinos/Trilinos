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

#include "Ifpack_PreconditionerFactory.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Ifpack_CrsRiluk.h"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Time.hpp"

// RAB: 2004/01/19: It turns out that as of this writting the Ifpack
// preconditioner class Ifpack_CrsRiluk does not allow reuse with
// InitValues(...) and Factor(...)  which is unfortunate!
#define IFPACK_PRECGENERATOR_CREATE_PREC_EVERY_TIME  1

namespace {
const std::string Op_iluk_graph_name = "Op_iluk_graph";
}

namespace Ifpack {

PreconditionerFactory::PreconditionerFactory(
  const int        levelFill
  ,const int       levelOverlap
  ,const double    absThreshold
  ,const double    relThreshold
	,const bool      calcCondEst
  )
  :levelFill_(levelFill)
  ,levelOverlap_(levelOverlap)
  ,absThreshold_(absThreshold)
  ,relThreshold_(relThreshold)
	,calcCondEst_(calcCondEst)
{}

void PreconditionerFactory::setupPrec(
  const Epetra_Operator                         &Op
  ,Teuchos::RefCountPtr<Epetra_Operator>        *Prec_in
	,std::ostream                                 *out
  ,const std::string                            &li
  ,const std::string                            &is
  ) const
{
  using std::endl;
  using Teuchos::dyn_cast;
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  using Teuchos::set_extra_data;
  using Teuchos::get_extra_data;

	Teuchos::Time timer("");

	if(out)
		*out <<endl<<li<<"*** Entering Ifpack::PreconditionerFactory::setupPrec(...) ...\n";

  RefCountPtr<Epetra_Operator> &Prec = *Prec_in;

	if(out)
		*out <<endl<<li<<"Trying to create preconditioner for Epetra_Operator of type \'"
				 << typeid(Op).name() << "\'.\n";
	
  // Determine what type of Epetra_Operator we have for DcDy
  const Epetra_RowMatrix *Op_rm = dynamic_cast<const Epetra_RowMatrix*>(&Op);
  if(Op_rm) {
		if(out)
			*out <<endl<<li<<"Detected that Op supports the Epetra_RowMatrix interface.\n";
    // Yea, We can use AztecOO or Ifpack preconditioners!  Now
    // determine if Op_rm is actually a a Crs or a Vbr matrix.
    const Epetra_CrsMatrix *Op_crs
      = dynamic_cast<const Epetra_CrsMatrix*>(Op_rm);
    const Epetra_VbrMatrix *Op_vbr
      = ( Op_crs ? (const Epetra_VbrMatrix*)NULL : dynamic_cast<const Epetra_VbrMatrix*>(Op_rm) );
    if( Op_crs || Op_vbr ) { 
			if(out)
				*out <<endl<<li<<"Detected that Op supports the "
						 << ( Op_crs ? "Epetra_CrsMatrix" : "Epetra_VbrMatrix"  )
						 << " interface.\n";
      //
      // Op is a Crs or a Vbr matrix!
      //
      // Create the preconditioner if it has not been created already.
			const bool createPrec
#ifdef IFPACK_PRECGENERATOR_CREATE_PREC_EVERY_TIME
				= true
#else
				= !Prec.get()
#endif
				;
			if( createPrec ) {
        RefCountPtr<Ifpack_IlukGraph>  Op_iluk_graph;
#ifdef IFPACK_PRECGENERATOR_CREATE_PREC_EVERY_TIME
				if(Prec.get()) {
					// Get the precreated graph
					if(out)
						*out <<endl<<li<<"Grabing a preexisting Ifpack_IlukGraph object out of preconditioner ...\n";
					Op_iluk_graph = get_extra_data<RefCountPtr<Ifpack_IlukGraph> >(Prec,Op_iluk_graph_name);
				}
				else {
#endif
					// Create the graph
					if(out)
						*out <<endl<<li<<"Creating a new Ifpack_IlukGraph object with:"
								 <<endl<<li<<"  levelFill = " << levelFill()
								 <<endl<<li<<"  levelOverlap = " << levelOverlap()
								 << std::endl;
					timer.start(true);
					Op_iluk_graph = rcp(
            new Ifpack_IlukGraph(
              ( Op_crs
                ? Op_crs->Graph()
                : Op_vbr->Graph() )
              ,levelFill()
              ,levelOverlap()
              )
            );
					Op_iluk_graph->ConstructFilledGraph();
					timer.stop();
					if(out)
						*out <<endl<<li<<"  => time = " << timer.totalElapsedTime() << " sec\n";
#ifdef IFPACK_PRECGENERATOR_CREATE_PREC_EVERY_TIME
				}
#endif
				// Create the preconditioner
				if(out)
					*out <<endl<<li<<"Creating a new Ifpack_CrsRiluk object ...\n";
        Prec = rcp(new Ifpack_CrsRiluk(*Op_iluk_graph));
        // Give ownership of the graph to the preconditioner.  Note,
        // if multiple operator and multiple preconditioner objects
        // that have the same structure are to be used, this
        // implementation will create a new graph for each of these
        // preconditioner objects.  However, as long as these
        // preconditioner objects are reused over and over again then
        // this approach should not be a performance problem and this
        // greatly simplifies how this class is used.
        set_extra_data( Op_iluk_graph, Op_iluk_graph_name, &Prec );
			}
      // Get a Ifpack_CrsRiluk subclass pointer for prec
      Ifpack_CrsRiluk *Prec_crs_riluk = &dyn_cast<Ifpack_CrsRiluk>(*Prec);
      // Set diagonal perturbations.
			// Note: You must do this *before* you set the matrix values
			// (this is done wrong in NOX).
			if(out)
				*out <<endl<<li<<"Setting Ifpack_CrsRiluk object options:"
						 <<endl<<li<<"  absThreshold = " << absThreshold()
						 <<endl<<li<<"  relThreshold = " << relThreshold()
						 << std::endl;
      Prec_crs_riluk->SetAbsoluteThreshold(absThreshold());
      Prec_crs_riluk->SetRelativeThreshold(relThreshold());
      // Now initialize the values
			if(out)
				*out <<endl<<li<<"Initializing matrix values of the preconditioner ...\n";
			timer.start(true);
      if(Op_crs) {
				TEST_FOR_EXCEPTION(0!=Prec_crs_riluk->InitValues(*Op_crs),std::runtime_error,"Error!");
			}
      else if(Op_vbr) {
        TEST_FOR_EXCEPTION(0!=Prec_crs_riluk->InitValues(*Op_vbr),std::runtime_error,"Error!");
			}
      else
        assert(0); // Should never get here!
			timer.stop();
			if(out)
				*out <<endl<<li<<"  => time = " << timer.totalElapsedTime() << " sec\n";
      // Complete the factorization
			if(out)
				*out <<endl<<li<<"Factoring the preconditoner matrix ...\n";
			timer.start(true);
			TEST_FOR_EXCEPTION(0!=Prec_crs_riluk->Factor(),std::runtime_error,"Error!");
			timer.stop();
			if(out)
				*out <<endl<<li<<"  => time = " << timer.totalElapsedTime() << " sec\n";
			// Compute condition number estimate
			if(out && calcCondEst()) {
				double fwdCondEst = 0.0, adjCondEst = 0.0;
				TEST_FOR_EXCEPTION(0!=Prec_crs_riluk->Condest(false,fwdCondEst),std::runtime_error,"Error!");
				TEST_FOR_EXCEPTION(0!=Prec_crs_riluk->Condest(true,adjCondEst),std::runtime_error,"Error!");
				if(out)
					*out <<endl<<li<<"Computing condition number estimate (calcCondEst()==true):"
							 <<endl<<li<<"  Prec_crs_riluk->Condest(false,fwdCondEst): fwdCondEst = " << fwdCondEst
							 <<endl<<li<<"  Prec_crs_riluk->Condest(true,adjCondEst) : adjCondEst = " << adjCondEst
							 << std::endl;
			}
    }
    else {
      // It turns out that Op only supports the Epetra_RowMatrix interface!
			const char errMsg[]
				= "PreconditionerFactory::setupPrec(...): Have not implemented support "
				"for only Epetra_RowMatrix interface yet!";
			if(out)	*out << std::endl << errMsg << std::endl;
      TEST_FOR_EXCEPTION(true,std::logic_error,errMsg);
    }
  }
  else {
		const char errMsg[]
			= "PreconditionerFactory::setupPrec(...):Error, can't handle operator "
			"that does not support Epetra_RowMatrix!";
		if(out)	*out << std::endl << errMsg << std::endl;
    TEST_FOR_EXCEPTION(true,std::logic_error,errMsg);
  }
	
	if(out)
		*out <<endl<<li<<"*** Leaving Ifpack::PreconditionerFactory::setupPrec(...) ...\n";
	
}

} // namespace Ifpack
