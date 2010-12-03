/*
// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
//  
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
// 
// ***********************************************************************
// 
// @HEADER

*/

// Teko includes
#include "Teko_DiagnosticPreconditionerFactory.hpp"

#include "Teko_PreconditionerInverseFactory.hpp"
#include "Teko_DiagnosticLinearOp.hpp"

#include "Teuchos_TimeMonitor.hpp"

namespace Teko {

//! Default constructor, for use with the AutoClone class.
DiagnosticPreconditionerFactory::DiagnosticPreconditionerFactory()
   : outputStream_(Teko::getOutputStream()), invFactory_(Teuchos::null), diagString_("<label me!>"), printResidual_(false)
{ }

/** Construct a preconditioner factory that applies a specified
  * preconditioner, a fixed number of times.
  */
DiagnosticPreconditionerFactory::DiagnosticPreconditionerFactory(const Teuchos::RCP<Teko::InverseFactory> & invFactory, const std::string & label,
                                                                 const Teuchos::RCP<std::ostream> & os,bool printResidual)
   : outputStream_(Teko::getOutputStream()), invFactory_(invFactory), diagString_(label), printResidual_(printResidual)
{ 
   initTimers(diagString_);

   if(os!=Teuchos::null) 
      outputStream_ = os;
}

DiagnosticPreconditionerFactory::~DiagnosticPreconditionerFactory()
{ 
   // check timers for null
   if(buildTimer_==Teuchos::null || rebuildTimer_==Teuchos::null) {
      // (*outputStream_) << "DiagnosticPreconditionerFactory \"" << diagString_ << "\": "
      //                  << "Timers not initialized" << std::endl;

      return;
   }

   double initBuildTime = totalInitialBuildTime();
   int initBuilds = numInitialBuilds();

   double initRebuildTime = totalRebuildTime();
   int initRebuilds = numRebuilds();

   (*outputStream_) << "DiagnosticPreconditionerFactory \"" << diagString_ << "\":\n";

   // print build string
   (*outputStream_) << "   build elapsed = " << initBuildTime << ", "
                    << "num builds = " << initBuilds << ", ";
   if(initBuilds>0)
      (*outputStream_) << "timer/app = " << initBuildTime / double(initBuilds) << "\n";
   else
      (*outputStream_) << "timer/app = " << "none" << "\n";

   // print rebuild string
   (*outputStream_) << "   rebuild elapsed = " << initRebuildTime << ", "
                    << "num rebuilds = " << initRebuilds << ", ";
   if(initRebuilds>0)
      (*outputStream_) << "timer/app = " << initRebuildTime / double(initRebuilds) << "\n";
   else
      (*outputStream_) << "timer/app = " << "none" << "\n";

   // print total string
   (*outputStream_) << "   total elapsed = " << initRebuildTime+initBuildTime << ", "
                    << "num rebuilds = " << initRebuilds+initBuilds << ", ";
   if(initBuilds+initRebuilds>0)
      (*outputStream_) << "timer/app = " << (initRebuildTime+initBuildTime) / double(initRebuilds+initBuilds) << std::endl;
   else
      (*outputStream_) << "timer/app = " << "none" << std::endl;
}

/** \brief Function that is called to build the preconditioner
  *        for the linear operator that is passed in.
  */
LinearOp DiagnosticPreconditionerFactory::buildPreconditionerOperator(LinearOp & lo,PreconditionerState & state) const
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;

   TEST_FOR_EXCEPTION(invFactory_==Teuchos::null,std::runtime_error,
                      "ERROR: Teko::DiagnosticPreconditionerFactory::buildPreconditionerOperator requires that an "
                   << "inverse factory has been set. Currently it is null!");

   TEST_FOR_EXCEPTION(buildTimer_==Teuchos::null || rebuildTimer_==Teuchos::null,std::runtime_error,
                      "ERROR: Teko::DiagnosticPreconditionerFactory::buildPreconditionerOperator requires that "
                   << "the timers be initialized. Currently they are null! (label = \"" << diagString_ << "\")");

   // build user specified preconditioner
   ModifiableLinearOp & diagOp_ptr = state.getModifiableOp("diagnosticOp");

   if(diagOp_ptr==Teuchos::null) {
      ModifiableLinearOp invOp;
      {
         // start timer on construction, end on destruction
         Teuchos::TimeMonitor monitor(*buildTimer_,false);   

         invOp = Teko::buildInverse(*invFactory_,lo);
      }

      // only printing residual requires use of forward operator
      if(printResidual_)
         diagOp_ptr = createDiagnosticLinearOp(outputStream_,lo,invOp,diagString_);
      else
         diagOp_ptr = createDiagnosticLinearOp(outputStream_,invOp,diagString_);
   }
   else {
      RCP<DiagnosticLinearOp> diagOp = rcp_dynamic_cast<DiagnosticLinearOp>(diagOp_ptr);

      // only printing residual requires use of forward operator
      if(printResidual_) 
         diagOp->setForwardOp(lo);

      ModifiableLinearOp invOp = diagOp->getModifiableOp();
      {
         // start timer on construction, end on destruction
         Teuchos::TimeMonitor monitor(*rebuildTimer_,false);  

         Teko::rebuildInverse(*invFactory_,lo,invOp);
      }
   }

   return diagOp_ptr.getConst();
}

/** \brief This function builds the internals of the preconditioner factory
  *        from a parameter list.
  */
void DiagnosticPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList & settings)
{
   TEST_FOR_EXCEPTION(not settings.isParameter("Inverse Factory"),std::runtime_error,
                      "Parameter \"Inverse Factory\" is required by a Teko::DiagnosticPreconditionerFactory");
   TEST_FOR_EXCEPTION(not settings.isParameter("Descriptive Label"),std::runtime_error,
                      "Parameter \"Descriptive Label\" is required by a Teko::DiagnosticPreconditionerFactory");
      
   // grab library and preconditioner name
   std::string invName = settings.get<std::string>("Inverse Factory");
   diagString_ = settings.get<std::string>("Descriptive Label");

   // build preconditioner factory
   Teuchos::RCP<const InverseLibrary> il = getInverseLibrary();
   invFactory_ = il->getInverseFactory(invName);
   TEST_FOR_EXCEPTION(invFactory_==Teuchos::null,std::runtime_error,
                      "ERROR: \"Inverse Factory\" = " << invName
                   << " could not be found");

   if(settings.isParameter("Print Residual"))
      printResidual_ = settings.get<bool>("Print Residual");

   // build timers to use
   initTimers(diagString_);
}

/** \brief Request the additional parameters this preconditioner factory
  *        needs. 
  */
Teuchos::RCP<Teuchos::ParameterList> DiagnosticPreconditionerFactory::getRequestedParameters() const
{
   TEST_FOR_EXCEPTION(invFactory_==Teuchos::null,std::runtime_error,
                      "ERROR: Teko::DiagnosticPreconditionerFactory::getRequestedParameters requires that a "
                   << "preconditioner factory has been set. Currently it is null!");

   return invFactory_->getRequestedParameters();
}

/** \brief Update this object with the fields from a parameter list.
  */
bool DiagnosticPreconditionerFactory::updateRequestedParameters(const Teuchos::ParameterList & pl)
{
   TEST_FOR_EXCEPTION(invFactory_==Teuchos::null,std::runtime_error,
                      "ERROR: Teko::DiagnosticPreconditionerFactory::updateRequestedParameters requires that a "
                   << "preconditioner factory has been set. Currently it is null!");

   return invFactory_->updateRequestedParameters(pl);
}

void DiagnosticPreconditionerFactory::initTimers(const std::string & str)
{
   buildTimer_ = Teuchos::rcp(new Teuchos::Time(str+" buildTimer")); 
   rebuildTimer_ = Teuchos::rcp(new Teuchos::Time(str+" rebuildTimer")); 
}

} // end namespace Teko
