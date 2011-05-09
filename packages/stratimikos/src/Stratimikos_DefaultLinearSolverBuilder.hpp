// @HEADER
// ***********************************************************************
// 
//         Stratimikos: Thyra-based strategies for linear solvers
//                Copyright (2006) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef STRATIMIKOS_DEFAULT_LINEAR_SOLVER_BUILDING_BASE
#define STRATIMIKOS_DEFAULT_LINEAR_SOLVER_BUILDING_BASE

#include "Stratimikos_ConfigDefs.hpp"
#include "Thyra_LinearSolverBuilderBase.hpp"
#include "Teuchos_AbstractFactory.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

// Include these to make all of the helpful decls appear
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_VectorBase.hpp"


namespace Teuchos { class CommandLineProcessor; }


namespace Stratimikos {


/** \brief . */
using Teuchos::RCP;
/** \brief . */
using Teuchos::Array;
/** \brief . */
using Teuchos::AbstractFactory;
/** \brief . */
using Teuchos::ParameterList;


/** \brief Concrete subclass of <tt>Thyra::LinearSolverBuilderBase</tt> for
 * creating <tt>LinearOpWithSolveFactoryBase</tt> objects and
 * <tt>PreconditionerFactoryBase</tt> object on demand for various Trilinos
 * linear solver packages.
 *
 * For an example of how to use this class see
 *  <a href="simple__stratimikos__example_8cpp-example.html">simple_stratimikos_example.cpp</a></tt>.
 *
 * The parameters this class accepts are shown below in different format:
 * <ul>
 * <li> \ref HumanReadableWithDocumentation "Human readable format (with documentation) for valid parameters accepted by this class"
 * <li> \ref HumanReadableWithoutDocumentation "Human readable format (without documentation) for valid parameters accepted by this class"
 * <li> \ref XmlFormat "XML format for valid parameters accepted by this class"
 * </ul>
 *
 * <b>\anchor HumanReadableWithDocumentation Human readable format (with documentation) for valid parameters accepted by this class</b>
 *
 * \verbinclude simple_stratimikos_example.options.readablewithdoc.out
 *
 * <b>\anchor HumanReadableWithoutDocumentation Human readable format (without documentation) for valid parameters accepted by this class</b>
 *
 * \verbinclude simple_stratimikos_example.options.readable.out
 *
 * <b>\anchor XmlFormat XML format for valid parameters accepted by this class</b>
 *
 * \verbinclude simple_stratimikos_example.options.xml.out
 * 
 */
class DefaultLinearSolverBuilder
  : public Thyra::LinearSolverBuilderBase<double>
{
public:

  /** @name Constructors/Initializers/Accessors */
  //@{

  /** \brief Construct with default parameters.
   *
   * <b>Warning!</b> Do not change the defaults by passing then into this
   * constructor.  Instead, use the member functions to set them after
   * <tt>*this</tt> is constructed.  This will help to avoid problems with
   * updates to the ordering of the arguments.
   */
  DefaultLinearSolverBuilder(
    const std::string    &paramsXmlFileName                = ""
    ,const std::string   &extraParamsXmlString             = ""
    ,const std::string   &paramsUsedXmlOutFileName         = ""
    ,const std::string   &paramsXmlFileNameOption          = "linear-solver-params-file"
    ,const std::string   &extraParamsXmlStringOption       = "extra-linear-solver-params"
    ,const std::string   &paramsUsedXmlOutFileNameOption   = "linear-solver-params-used-file"
    );

  /** \brief . */
  ~DefaultLinearSolverBuilder();
  
  /** \brief The name an XML file that will be read to get XML parameters (if
   * not "").
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS(std::string,paramsXmlFileName);
    
  /** \brief An XML string that will be used to update the parameters (if not
   * "").
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS(std::string,extraParamsXmlString);

  /** \brief The name of an XML file that will be written (if not "") for the
   * parameters actually used.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS(std::string,paramsUsedXmlOutFileName);

  /** \brief The name of the option that will be added the the commandline
   * processor that will set <tt>paramsXmlFileName()</tt> .
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS(std::string,paramsXmlFileNameOption);

  /** \brief The name of the option that will be added the the commandline
   * processor that will set <tt>extraParamsXmlString()</tt> .
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS(std::string,extraParamsXmlStringOption);

  /** \brief The name of the option that will be added the the commandline
   * processor that will set <tt>paramsUsedXmlOutFileName()</tt> .
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS(std::string,paramsUsedXmlOutFileNameOption);

  /** \brief Set a new linear solver strategy factory object. */
  void setLinearSolveStrategyFactory(
    const RCP<const AbstractFactory<Thyra::LinearOpWithSolveFactoryBase<double> > >
    &solveStrategyFactory,
    const std::string &solveStrategyName,
    const bool makeDefault = false
    );

  /** \brief Set the default linear solver factory name. */
  void setDefaultLinearSolveStrategyFactoryName(
    const std::string &solveStrategyName);

  /** \brief Set a new preconditioner strategy factory object. */
  void setPreconditioningStrategyFactory(
    const RCP<const AbstractFactory<Thyra::PreconditionerFactoryBase<double> > >
    &precStrategyFactory,
    const std::string &precStrategyName,
    const bool makeDefault = false
    );

  /** \brief Set the default linear solver factory name. */
  void setDefaultPreconditioningStrategyFactoryName(
    const std::string &precStrategyName);

  /** \brief Setup the command-line processor to read in the needed data to
   * extra the parameters from.
   *
   * Command-line options with names <tt>this->paramsXmlFileNameOption()</tt>,
   * <tt>this->extraParamsXmlStringOption()</tt>, and
   * <tt>this->paramsUsedXmlOutFileNameOption()</tt> will be set if they are
   * not empty.
   *
   * Then, when <tt>cpl->parse(...)</tt> is called, then the options set will
   * be read into <tt>this->paramsXmlFileName()</tt>,
   * <tt>this->extraParamsXmlString()</tt>, and
   * <tt>this->paramsUsedXmlOutFileName()</tt>.
   *
   * After this function is called, <tt>this->readParameters()</tt> can be
   * called to actually read in the parameters and fill the parameter list.
   */
  void setupCLP( Teuchos::CommandLineProcessor *clp );

  /** \brief Force the parameters to be read from a file and/or an extra XML
   * string.
   *
   * First, if <tt>this->getParameterList().get()==NULL</tt> and new parameter
   * list will be created.
   *
   * Second, if <tt>this->paramsXmlFileName()!=""</tt> then the file
   * <tt>this->paramsXmlFileName()</tt> will be read to get XML parameters
   * append/update those already in the parameter list.
   *
   * Third, if <tt>this->extraParamsXmlString()!=""</tt> then the XML string
   * <tt>this->extraParamsXmlString()</tt> will be read and used to
   * append/update the parameters already in the parameter list..
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->getParameterList().get()!=NULL</tt>
   * </ul>
   */
  void readParameters( std::ostream *out );

  /** \brief Write the parameters list for a
   * <tt>LinearOpWithSolveFactoryBase</tt> object to a file after the
   * parameters are read in order to show defaults and create a new list for
   * input the next time.
   *
   * If <tt>outputXmlFileName!=""</tt> then the parameter list with be written
   * to the file <tt>outputXmlFileName</tt> in XML format. If
   * <tt>outputXmlFileName==""</tt>, but
   * <tt>this->paramsUsedXmlOutFileNameOption()!=""</tt> then the parameter
   * list will be written to the file
   * <tt>this->paramsUsedXmlOutFileNameOption()</tt>.  If both
   * <tt>outputXmlFileName==""</tt> and
   * <tt>this->paramsUsedXmlOutFileNameOption()==""</tt> then no file is
   * written.
   */
  void writeParamsFile(
    const Thyra::LinearOpWithSolveFactoryBase<double> &lowsFactory,
    const std::string &outputXmlFileName  = "" 
    ) const;
  
  /** \brief Get the name of the linear solver strategy that will be created
   * on the next call to <tt>this->createLinearSolverStrategy()</tt>.
   */
  std::string getLinearSolveStrategyName() const;

  /** \brief Get the name of the preconditioner strategy that will be created
   * on the next call to <tt>this->createPreconditioningStrategy()</tt>.
   */
  std::string getPreconditionerStrategyName() const;

  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);
  /** \brief . */
  RCP<ParameterList> getNonconstParameterList();
  /** \brief . */
  RCP<ParameterList> unsetParameterList();
  /** \brief . */
  RCP<const ParameterList> getParameterList() const;
  /** \brief . */
  RCP<const ParameterList> getValidParameters() const;

  //@}
  
  /** \name Overridden from LinearSolverBuilderBase. */
  //@{

  /** \brief . */
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
  createLinearSolveStrategy(
    const std::string &linearSolveStrategyName
    ) const;
  /** \brief . */
  RCP<Thyra::PreconditionerFactoryBase<double> >
  createPreconditioningStrategy(
    const std::string &preconditioningStrategyName
    ) const;

  //@}

private:

  // //////////////////////////////////////
  // Private types

  typedef RCP<const AbstractFactory<Thyra::LinearOpWithSolveFactoryBase<double> > >
  lowsf_fcty_t;
  typedef RCP<const AbstractFactory<Thyra::PreconditionerFactoryBase<double> > >
  pf_fcty_t;

  // //////////////////////////////////////
  // Private data members
  
  RCP<ParameterList> paramList_;
  Array<std::string> validLowsfNames_;
  Array<lowsf_fcty_t> lowsfArray_;
  std::string defaultLOWSF_;
  Array<std::string> validPfNames_; // Contains "None" as the 0th entry!
  Array<pf_fcty_t> pfArray_;
  std::string defaultPF_;
  bool enableDelayedSolverConstruction_;
  mutable RCP<const ParameterList> validParamList_;
  mutable RCP<const Teuchos::StringToIntegralParameterEntryValidator<int> > lowsfValidator_;
  mutable RCP<const Teuchos::StringToIntegralParameterEntryValidator<int> > pfValidator_;

  // //////////////////////////////////////
  // Private member functions

  void initializeDefaults();
  void justInTimeInitialize() const;

};


} // namespace Stratimikos


#endif // STRATIMIKOS_DEFAULT_LINEAR_SOLVER_BUILDING_BASE
