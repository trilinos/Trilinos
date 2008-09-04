// @HEADER
// ************************************************************************
// 
//            Phalanx: A Partial Differential Equation Assembly 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#ifndef PHX_DOXYGEN_DOCUMENTATION_HPP
#define PHX_DOXYGEN_DOCUMENTATION_HPP

/*!

\mainpage

\section main_index Index

- \ref overview
- \ref user_guide
- \ref faq
- \ref bugs
- \ref history
- \ref authors
- \ref questions

\section overview Overview

Phalanx is an assembly kernel that provides fast and flexible
evaluations of field values with complex dependency chains for partial
differential equation systems.  It provides a flexible mechanism for
automatically switching field evaluation routines at run-time (and
thus the dependency chain).  A unique aspect of this library is that
it was written to allow for arbitrary data field types using template
metaprogramming concepts.  This opens the door for embedded technology
(embedded data types with overloaded mathematical operators) such as
arbitrary precision, sensitivity analysis (via embedded automatic
differentiation), and uncertainty quantification, all by swapping the
scalar type and reusing the same function evaluation code.  It can be
used with any cell-based discretization techniques including finite
element and finite volume.

Phalanx is a hammer.  It's use should be carefully considered.  It
should only be used for general PDE frameworks where the evaluation
routines (and thus dependencies) are changed at run-time without
recompiling the code.  Swapping both material models (such as stress
tensors in the Navier-Stokes equations) and material properties
between runs are two key examples.  Phalanx should be used for complex
dependency chains where those dependencies can change between runs.
It should not be used for simple sets of PDEs where the equations
rarely change.  Below we discuss advantages, disadvantages and some
use cases to guide whether Pahalnx is right for you.

\subsection advantages Advantages of using Phalanx

- Increased flexibility because each simpler piece becomes an
extension point that can be swapped out with different
implementations.

- Easier to implement new code because each piece is simpler, more
focused and easier to test in isolation.

- Ensures consistent evaluations of fields according to their
dependencies.  For example, if density were a field to be evaluated,
many models that have different dependencies.  Model A could depend on
temperature and pressure, whille model B could depend on temperature,
pressure, and species mass fractions of a multicomponent mixture.  The
order of evaluation of fields could change based on the change in
dependencies.  Using Model B will require that the species mass
fractions be evaluated and available for the density evaluation while
Model A does not even require them.  This is a simple example, but the
dependency chain can become quite complex when solving thousands of
coupled PDE equations.  Phalanx aill automatically perform the sorting
and ordering of the evaluator routines.

- Fast evaluation of cells.  Phalanx provides for an arbitrary
size block of cells to be evaluated on the local processor.  By using
the contiguous allocator and correctly sizing the number of cells to
fit in processor cache, one can arrange all fields to exist in a
contiguous block of memory, independent of the data type objects used
to store the data (one must still be careful of data alignement
issues).  By keeping all fields in cache, the code will run much
faster.

- Allows the reuse of the same code used to evaluate a function for
computing sensitivities and for uncertainty quantification.

- Automatically extended to arbitrary data types.

\subsection disadvantages Disadvantages of using Phalanx

- A potential performance loss due to fragmentation of the over-all
algorithm (e.g., several small loops instead of one large loop).

- A potential loss of visibility of the original, composite problem
 (since the code is scattered into multiple places).  

Managing these trade-offs can result in application code that both performs well and supports rapid development and extensibility.  

\subsection simple_examples Some Use Cases

This section outlines some more complex use cases where Phalanx should be used.

- 


\section bugs Reporting Bugs and Making Enhancement Requests

  To reports bugs or make enhancement requests, visit the <A HREF="http://software.sandia.gov/bugzilla/">Trilinos Bugzilla (Bug Tracking) Database</A>, and use the following instructions.
      <UL>
      <LI>Click on "Enter a new bug report"
      <LI>Choose "Phalanx"
      <LI>Either login or create a new account
      <LI>Submit your bug report
      </UL>

\section history History

Phalanx grew out of the Variable Manger in the Charon code and the Expression Manager in the Aria code.  It is an attempt at merging the two capabilities and is slated to provide nonlinear function evaluation to the Intrepid discretiztion package.

\section authors Authors

The following have contributed to the design and/or development of Phalanx:

  - Roger Pawlowski (PI), SNL 01414
  - Eric Phipps, SNL 01411
  - Pat Notz, SNL 01541 

\section questions For All Other Questions and Comments...
  
   Please contact Roger Pawlowski (rppawlo@sandia.gov).

*/

/* ************************************************************************ */
/* ************************************************************************ */

/*! \page user_guide Users Guide

\section user_guide_index Index

- \ref user_guide_getting_started
- \ref user_guide_domain_model
- \ref user_guide_mdarray_domain_model
- \ref user_guide_step1
- \ref user_guide_step2
- \ref user_guide_step3
- \ref user_guide_step4
- \ref user_guide_step5
- \ref user_guide_step6

\section user_guide_getting_started Getting Started

\subsection ug_dummy_1 Understand Templates
Phalanx is a complex package that make heavy use of the C++ templating mechanism.  We recommend that users first familiarize themselves with C++ templates.  An excellent reference is "C++ Templates: The Complete Guide" by Vandevoorde and Josuttis, 2003.  While users do not need to understand template metaprogramming to use Phalanx, the concepts underlying many operations in Phalanx can be found in "C++ Template Metaprogramming" by Abrahams and Gurtovoy, 2004.

\subsection ug_dummy_2 Learn the Phalanx Nomenclature
Users should then learn the nomeclature used in the package defined in the \ref user_guide_domain_model.

\subsection ug_dummy_3 Tutorial

The main concept of Phalanx is to evaluate fields for solving PDEs.  Let's start with a simple example found in Trilinos/packages/phalanx/example/EnergyFlux.  Suppose that we want to solve the heat equation over the physical space \f$ \Omega \f$:

\f[
  \nabla \cdot (-\rho C_p \nabla T) + s = 0
\f] 

 \f$s\f$ is a nonlinear source term.  The specific discretization technique whether finite element (FE) of finite volume (FV) will ask for \f$\mathbf{q}\f$ and \f$s\f$ at points on the cell.  Phalanx will evaluate \f$\mathbf{q}\f$ and \f$s\f$ at those points and return them to the discretization driver.

Using finite elements, we pose the problem in variational form:

Find \f$ u \in {\mathit{V^h}} \f$ and \f$ \phi \in {\mathit{S^h}} \f$ such that:

\f[
  - \int_{\Omega} \nabla \phi \cdot \mathbf{q} d\Omega 
  + \int_{\Gamma} \phi \mathbf{n} \cdot \mathbf{q} 
  + \int_{\Omega} \phi s d\Omega = 0 
\f]

\f[
  \mathbf{q} = -\rho C_p \nabla T
\f]

Follow the steps below to integrate Phalanx into your application.  The example code shown in the steps comes from the energy flux example in the directory "phalanx/example/EnergyFlux".  
- \ref user_guide_step1
- \ref user_guide_step2
- \ref user_guide_step3
- \ref user_guide_step4
- \ref user_guide_step5
- \ref user_guide_step6

\section user_guide_domain_model Phalanx Domain Model

<ul>
<li><h3>%Cell</h3>
Partial differential equations are solved in a domain.  This domain is discretized into cells (also called elements for the finite element method).  This library assumes that the block of cells being iterated over is of the same type!  If different evaluators (i.e. different material properties) are required in different blocks of cells, a new FieldMangager must be used for each unique block of elements.  This is required for efficiency.  

<li><h3>Scalar Type</h3>
A scalar type, typically the template argument ScalarT in Phalanx code, is the type of scalar used in an evaluation.  It is typically a double, but can be special object types for embedded methods such as sensitivity analysis.  for example, for sensitivity analysis, a double scalar type is replaced with a foward automatic differentiation object (FAD) or a reverse automatic differentaion object (RAD) to produce sensitivity information.  Whatever type is used, the standard mathematical operators are overloaded for the particular embedded technology.  For an example of this, see the <a href="http://trilinos.sandia.gov/packages/sacado">Sacado Automatic Differentiation Library</a>.

<li><h3>Algebraic Type</h3>
An algebraic type is the type of object that hold data.  It is usually a rank n tensor.  Simple examples include a scalar (rank-0 tensor), a vector (rank 1 tensor) or a matrix (rank 2 tensor).  It is not actually restircted to tensors, but can be any struct/class.

<li><h3>Data Type</h3>
A data type, typically the template argument DataT in Phalanx code, is an actual type used for storing fields.  It is the combination of a specific scalar type and an algebraic type.  For example it can be a vector object of doubles: MyVector<double>, or it could be a tensor object of FAD type: MyTensor<Scadao::FAD<double> >, or it could just be a double.

<li><h3>Evaluation Type</h3>
The evaluation type, typically the template argument EvalT in Phalanx code, defines a unique type of evaluation to perform.  An EvaluationContainer is allocated for each evaluation type specified in the traits class.  Examples include a residual type, a Jacobian type, and a parameter sensitivity type.  The EvaluationContainer is associated with at least one scalar type and possibly more.  The scalar type usually determines what is being evaluated. For example, to evaluate the equation residuals, the scalar type is usually a double.  To evaluate a Jacobian, the scalar type could be a forward automatic differentiation object, Sacado::FAD<double>.  By using an evaluation type, the same scalar type can be used for different evaluation types.  For example computing the Jacobian and computing parameter sensitivities both could use the Sacado::FAD scalar type.

<li><h3>Storage</h3>
A DataContainer object stores all fields of a particular data type.  Each EvaluationContainer holds a vector of DataContainers, one DataContainer for each vaid data type that is associated with that particular evaluation type.  One EvaluationContainer is constructed for each evaluation type.

<li><h3>Data Layout</h3>
The DataLayout object is used to distinguish fields with the same name, but exist at different discretization points in the cell.  For example, supposed we have written an evaluator the computes the "Density" field at a point.  Now we want to evaluate the density at multiple points in the domain.  We might have a "Density" field in a cell associated with a set of integration points (quadrature points in finite elements) and another field associated with the nodes (nodal basis points in finite elements).  We use the same field name (so we can reuse the same Evaluator), "Density", but use two different DataLayouts, one for integration points and one for nodes.  Now a FieldTag comparison will differentiate the fields due to the different DataLayout.  Additionally, the DataLayout contains the number of DataT objects associated with the field in a single cell.  This size() parameter is not needed to distinguish uniqueness, since the number of objects can be the same for different fields.  It is stored here for convenience when figuring out the size of field arrays.

<li><h3>Field Tag</h3>
The FieldTag is a description of a field.  It is templated on the data type, DataT.  It is used to describe every field that is stored in the field manager.  It contains a unique identifier (a stl std::string) and a pointer to a data layout object.  If two FieldTags are equal, the DataLayout, the data type, and the string name are exactly the same.

<li><h3>Field</h3>
A Field is a set of values for a particular quantity of interest that must be evaluated.  It is templated on the data type, DataT.  It consists of a FieldTag and a reference counted smart pointer (Teuchos::RCP<DataT>) to the data array where the values are stored.  

<li><h3>Evaluator</h3>
An Evaluator is an object that evaluates a set of Fields.  It contains two vectors of Fields, one set is the fields it will evaluate and one set is the fields it depends on for the evaluation.  For example to evaluate a density field that is a function of temperature and pressure, the field the evaluator will evaluate (the first field set) is a density field, and the set of fields it requires to meet its dependencies are the temperature and pressure fields (the second field set).  The evaluator is templated on the evaluation type and the traits class.  

<li><h3>Evaluator Manager</h3>
The main object that stores all fields and evauation routines. Users can get access to all fields through this object. 

</ul>

\section user_guide_mdarray_domain_model Multi-Dimensional Array Domain Model

Document is under development.  Will have reference shortly.

\section user_guide_step1 Step 1: Configuring, Building, and installing Phalanx

Phalanx is distributed as a package in the <a href="http://trilinos.sandia.gov">Trilinos Framework</a>.  It can be enabled as part of a trilinos build with the configure option "--enable-phalanx".  Phalanx currently has direct dependencies on the following third party libraries:

 - Requires the <a href="http://trilinos.sandia.gov/packages/teuchos">Teuchos</a> utilities library, part of the <a href="http://trilinos.sandia.gov/">Trilinos Framework</a>.
 
 - Requires the <a href="http://trilinos.sandia.gov/packages/sacado">Sacado Automatic Differentiation Library</a>, part of the <a href="http://trilinos.sandia.gov/">Trilinos Framework</a>.

 - Requires the <a href="http://www.boost.org">Boost Template Metaprogramming (MPL) Library</a>

 - Optional: Some performance tests run comparisons against <a href="http://tvmet.sourceforge.net/">TVMET: Tiny Vector Matrix library using Expression Templates</a>.  This is to get a feel for how our "dumb" vector matrix objects perform compared to expression templates.  Use the configure option "--with-tvmet" to enable the tvmet functionality in the performance tests.

We plan to remove the Boost dependencies in a future release.  

\subsection da A. Install Boost
You must have boost installed on your system.  As we only require the MPL library headers, you only need to untar the source code and should not even need to run the configure script.

\subsection db B. Configure Trilinos/Phalanx
The general instructions for building trilinos can be found at <a href="http://trilinos.sandia.gov/documentation.html">Trilinos Documentation Page</a>.  Of particular importance are the Overview, User Guide, and Tutorial documents.  At a minimum you must enable the Teuchos, Sacado, and Phalanx packages.  An example configure script is:

\verbinclude reconfigure

Once configure is run, build and install the library with the command:

\code
make install
\endcode

\section user_guide_step2 Step 2: Determine Types

Users must next determine the evaluation types and data types they will require in their simulation.  Please see the section \ref user_guide_domain_model for detailed explanation of types.  Following the example in the phalanx/example/EnergyFlux directory, we will be requiring two evaluation types, one for the residual evaluation of the discretized PDE equation and one for the corresponding Jacobian.  Additional evaluation types might be for the parameter sensitivites and uncertainty quantification.  

Once the evaluation types are chosen, users must decide on a default scalar type and on all data types that are valid for each evaluation type.  The data types should all be templated on the default scalar type for the particular evaluation type, but this is not a requirement (expert users can violate this for performance reasons).  Uses must implement their own data types or get them from a separate library.  For sensitivities, the trilinos package <a href="http://trilinos.sandia.gov/packages/sacado">Sacado</a> should be used.  For uncertainty quantification, the Trilinos package <a href="http://trilinos.sandia.gov/packages/stokhos">Stokhos</a> should be used.

In our example, the user has written an implementation of vector (MyVector) and matrix (MyTensor) classes found in the file AlgebraicTypes.hpp.  They are templated on a scalar type and can be used for all evaluation types. 

<ul>
<li> Residual - Default scalar type will be "double".  The data types will be:
<ul>
<li> double
<li> MyVector<double>
<li> MyTensor<double>
</ul>
<li> Jacobian - Default scalar type will be "Sacado::Fad::DFad<double>".  This scalar type will carry both the residual information as well as sensitivity information.  The data types will be:
<ul>
<li> Sacado::Fad::DFad<double>
<li> MyVector< Sacado::Fad::DFad<double> >
<li> MyTensor< Sacado::Fad::DFad<double> >
</ul>
</ul>

Typical examples of algebraic types include vectors, matrices, or higher order tensor objects, but in reality can be any struct/class that the user desires.  Remember to template the objects on the scalar type if possible.  An example of a very inefficient Vector and Matrix implementation (operator overloading without expression templates) can be found in AlgebraicTypes.hpp.

\section user_guide_step3 Step 3: Write the Traits Object

The Traits object is a struct that defines the evaluation types, the data types for each evaluation type, the allocator type, and the user defined types that will be passed through evaluator calls.  We will go through each of these defninitions.

The basic class should derive from the PHX::TraitsBase object.

The Traits struct must define each of the following typedef members of the struct:

- EvalTypes - an mpl::vector of user defined evaluation types.  Each evaluation type must have a typedef member called ScalarT that provides the default scalar type.  This is used to automate the building of evaluators for each evaluation type using the EvaluatorFactory.
      
- EvalToDataMap - an mpl::map.  The key is an evaluation type and the value is an mpl::vector of valid data types for that particular evaluation type.
      
- Allocator type - type that defines the allocator class to use to allocate the memory for data storage.
      
- EvalData - A user defined type to be passed in to the evaluateFields() call.  Allows users to pass in arbitrary data on the cells.

- PreEvalData - A user defined type to be passed in to the preEvaluate() call.  Allows users to pass in arbitrary data on the cells.

- EvalData - A user defined type to be passed in to the postEvaluate() call.  Allows users to pass in arbitrary data on the cells.

\subsection td1 A. Basic Class
the basic outline of the traits struct is:
\code
struct MyTraits : public PHX::TraitsBase {
    .
    .
    .
};
\endcode

Inside this struct we need to implement all the typedefs listed above.  The example we will follow is in the file Traits.hpp in "phalanx/example/EnergyFlux" directory.

\subsection td2 B. EvalTypes
First we need to know the evaluation types.  Each evaluation type must include at least a typedef'd default scalar type argument as a public member called ScalarType.  Here is the example code:

\code
  struct MyTraits : public PHX::TraitsBase {
    
    // ******************************************************************
    // *** Scalar Types
    // ******************************************************************
    
    // Scalar types we plan to use
    typedef double RealType;
    typedef Sacado::Fad::DFad<double> FadType;
    
    // ******************************************************************
    // *** Evaluation Types
    // ******************************************************************
    struct Residual { typedef RealType ScalarT; };
    struct Jacobian { typedef FadType ScalarT;  };
    typedef Sacado::mpl::vector<Residual, Jacobian> EvalTypes;

     .
     .
     .
\endcode 

\subsection td3 C. EvaltoDataMap
Next we need to link the data types to the evaluation type:

\code
     .
     .
     .
    // Residual (default scalar type is RealType)
    typedef Sacado::mpl::vector< RealType, 
				 MyVector<RealType>,
				 MyTensor<RealType> 
    > ResidualDataTypes;
  
    // Jacobian (default scalar type is Fad<double>)
    typedef Sacado::mpl::vector< FadType,
				 MyVector<FadType>,
				 MyTensor<FadType> 
    > JacobianDataTypes;

    // Maps the key EvalType a vector of DataTypes
    typedef boost::mpl::map<
      boost::mpl::pair<Residual, ResidualDataTypes>,
      boost::mpl::pair<Jacobian, JacobianDataTypes>
    >::type EvalToDataMap;
     .
     .
     .
\endcode

\subsection td4 D. Allocator
Define the Allocator type to use.  Phalanx comes with two allocators, but the user can write their own allocator class if these aren't sufficient.  The Phalanx allocator classes are:
- PHX::NewAllocator: uses the C++ "new" command to allocate each field on the heap separately.
- PHX::ContiguousAllocator: allocates a single contiguous block of memory on the heap for all fields regardless of the type.  This allows us to fit a subset of elements into cache to speed up the evaluation.

Code using the NewAllocator is:
\code
     .
     .
     .
    // ******************************************************************
    // *** Allocator Type
    // ******************************************************************
    typedef PHX::NewAllocator Allocator;
     .
     .
     .
\endcode

\subsection td5 E. EvalData,PreEvalData,PostEvalData
Users can pass their own data to the evaluate, preEvaluate and PostEvaluate methods of the PHX::FiledManager class.  In this example, the user passes in a vector of cell data objects that contain used defined auxiliary data fo each cell in the evaluation loop.  The pre and post evalution data is just to void - we are not using it in this example.  
\code
     .
     .
     .
    // ******************************************************************
    // *** User Defined Object Passed in for Evaluation Method
    // ******************************************************************
    typedef std::vector<CellData>& EvalData;
    typedef void* PreEvalData;
    typedef void* PostEvalData;

  };
\endcode

\section user_guide_step4 Step 4: Specialize the PHX::TypeString Object
For debugging information, Phalanx makes a forward declaration of the PHX::TypeString object.  This must be specialized for each evaluation type and each data type so that if there is a run-time error, phalanx can report detailed information on the problem.  We could have used the typeinfo from the stl, but the name() method is not demangled on every platform, so it can make debugging a challenge.  The specialized classes can go into their own file or can be added to the traits file above depending on how you use the Traits class.  During linking, if the compiler complains about multiple defninitions of your specialized traits classes, separate the traits implementation into their own file.

\code
namespace PHX {

  // ******************************************************************
  // ******************************************************************
  // Debug strings.  Specialize the Evaluation and Data types for the
  // TypeString object in the phalanx/src/Phalanx_TypeString.hpp file.
  // ******************************************************************
  // ******************************************************************

  // Evaluation Types
  template<> struct TypeString<MyTraits::Residual> 
  { static const std::string value; };

  template<> struct TypeString<MyTraits::Jacobian> 
  { static const std::string value; };

  const std::string TypeString<MyTraits::Residual>::value = 
    "Residual";

  const std::string TypeString<MyTraits::Jacobian>::value = 
    "Jacobian";

  // Data Types
  template<> struct TypeString<double> 
  { static const std::string value; };

  template<> struct TypeString< MyVector<double> > 
  { static const std::string value; };

  template<> struct TypeString< MyTensor<double> > 
  { static const std::string value; };

  template<> struct TypeString< Sacado::Fad::DFad<double> > 
  { static const std::string value; };

  template<> struct TypeString< MyVector<Sacado::Fad::DFad<double> > > 
  { static const std::string value; };

  template<> struct TypeString< MyTensor<Sacado::Fad::DFad<double> > > 
  { static const std::string value; };

  const std::string TypeString<double>::value = 
    "double";

  const std::string TypeString< MyVector<double> >::value = 
    "MyVector<double>";

  const std::string TypeString< MyTensor<double> >::value = 
    "MyTensor<double>";

  const std::string TypeString< Sacado::Fad::DFad<double> >::
  value = "Sacado::Fad::DFad<double>";
  
  const std::string TypeString< MyVector<Sacado::Fad::DFad<double> > >::
  value = "Sacado::Fad::DFad< MyVector<double> >";

  const std::string TypeString< MyTensor<Sacado::Fad::DFad<double> > >::
  value = "Sacado::Fad::DFad< MyTensor<double> >";
}
\endcode
  
\section user_guide_step5 Step 5: Write your Evaluators

For each field, you will need an evaluator.  Evaluators can evlauate multiple fields at the same time.  You can derive from the base class PHX::Evaluator, or you can derive from class PHX::EvaluatorWithBaseImpl that has most of the methods already implemented so that a no code is repeated in each evaluator.

An example for evaluating the density field is:

\code
#ifndef PHX_EXAMPLE_VP_DENSITY_HPP
#define PHX_EXAMPLE_VP_DENSITY_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_Field.hpp"

template<typename EvalT, typename Traits>
class Density : 
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits> {
  
public:
  
  Density(const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData ud);
  
private:
  
  typedef typename EvalT::ScalarT ScalarT;

  double constant;

  PHX::Field<ScalarT> density;
  PHX::Field<ScalarT> temp;

  std::size_t data_layout_size;

};

#include "Evaluator_Density_Def.hpp"

#endif
\endcode

Note that if you want to use the automated EvaluatorFactor builder to build an object of each evaluation type, you must derive from the PHX::EvaluatorDerived class.  This allows the variable manager to store a vector of base object pointers for each evaluation type in a single stl vector.

Also note that we pull the scalar type out of the evaluation type.

The implementation is just as simple:
\code
// **********************************************************************
template<typename EvalT, typename Traits> Density<EvalT, Traits>::
Density(const Teuchos::ParameterList& p) :
  density("Density", p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") ),
  temp("Temperature", p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") )
{ 
  this->addEvaluatedField(density);
  this->addDependentField(temp);
  this->setName("Density");
}

// **********************************************************************
template<typename EvalT, typename Traits>
void Density<EvalT, Traits>::
postRegistrationSetup(PHX::FieldManager<Traits>& vm)
{
  this->utils.setFieldData(density,vm);
  this->utils.setFieldData(temp,vm);

  data_layout_size = density.fieldTag().dataLayout().size();
}

// **********************************************************************
template<typename EvalT, typename Traits>
void Density<EvalT, Traits>::evaluateFields(typename Traits::EvalData d)
{ 
  std::size_t size = d.size() * data_layout_size;
  
  for (std::size_t i = 0; i < size; ++i)
    density[i] =  temp[i] * temp[i];
}
\endcode

The constructor pulls out data from the parameter list to set the correct data layout.  Additionally, it tells the FieldManager what fields it will evaluate and what fields it requires/depends on to perform the evaluation.

The postRegistrationSetup gets pointers from the FieldManager to the start of the array for each particular field.

\section user_guide_step6 Step 6: Implement the FieldManager in your code
  Adding the FieldManager to your code is broken into steps.  You must build each Evaluator for each field type, register the evaluators with the FieldManager, and then call the evaluate routines.  Continuing from our example:

\subsection fmd1 Building your Evaluators 
  Phalanx provides an automated builder PHX::EvaluatorFactory that will create an object of each evaluation type.  The only requirement is that the constructor of your object take one argument that is a Teuchos::ParamterList object.  This allows you to pass an arbitrary number of objects with arbitrary types.  If you use this factory, then you must provide a FactoryTraits class that gives the factory a list of it's evaluator types.  An example of the factory traits class is in the file FactoryTraits.hpp:

\code
#ifndef EXAMPLE_FACTORY_TRAITS_HPP
#define EXAMPLE_FACTORY_TRAITS_HPP

// mpl (Meta Programming Library) templates
#include "Sacado_mpl_vector.hpp"

// User Defined Evaluator Types
#include "Evaluator_Constant.hpp"
#include "Evaluator_Density.hpp"
#include "Evaluator_EnergyFlux_Fourier.hpp"
#include "Evaluator_FEInterpolation.hpp"
#include "Evaluator_NonlinearSource.hpp"


#include "boost/mpl/placeholders.hpp"
using namespace boost::mpl::placeholders;

/*! \brief Struct to define Evaluator objects for the EvaluatorFactory.
    
    Preconditions:
    - You must provide a Sacado::mpl::vector named EvaluatorTypes that contain all Evaluator objects that you wish the factory to build.  Do not confuse evaluator types (concrete instances of evaluator objects) with evaluation types (types of evaluations to perform, i.e., Residual, Jacobian). 

*/
template<typename Traits>
struct MyFactoryTraits {
  
  static const int id_constant = 0;
  static const int id_density = 1;
  static const int id_fourier = 2;
  static const int id_feinterpolation = 3;
  static const int id_nonlinearsource = 4;

  typedef Sacado::mpl::vector< Constant<_,Traits>,             // 0
 			       Density<_,Traits>,              // 1
 			       Fourier<_,Traits>,              // 2
 			       FEInterpolation<_,Traits>,      // 3
 			       NonlinearSource<_,Traits>       // 4
  > EvaluatorTypes;

};

#endif

\endcode 

Since the factory is built at compile time, we need to link a run-time choice to the compile-time list of object types.  Thus the user must provide the static const int identifiers that are unique for each type that can be constructed.  Users can ignore the "_" argument in the mpl vector.  This is a placeholder argument that allows us to iterate over each evaluation type.

Now let's build the evaluators.  We create a ParameterList for each evaluator and call the factory constructor.  This code can be found in the Example.cpp file in the directory "phalanx/example/EnergyFlux".

\code
      RCP<DataLayout> qp = rcp(new FlatLayout("QP", 4));
      RCP<DataLayout> node = rcp(new FlatLayout("NODE", 4));

      // Parser will build parameter list that determines the field
      // evaluators to build
      map<string, RCP<ParameterList> > evaluators_to_build;
      
      { // Temperature
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_constant;
	p->set<int>("Type", type);
	p->set<string>("Name", "Temperature");
	p->set<double>("Value", 2.0);
	p->set< RCP<DataLayout> >("Data Layout", node);
	evaluators_to_build["DOF_Temperature"] = p;
      }
      { // Density
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_density;
	p->set<int>("Type", type);
	p->set< RCP<DataLayout> >("Data Layout", qp);
	evaluators_to_build["Density"] = p;
      }

      { // Constant Heat Capacity
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_constant;
	p->set<int>("Type", type);
	p->set<string>("Name", "Heat Capacity");
	p->set<double>("Value", 2.0);
	p->set< RCP<DataLayout> >("Data Layout", qp);
	evaluators_to_build["Heat Capacity"] = p;
      }
      
      { // Nonlinear Source
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_nonlinearsource;
	p->set<int>("Type", type);
	p->set< RCP<DataLayout> >("Data Layout", qp);
	evaluators_to_build["Nonlinear Source"] = p;
      }

      { // Fourier Energy Flux
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_fourier;
	p->set<int>("Type", type);
	p->set< RCP<DataLayout> >("Data Layout", qp);
	evaluators_to_build["Energy Flux"] = p;
      }

      { // FE Interpolation
	RCP<ParameterList> p = rcp(new ParameterList);

	int type = MyFactoryTraits<MyTraits>::id_feinterpolation;
	p->set<int>("Type", type);

	p->set<string>("Node Variable Name", "Temperature");
	p->set<string>("QP Variable Name", "Temperature");
	p->set<string>("Gradient QP Variable Name", "Temperature Gradient");

	p->set< RCP<DataLayout> >("Node Data Layout", node);
	p->set< RCP<DataLayout> >("QP Data Layout", qp);

	evaluators_to_build["FE Interpolation"] = p;
      }

      // Build Field Evaluators for each evaluation type
      EvaluatorFactory<MyTraits,MyFactoryTraits<MyTraits> > factory;
      RCP< vector< RCP<Evaluator_TemplateManager<MyTraits> > > > 
	evaluators;
      evaluators = factory.buildEvaluators(evaluators_to_build);
 
      // Create a FieldManager
      FieldManager<MyTraits> fm;

      // Register all Evaluators 
      registerEvaluators(evaluators, fm);

\endcode

Next tell the Field Manager which quantities it need to provide to your code:

\code
      // Request quantities to assemble RESIDUAL PDE operators
      {
	typedef MyTraits::Residual::ScalarT ResScalarT;
	Tag< MyVector<ResScalarT> > energy_flux("Energy_Flux", qp);
	fm.requireField<MyTraits::Residual>(energy_flux);
	Tag<ResScalarT> source("Nonlinear Source", qp);
	fm.requireField<MyTraits::Residual>(source);
      }

      // Request quantities to assemble JACOBIAN PDE operators
      {
	typedef MyTraits::Jacobian::ScalarT JacScalarT;
	Tag< MyVector<JacScalarT> > energy_flux("Energy_Flux", qp);
	fm.requireField<MyTraits::Jacobian>(energy_flux);
	Tag<JacScalarT> source("Nonlinear Source", qp);
	fm.requireField<MyTraits::Jacobian>(source);
      }
\endcode

Once the evaluators are registered with the FieldManager and it knows which field it needs to provide, call the postRegistrationSetup method to figure out the total amount of memory to allocate for all fields.  It will also figure out which evaluators to call and the order to call them in.  This will require the user to specify the maximum number of cells for each evaluation.  This number needs to be selected so that all fields can fit in the cache of the processor if possible.  

\code
      const std::size_t max_num_cells = 100;
      fm.postRegistrationSetup(max_num_cells);
\endcode


Finally, users can call the evlauate routines and the pre/post evaluate routines if reuqired.

\code
      std::vector<CellData> cells(max_num_cells);

      fm.preEvaluate<MyTraits::Residual>(NULL);
      fm.evaluateFields<MyTraits::Residual>(cells);
      fm.postEvaluate<MyTraits::Residual>(NULL);
\endcode

Accessing field data is achieved as follows:

\code
      Field< MyVector<double> > ef("Energy_Flux", qp);
      fm.getFieldData<MyVector<double>,MyTraits::Residual>(ef);
\endcode


*/

/* ************************************************************************ */
/* ************************************************************************ */

/*! \page faq Frequently Asked Questions

\section Questions
\ref faq1

\ref faq2

\ref faq3

\ref faq4

\section faq1 1. Why name it Phalanx?  
The phalanx was one of the most dominant military formations of the Greek armies in the classical period.  It was a strictly ordered formation.  The Phalanx software library figures out ordered dependencies of field evaluators.   A second more obscure reference relates to the US Navy.  The Phalanx software package was designed to provide nonlinear functionality for the <a href="http://trilinos.sandia.gov/packages/intrepid">Intrepid</a> itegration library.  Intrepid was the name of an aircraft carrier during in WW II.  Modern US aircraft carriers are protected by a Close-In Weapons System (CIWS) named the Phalanx.  Finally, the PI of this project is an avid strategy warfare gamer and leans towards military references.

\section faq2 2. Why does the evaluation type exist?
This allows for multiple evaluations that used the same Scalar type, but different evaluators specific to the CalculationType.  For example, evaluating the Jacobian (derivative of the residual with respect to the degrees of freedom) and the parameter sensitivities (derivative of the residual with respect to the parameters) both use the same Sacado data type.  They both want derivatives but with respect to a different set of entities.

\section faq3 3. Why do evaluation types have arbitrary data types.  Shouldn't all fields of an evaluation type require the same scalar type?
In 99% of the use cases, the answer is yes!  By changing the scalar type for a specific type of computation, you can break the dependency chain and get incorrect answers.  For example, if computing a Jacobian, you would use a FAD scalar type instead of a double scalar type for the evaluation.  By mixing FAD and double fields in an evaluation type, you will lose derivative components and get incorrect sensitivities.  But there are some use cases where one might want to use different data types.  If the user knows that certain evaluation routines have no derivative components, they can bypass the AD chain by switching to doubles for that part of the computation.  The is an expert only user mode, but we allow for it.

\section faq4 4. Why are the evaluators templated on the evaluation type?
Each evaluator must be templated on the evaluation type to be able to (1) allow for the expert use mode described in FAQ \ref faq3; (2) use the automated factory class (PHX::EvaluatorFactory) to build the evaluator for each evalution type.  

*/

/* ************************************************************************ */
/* ************************************************************************ */

/*! \page developer_notes Developer Notes

\section dnotes_notz Description

This is a description Pat Notz called his elevator speech on Phalanx:

\verbatim
The {NAME} package is a library to help decompose a complex problem
into a set of simpler problems with managed dependencies. Two benefits
of decomposing a problem using {NAME} are (1) increased flexibility
because each simpler piece becomes an extension point that can be
swapped out with different implementations and (2) easier to craft
code because each piece is simpler, more focused and easier to test in
isolation.  The counter-benefits to these are (1) a potential
performance loss due to fragmentation of the over-all algorithm (e.g.,
several small loops instead of one large loop) and (2) a potential
loss of visibility of the original, composite problem (since the code
is scattered into multiple places).  Managing these trade-offs can
result in application code that both performs well and supports rapid
development and extensibility.  
\endverbatim

\section dnotes_contiguous_allocator Contiguous Allocator

Simple concept for chunk-layout of heterogeneous members with proper
data alignment.

A 'malloc'ed chunk of memory is always suitable to hold any data type,
or array of data types (do a 'man malloc').

\verbatim
        unsigned char * chunk = malloc( chunk_size );
\endverbatim

If a 'double' is to occupy a portion of the chunk then is must be
aligned such that:

\verbatim
        double * d = (double *)( chunk + offset_d );
\endverbatim

where:

\verbatim
        0 == offset_d % sizeof(double);
\endverbatim

The key point is for the offset of each member within the chunk to
have the correct alignment, and of course not overlap any other
member.

\section dnotes_old_description Old email with description of variable manager 

Email description of Phalanx (formerly known as the variable manager):

Phalanx is a tool for controlling the evaluation of
data used in the assembly of a residual/Jacobian fill for PDE solvers.
I implemented the one in Charon and one also exists in Aria (although it
is used for a different purpose).  The idea behind the VM is to allow
for fast changes to dependency trees.  The VM provides variable
quantities.  A variable consists of a name (string), a type (scalar,
vector, tensor - all templated for ad), and the place that it lives
(quadrature point, node, edge, face, element, etc...).  A variable is
evaluated using a provider.  The provider contains lists of variables it
provides (evlauates) and a list of variables it requires to perform the
evaluation (its dependencies).  The VM figures out the correct order to
call the providers to get the dependency chain correct.

I developed this for Charon to address the following issues:

1. Determine the order of variable evaluation so that the correct
evaluation is performed: For example, to compute the density, it could
be a constant, it could be a function of temperature or it could be a
function of temperature and pressure.  A factory creates the density
provider and from this, the manager automatically figures out the
correct evaluation order - whether T and/or P will have to be
evaluated before the density.  This is an easy example, but there are
much more complex dependency chains.

2. Provides a single place to store all evaluation variables so that the
procedure by which the variable is evaluated is not known by the objects
using the variable.  For example, in Charon a quantity of interest is
the velocity components.  The velocity could be a degree of freedom or
it could be a user specified constant, or it could be a user specified
field (read from an exodus file).  Before we had the variable manager,
everywhere in Charon that used velocity, there were repeated loops that
looked like:

\verbatim
if (velocity is DOF)
  // get values from dof structure
else if (velocity is auxiliary data)
  // go get from aux data structure
else if (...)
\endverbatim

There were hard coded loops in multiple places that had to be updated
every time we added a new way to specify velocity.   Now, when velocity
is required, you just go to the VM and get it.  No questions asked.

3. By having a single place to store the variables, we eliminated the
reevaluation of temporary quantities that are used in different PDE
operators and different equations.  Prior to using the variable manager,
the reuse of variables like density was ignored. It was recomputed for
each operator evaluation.  After implementing the variable manager, my
residual evaluation times were reduced by a factor of 7 for reacting
flow problems!

4. By having the evaluation order, one can perform automatic
differentiation without any library like saccado.  When the user writes
a provider for a function, they can also write a local Jacobian for the
individual function.  Then by using the chain rule, the local element
Jacobian can be assembled by traversing the dependency tree generated by
the variable manager.  Ariawrote their VM to do this, since Sacado was
not around when they started.  Aria does use Sacado for user defined
functions, but uses the VM for internal functions.

There are other benefits, but you get the idea.  A primary goal of my
current ASC funding is to generalize the variable manager in Charon
(there are many things I can do better now to make it more efficient and
easier to use) as a trilinos package to be used with Intrepid to aid in
element assembly for nonlinear equation sets.

  */

/* ************************************************************************ */
/* ************************************************************************ */

/*! \page junk Junk

  \todo Add a configure check for BOOST.

*/

/* ************************************************************************ */
/* ************************************************************************ */

#endif
