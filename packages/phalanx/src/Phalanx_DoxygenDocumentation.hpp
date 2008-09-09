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

Phalanx is a local assembly kernel specifically designed for general
partial differential equation (PDE) solvers.  It can be used with any
cell-based discretization techniques including finite element and
finite volume.  The main goal of Phalanx is to decompose a complex
problem into a number of simpler problems with managed dependencies to
support rapid development and extensibility of the PDE code.  This
approach, coupled with the template capabilities of C++ offers a
number of unique and powerful capabilities:

<ul>
<li> Fast Integration with Flexible and Extensible Models:

     <ul> 
     
     <li> Increased flexibility because each simpler piece of the
     decomposed problem becomes an extension point that can be swapped
     out with different implementations.

     <li> Easier to implement new code because each piece is simpler,
     more focused and easier to test in isolation.

     <li> Easier for users to add their own models.  While Phalanx is
     designed for maximum flexibility and efficiency, the interfaces
     exposed to users are very simple, requiring minimal training
     and/or knowledge of the C++ language to use.

     <li> Through the use of template metaprogramming concepts,
     Phalanx supports arbitrary user defined data types.  This allows
     for unprecedented flexibility for direct integration with user
     applications and opens the door to embedded technology.
   
     </ul>

<li> Support for Advanced Embedded Technology: Phalanx is fully
compatible with advanced embedded technologies. Embedded technology
consists of replacing the default scalar type (i.e., "double" values)
in a code with an object that overloads the typical mathematical
operators.  By replacing the scalar type, we can reuse the same code
base to produce different information.  For example, if we were to
compute a function \f$ F(x) \in R^n\f$, where \f$ F: R^n \rightarrow
R^n \f$ and \f$ x \in R^n \f$.  We could automatically compute the
sensitivities of the function, such as the Jacobian, \f$ J_{ij} =
\frac{\partial F_i}{\partial x_j} \f$ merely by replacing the scalar
type of double with an automatic differentiation (AD) object.  With
the C++ template mechanism, this is a trivial task.  Embedded
technology provides:

     <ul> 

     <li> Applications can reuse the exact same code base written to
     evaluate a function yet produce new information including
     arbitrary precision, sensitivities (via embedded automatic
     differentiation), and uncertainty quantification.  

     <li> By reusing the same code for both function and sensitivity
     evaluation, developers avoid having to hand code time consuming
     and error prone analytic derivatives and ensure that the equation
     sets are consistent.

     </ul>

<li> Consistent field evaluations via dependency chain management:
When users switch models, the dependencies for fields evaluated by the
model can change.  This can force field to be evaluated in a different
order based on the new dependencies.  Phalanx will automatically
perform the sorting and ordering of the evaluator routines.  For
example, if density were a field to be evaluated, many models that
have different dependencies.  Model A could depend on temperature and
pressure, whille model B could depend on temperature, pressure, and
species mass fractions of a multicomponent mixture.  The order of
evaluation of fields could change based on the change in dependencies.
Using Model B will require that the species mass fractions be
evaluated and available for the density evaluation while Model A does
not even require them.  This is a simple example, but the dependency
chain can become quite complex when solving thousands of coupled PDE
equations.

<li> Efficient evaluation of field data: Phalanx provides for an
arbitrary size block of cells to be evaluated on the local processor.
It evaluates all fields of interest at once for each block of cells.
By using the contiguous allocator and correctly sizing the number of
cells to fit in processor cache, one can arrange all fields to exist
in a contiguous block of memory, independent of the data type objects
used to store the data (one must still be careful of data alignement
issues).  By keeping all fields in cache, the code will run much
faster.  This blocking will also, in the future, allow for multi-core
distrubution of the cell evaluations.

</ul>

Phalanx is a hammer.  It's use should be carefully considered.  We
recommend its use when writing a general PDE framework where one needs
support for flexibility in equation sets and material properties in
those sets.  It should not be used for simple sets of PDEs where the
equations rarely change.  There are some drawbacks to using Phalanx
that should be considered:

- A potential performance loss due to fragmentation of the over-all
algorithm (e.g., several small loops instead of one large loop).  A
judicous choice of field variables can alleviate this problem.

- A potential loss of visibility of the original, composite problem
 (since the code is scattered into multiple places).  

Managing these trade-offs can result in application code that both
performs well and supports rapid development and extensibility.

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

\section authors Authors and Contributors

The following have contributed to the design through ideas, discussions, and/or code development of Phalanx:

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

\subsection ug_dummy_1 A. Understand Templates
Phalanx is a complex package that make heavy use of the C++ templating
mechanism.  We recommend that users of the Phalanx package first
familiarize themselves with C++ templates.  An excellent reference is
"C++ Templates: The Complete Guide" by Vandevoorde and Josuttis, 2003.
While users do not need to understand template metaprogramming to use
Phalanx, the concepts underlying many operations in Phalanx can be
found in "C++ Template Metaprogramming" by Abrahams and Gurtovoy,
2004 and the Boost template metaprogramming library (MPL).

Once Phalanx is integrated into a code, the actual addition of
Evalators requires very little knowledge of templates.  In fact if the
user does a cut-and-paste of a current evaluator and makes the simple
modifications, they won't have to know anything about templates.  So
once Phalanx is integrated into a code, the users writing evaluators
(i.e., new material models) need not know anything about templates.

\subsection ug_dummy_2 B. Learn the Phalanx Nomenclature
Users should then learn the nomeclature used in the package defined in
the \ref user_guide_domain_model.

\subsection ug_dummy_3 C. Tutorial

The main concept of Phalanx is to evaluate fields for solving PDEs.  We demonstrate the integration process using the simple example found in phalanx/example/EnergyFlux.  Suppose that we want to solve the heat equation over the physical space \f$ \Omega \f$:

\f[
  \nabla \cdot (-\rho C_p \nabla T) + s = 0
\f] 

where \f$s\f$ is a nonlinear source term.  The specific discretization technique whether finite element (FE) of finite volume (FV) will ask for \f$\mathbf{q}\f$ and \f$s\f$ at points on the cell.  Phalanx will evaluate \f$\mathbf{q}\f$ and \f$s\f$ at those points and return them to the discretization driver.

Using finite elements, we pose the problem in variational form:

Find \f$ u \in {\mathit{V^h}} \f$ and \f$ \phi \in {\mathit{S^h}} \f$ such that:

\f[
  - \int_{\Omega} \nabla \phi \cdot \mathbf{q} d\Omega 
  + \int_{\Gamma} \phi \mathbf{n} \cdot \mathbf{q} 
  + \int_{\Omega} \phi s d\Omega = 0 
\f]

where \f$ \mathbf{q} = -\rho C_p \nabla T \f$.






Follow the steps below to integrate Phalanx into your application.  The example code shown in the steps comes from the energy flux example in the directory "phalanx/example/EnergyFlux".  
- \ref user_guide_step1
- \ref user_guide_step2
- \ref user_guide_step3
- \ref user_guide_step4
- \ref user_guide_step5
- \ref user_guide_step6

\section user_guide_domain_model Phalanx Domain Model

<ul>
<li><b>%Cell</b>

Partial differential equations are solved in a domain.  This domain is
discretized into cells (also called elements for the finite element
method).  This library assumes that the block of cells being iterated
over is of the same type!  If different evaluators (i.e. different
material properties) are required in different blocks of cells, a new
FieldMangager must be used for each unique block of elements.  This is
required for efficiency.

<li><b>Scalar Type</b>

A scalar type, typically the template argument ScalarT in Phalanx
code, is the type of scalar used in an evaluation.  It is typically a
double, but can be special object types for embedded methods such as
sensitivity analysis.  For example, for sensitivity analysis, a double
scalar type is replaced with a foward automatic differentiation object
(FAD) or a reverse automatic differentaion object (RAD) to produce
sensitivity information.  Whatever type is used, the standard
mathematical operators are overloaded for the particular embedded
technology.  For an example of this, see the
<a href="http://trilinos.sandia.gov/packages/sacado">Sacado Automatic
Differentiation Library</a>. Some sample scalar types include:
<ul>
<li> float
<li> double
<li> Sacado::Fad::DFad<double> (for sensitivity analysis)
</ul>

<li><b>Algebraic Type</b>

An algebraic type is the type of objects that hold data.  It is usually a rank n tensor.  Simple examples include a scalar (rank-0 tensor), a vector (rank-1 tensor) or a matrix (rank-2 tensor).  It is not actually restircted to tensors, but can be any struct/class that a user implements.  The algebraic type is a description of how data is stored but does NOT have a corresponding type in the Phalanx code.  It is a notion or idea we use to describe a data type without specifying the actual scalar type (See "Data Type" for more information).  These types are defined by the user.  In the example in the directory "phalanx/example/EnergyFlux", the user selects three algebraic types to represent scalars, vectors and tensors.  The scalar algebraic type is equivalent to the scalar type used in the evaluation.  The vector and tensor objects are objects templated on the scalar type:

<ul>
<li> template<typename ScalarT> class MyVector { ... };
<li> template<typename ScalarT> class MyTensor { ... };
</ul>

In in a function evaluation routine templated on the scalar type, the code would look something like:
\code
template<typename ScalarT>
void myFunction() {

  ScalarT scalar_value;
  MyVector<ScalarT> vector_value;
  MyTensor<ScalarT> matrix_value;
    .
    .
    .
}
\endcode

<li><b>Data Type</b>

A data type, typically the template argument DataT in Phalanx code, is an actual type used for storing fields.  It is the combination of a scalar type and an algebraic type.  Some examples include:
<ul>
<li> double
<li> MyVector<double>
<li> MyTensor<double>
<li> Sacado::Fad::DFad<double>
<li> MyVector< Sacado::Fad::DFad<double> >
<li> MyTensor< Sacado::Fad::DFad<double> >
</ul>


<li><b>Evaluation Type</b>

The evaluation type, typically the template argument EvalT in Phalanx code, defines a unique type of evaluation to perform.  The user is free to choose the evaluation types and actually creates their own evaluation types.  An EvaluationContainer is allocated for each evaluation type specified in the users traits class.  Examples include:
<ul>
<li> Residual
<li> Jacobian
<li> ParameterSensitivity
</ul>

The evaluation type must be associated with one default scalar type and can optionally additional scalar types.  The scalar type usually determines what is being evaluated. For example, to evaluate the equation residuals, the scalar type is usually a double or float.  To evaluate a Jacobian, the scalar type could be a forward automatic differentiation object, Sacado::Fad::DFAD<double>.  By introducing the evaluation type in Phalanx, the same scalar type can be used for different evaluation types and can be specialized accordingly.  For example computing the Jacobian and computing parameter sensitivities both could use the Sacado::Fad::DFAD<double> scalar type.

<li><b>Storage</b>

A DataContainer object stores all fields of a particular data type.  Each EvaluationContainer holds a vector of DataContainers, one DataContainer for each vaid data type that is associated with that particular evaluation type.  One EvaluationContainer is constructed for each evaluation type.

<li><b>Data Layout</b>

The DataLayout object is used to distinguish fields with the same name, but exist at different discretization points in the cell.  For example, supposed we have written an evaluator the computes the "Density" field for a set of points in the cell.  Now we want to evaluate the density at a different set of points in the cell.  We might have a "Density" field in a cell associated with a set of integration points (quadrature points in finite elements) and another field associated with the nodes (nodal basis points in finite elements).  We use the same field name (so we can reuse the same Evaluator), "Density", but use two different DataLayouts, one for integration points and one for nodal point.  Now a FieldTag comparison will differentiate the fields due to the different DataLayout.  

Additionally, the DataLayout contains the number of DataT objects associated with the field in a single cell.  This size() parameter is not needed to distinguish uniqueness, since the number of objects can be the same for different fields.  It is stored here for convenience when figuring out the size of field arrays.

<li><b>Field Tag</b>

The FieldTag is a description of a field.  It is templated on the data type, DataT.  It is used to identify a field stored in the field manager.  It contains a unique identifier (an stl std::string) and a pointer to a data layout object.  If two FieldTags are equal, the DataLayout, the data type, and the string name are exactly the same.

<li><b>Field</b>

A Field is a set of values for a particular quantity of interest that must be evaluated.  It is templated on the data type, DataT.  It consists of a FieldTag and a reference counted smart pointer (Teuchos::RCP<DataT>) to the data array where the values are stored.  

<li><b>Evaluator</b>

An Evaluator is an object that evaluates a set of Fields.  It contains two vectors of Fields, one set is the fields it will evaluate and one set is the fields it depends on for the evaluation.  For example to evaluate a density field that is a function of temperature and pressure, the field the evaluator will evaluate (the first field set) is a density field, and the set of fields it requires to meet its dependencies are the temperature and pressure fields (the second field set).  The evaluator is templated on the evaluation type and the traits class.  

<li><b>Evaluator Manager</b>

The main object that stores all Fields and Evaluators.  The evaluator manager (EM) sorts the evaluators and determines which evaluators to call and the order necessary to ensure consistency in the fields.  The EM also allocates the memory for storage of all fields so that if we want to force all fields to be in a contiguous block, we have that option available.  Users can write their own allocator for memory management.

<li><b>Local</b> 

The term local refers to having all information on a single processor in a multiprocessor run.
Phalanx does local evaluation of fields, meaning that all information
  required for the evaluation is stored locally on that processor (The user must handle pulling the info into the processor if distributed).  This
does NOT in any way limit a user to serial runs.  Phalanx is used
routinely for large-scale parallel distributed architecure codes.  The
design of Phalanx will also allow for users to take advantage of multi-core
architectures as well through a variety of implementations.

</ul>

\section user_guide_mdarray_domain_model Multi-Dimensional Array Domain Model

Document has not been publicly released.  Will add a reference when available.

\section user_guide_step1 Step 1: Configuring, Building, and installing Phalanx

Phalanx is distributed as a package in the <a href="http://trilinos.sandia.gov">Trilinos Framework</a>.  It can be enabled as part of a trilinos build with the configure option "--enable-phalanx".  Phalanx currently has direct dependencies on the following third party libraries:

- Requires the <a href="http://trilinos.sandia.gov/packages/teuchos">Teuchos</a> utilities library, part of the <a href="http://trilinos.sandia.gov/">Trilinos Framework</a>.  This will automatically be enabled when you enable the phalanx library.
 
 - Requires the <a href="http://trilinos.sandia.gov/packages/sacado">Sacado Automatic Differentiation Library</a>, part of the <a href="http://trilinos.sandia.gov/">Trilinos Framework</a>.  This will automatically be enabled when you enable the phalanx library.

 - Requires the <a href="http://www.boost.org">Boost Template Metaprogramming (MPL) Library</a>.  You must add the path to the Boost library during Trilinos configuration using the flag "--withincdirs=<path>".

 - Optional: Some performance tests run comparisons against <a href="http://tvmet.sourceforge.net/">TVMET: Tiny Vector Matrix library using Expression Templates</a>.  This is to get a feel for how our "dumb" vector matrix objects perform compared to expression templates.  Use the configure option "--with-tvmet" to enable the tvmet functionality in the performance tests.  You must add the path to the TVMET library during Trilinos configuration using the flag "--withincdirs=<path>".

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

<b> The Traits struct must define each of the following typedef members of the struct:</b>

- \b EvalTypes - an mpl::vector of user defined evaluation types.  Each evaluation type must have a typedef member called ScalarT that provides the default scalar type.  This is used to automate the building of evaluators for each evaluation type using the EvaluatorFactory.
      
- \b EvalToDataMap - an mpl::map.  The key is an evaluation type and the value is an mpl::vector of valid data types for that particular evaluation type.
      
- \b Allocator type - type that defines the allocator class to use to allocate the memory for data storage.
      
- \b EvalData - A user defined type to be passed in to the evaluateFields() call.  Allows users to pass in arbitrary data.

- \b PreEvalData - A user defined type to be passed in to the preEvaluate() call.  Allows users to pass in arbitrary data.

- \b PostEvalData - A user defined type to be passed in to the postEvaluate() call.  Allows users to pass in arbitrary data.

\subsection td1 A. Basic Class
The basic outline of the traits struct is:
\code
struct MyTraits : public PHX::TraitsBase {
    .
    .
    .
};
\endcode

Inside this struct we need to implement all the typedefs listed above.  The example we will follow is in the file Traits.hpp in "phalanx/example/EnergyFlux" directory.

\subsection td2 B. EvalTypes
First we need to know the evaluation types.  Each evaluation type must include at least a typedef'd default scalar type argument as a public member called ScalarT.  Here is the example code:

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

The typedefs ReatType and FadType are done only for convenience.  They are not actually required.  Only the EvalTypes typedef is required in the code above.

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

Before Writing your evaluators, you must decide on how you plan to build the evaluators.  In most cases you will want to build one Evaluator for each evaluation type.  Phalanx provides an automated factory called the PHX::EvaluatorFactory that will build an evaluator for each evaluation type automatically, but this places a restriction on the constructor of all evaluators built this way.  If you plan to use the automated builder, the constructor for the Evaluator must contain only one argument - a Teuchos::ParameterList.  This paramter list must contain a key called "Type" with an integer value corresponding to the type of Evaluator object to build.  The parameterlist can contain any other information the user requires for proper construction of the Evaluator.  You are not restricted to using the automated factory for every Evalautor.  You can selectively use the automated factory where convenient.

For each field, you will need an evaluator.  Evaluators can evlauate multiple fields at the same time.  You can derive from the base class PHX::Evaluator, or you can derive from class PHX::EvaluatorWithBaseImpl that has most of the methods already implemented so that the same support code is not replicted in each evaluator.  We STRONGLY recommend deriving from the class PHX::EvaluatorWithBaseImpl.

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

Note that if you want to use the automated factory PHX::EvaluatorFactory to build an object of each evaluation type, you must derive from the PHX::EvaluatorDerived class as shown in the example above.  This allows the variable manager to store a vector of base object pointers for each evaluation type in a single stl vector.

  Also note that we pull the scalar type, ScalarT, out of the evaluation type.

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

The postRegistrationSetup gets pointers from the FieldManager to the array for storing data for each particular field.

  Writing evaluators can be tedious.  We have invested much time in minimizing the amount of code a user writes for a new evaluator.  Our experience is that you can literally have hundreds of evaluators.  So we have added macros to hide the boilerplate code in each evaluator.  Not only does this streamline/condense the code, but it also hides much of the templating.  So if your userbase is uncomfortable with C++ templates, the macro definitions could be very helpful.  The definitions are found in the file Phalanx_Evaluator_Macros.hpp.  The same evaluator shown above is now implemented using the macro definitions:

Class declaration:
\code
#ifndef PHX_EXAMPLE_VP_DENSITY_HPP
#define PHX_EXAMPLE_VP_DENSITY_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"

PHX_EVALUATOR_CLASS(Density)

  double constant;

  PHX::Field<ScalarT> density;
  PHX::Field<ScalarT> temp;

  std::size_t data_layout_size;

PHX_EVALUATOR_CLASS_END

#include "Evaluator_Density_Def.hpp"

#endif
\endcode

Class definition:
\code
//**********************************************************************
PHX_EVALUATOR_CTOR(Density,p) :
  density("Density", p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") ),
  temp("Temperature", p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") )
{ 
  this->addEvaluatedField(density);
  this->addDependentField(temp);
  this->setName("Density");
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(Density,fm)
{
  this->utils.setFieldData(density,fm);
  this->utils.setFieldData(temp,fm);

  data_layout_size = density.fieldTag().dataLayout().size();
}

//**********************************************************************
PHX_EVALUATE_FIELDS(Density,d)
{ 
  std::size_t size = d.size() * data_layout_size;
  
  for (std::size_t i = 0; i < size; ++i)
    density[i] =  temp[i] * temp[i];
}
\endcode

The evaluators for the example problem in "phalanx/example/EnergyFlux" have been rewritten using the macro definitions in the directory "phalanx/test/Utilities/evaluators".  

\section user_guide_step6 Step 6: Implement the FieldManager in your code
  Adding the FieldManager to your code is broken into steps.  You must build each Evaluator for each field type, register the evaluators with the FieldManager, and then call the evaluate routines.  Continuing from our example:

\subsection fmd1 A. Building your Evaluators 
Users can build their own Evaluators and register them with the FieldManager or they can use our automated factory, PHX::EvaluatorFactory to handle this for them.  Normally, users will want to build an evaluator for each evaluation type.  The factory makes this very easy.  Additionally, you are not restricted to using the automated factory for every Evalautor.  You can selectively use the automated factory where convenient.  The next two sections describe how to build the evaluators.

\subsubsection fmd1s1 A.1 Building and Registering Evaluators Manually.
  To build an Evaluator manually, all you need to do is allocate the Evaluator on the heap using a reference counted smart pointer (Teuchos::RCP) to point ot the object.  This will ensure proper memory management of the object.  Here is an example of code to build a Density evaluator for each evaluation type and register it with the corresponding manager.

\code
// Create a FieldManager
FieldManager<MyTraits> fm;

// Constructor requirements
RCP<ParameterList> p = rcp(new ParameterList);
p->set< RCP<DataLayout> >("Data Layout", qp);
     
// Residual
Teuchos::RCP< Density<MyTraits::Residual,MyTraits> > residual_density = 
  Teuchos::rcp(new Density<MyTraits::Residual,MyTriats>(p));

fm.registerEvaluator<MyTraits::Residual>(residual_density);

// Jacobian
Teuchos::RCP< Density<MyTraits::Jacobian,MyTraits> > jacobian_density = 
  Teuchos::rcp(new Density<MyTraits::Jacobian,MyTriats>(p));

fm.registerEvaluator<MyTraits::Residual>(jacobian_density);
\endcode

As one can see, this becomes very tedious if there are many evaluation types.  It is much better to use the automated factory to build one for each evalaution type.  Where this method is useful is if you are in a class already templated on the evaluation type, and would like to build and register an evaluator in that peice of code.

\subsubsection fmd1s2 A.2 Using the Automated Factory

  Phalanx provides an automated builder PHX::EvaluatorFactory that will create an object of each evaluation type.  The following requirements are placed on each and every Evaluator that a user writes if they plan to use the automated factory:
<ol>

<li>  The constructor of your Evaluator must take exactly one argument that is a Teuchos::ParamterList object.  A Teuchos::ParameterList allows you to pass an arbitrary number of objects with arbitrary types.

<li>  In the Teuchos::ParamterList, you must have one key called "Type" that is associated with an integer value that uniquely corresponds to an Evaluator object written by the user.  This integer is defined in the users factory traits object described below.

<li>  The Evaluator must derive from the PHX::EvaluatorDerived class as shown in the example above.  This allows the variable manager to store a vector of base object pointers for each evaluation type in a single stl vector yet be able to return the derived class.

</ol>


To build an PHX::EvaluatorFactory, you must provide a factory traits class that gives the factory a list of its evaluator types.  An example of the factory traits class is the FactoryTraits object in the file FactoryTraits.hpp:

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

Since the factory is built at compile time, we need to link a run-time choice to the compile-time list of object types.  Thus the user must provide the static const int identifiers that are unique for each type that can be constructed.  Users can ignore the "_" argument in the mpl vector.  This is a placeholder argument that allows us to iterate over and instert each evaluation type into the factory traits.

Now let's build the evaluators.  The following code can be found in the Example.cpp file in the directory "phalanx/example/EnergyFlux".  We create a ParameterList for each evaluator and call the factory constructor.  Since we are using the factory, we need to specify the "Type" argument set to the integer of the corresponding Evaluator in the factory traits that we wish to build:

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

The map "evaluators_to_build" has a key of type "std::string".  This key is irrelevant to the execution of the code.  It is an identifer that can help users debug code or search for a specific provider in the list.  What you put in the key is for your own use.

You are free to register evaluators until the time you call the method postRegistrationSetup() on the field manager.  Once this is called, you can not register any more evaluators.  You can make multiple registration calls to add more providers - there is no limit on the number of calls.

\subsection fmd2 B. Request which Fields to Evaluate

You must tell the FieldManager which quantities it should evaluate.  This step can occur before, after, or in-between the registration of Evaluators.  You must request variables separately for each evaluation type.  This is required since since data types can exist in multiple evaluation types.

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

\subsection fmd3 C. Call FieldManager::postRegistrationSetup()

Once the evaluators are registered with the FieldManager and it knows which field it needs to provide, call the postRegistrationSetup() method on the FieldManager.  This method requires that the user specify the maximum number of cells for each evaluation.  This number should be selected so that all fields can fit in the cache of the processor if possible.  This method causes the following actions to take place in the FieldManager:

<ol>
<li> Based on the requested fields in \ref fmd2, the FieldManager will trace through the evaluators to determine which evaluators to call and the order in which they need to be called to achieve a consistent evaluation. Not all evaluators that are registered will be used.  They will only be called to satisfy dependencies of the required fields.
<li> Once the dependency chain is known, we can pull together a flat list of all fields that will be used.  Now the FieldManager will allocate memory to store all fields.  It will use the Allocator object from the users traits class to accomplish this.  By unifying the allocation into a single object, we can force all fields to be allocated in a single contiguous block of memory if desired.
<li> Once the memory for field data is allocated, we must set the pointer to that memory block inside each Field or MDField object in each evaluator.  The FieldManager does this by calling the postRegistrationSetup() method on each evaluator that is required for the computation. 
</ol>

\code
      const std::size_t max_num_cells = 100;
      fm.postRegistrationSetup(max_num_cells);
\endcode


\subsection fmd4 D. Call evaluate()

Finally, users can call the evlauate routines and the pre/post evaluate routines if required.

\code
      std::vector<CellData> cells(max_num_cells);

      fm.preEvaluate<MyTraits::Residual>(NULL);
      fm.evaluateFields<MyTraits::Residual>(cells);
      fm.postEvaluate<MyTraits::Residual>(NULL);
\endcode

If you have more cells on a processor than you can fit into cache (the number of cells you set in the postRegistrationSetup call in \ref fmd4), you can perform multiple evaluate calls for each set:

\code
      int num_cells = 150;
      int max_num_cells_in_cache = 30;

      std::vector< std::vector<CellData> > cells(5);
      std::vector<CellData> cache_cells(max_num_cells_in_cache);

      fm.preEvaluate<MyTraits::Residual>(NULL);

      for (
      fm.evaluateFields<MyTraits::Residual>(cells);

      fm.postEvaluate<MyTraits::Residual>(NULL);
\endcode



\subsection fmd5 E. Accessing Data

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
