// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
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

- \ref elevator_speech
- \ref overview
- \ref news
- \ref user_guide
- \ref faq
- \ref bugs
- \ref history
- \ref authors
- \ref copyright
- \ref questions

\section elevator_speech The Elevator Speech

Phalanx is a local field evaluation kernel specifically designed for
general partial differential equation solvers. The main goal of
Phalanx is to decompose a complex problem into a number of simpler
problems with managed dependencies to support rapid development and
extensibility of the PDE code. Through the use of template
metaprogramming concepts, Phalanx supports arbitrary user defined data
types and evaluation types. This allows for unprecedented flexibility
for direct integration with user applications and provides extensive
support for embedded technology such as automatic differentiation for
sensitivity analysis, optimization, and uncertainty quantification.

\section overview Overview

Phalanx is a local field evaluation kernel specifically designed for
general partial differential equation (PDE) solvers.  It can be used
with any cell-based discretization techniques including finite element
and finite volume.  The main goal of Phalanx is to decompose a complex
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
example, if density were a field to be evaluated, different models
could have different dependencies.  Model A could depend on
temperature and pressure, whille model B could depend on temperature,
pressure, and species mass fractions of a multicomponent mixture.  The
order of evaluation of fields could change based on the change in
dependencies.  Using Model B will require that the species mass
fractions be evaluated and available for the density evaluation while
Model A does not even require them.  This is a simple example, but the
dependency chain can become quite complex when solving thousands of
coupled PDE equations.

<li> Efficient evaluation of field data: Phalanx was designed on the
idea of worksets.  A workset is an arbitrarily sized block of cells to
be evaluated on the local processor.  Phalanx evaluates all fields of
interest at once for each block of cells.  By using the contiguous
allocator and correctly sizing the workset to fit in processor cache
(if possible), one can arrange all fields to exist in a contiguous
block of memory, independent of the data type objects used to store
the data (one must still be careful of data alignement issues).  By
keeping all fields in cache, the code should run much faster.  This
workset idea will also, in the future, allow for multi-core
distrubution of the cell evaluations.

</ul>

Phalanx is a hammer.  It's use should be carefully considered.  We
recommend its use when writing a general PDE framework where one needs
support for flexibility in equation sets.  It should not be used for
simple sets of PDEs where the equations rarely change.  There are some
drawbacks to using Phalanx that should be considered:

- A potential performance loss due to fragmentation of the over-all
algorithm (e.g., several small loops instead of one large loop).  A
judicous choice of field variables can alleviate this problem.

- A potential loss of visibility of the original, composite problem
 (since the code is scattered into multiple places).  

Managing these trade-offs can result in application code that both
performs well and supports rapid development and extensibility.

\section news News

Two new capabilities are in development, but not yet released:
<ul>
<li> Expression templates:  We are exploring expression templates for the field objects.
<li> We are exploring hybrid CPU/GPU support.  
</ul>

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

  - Roger Pawlowski (Lead Developer), SNL 01414
  - Eric Phipps, SNL 01411
  - Pat Notz, SNL 01541 

\section copyright Copyright

\verbinclude copyright.txt

\section questions For All Questions and Comments...
  
   Please contact Roger Pawlowski (rppawlo@sandia.gov).

*/

/* ************************************************************************ */
/* ************************************************************************ */

/*! \page user_guide Users Guide

\section user_guide_index Index

- \ref user_guide_getting_started
- \ref user_guide_domain_model
- \ref user_guide_mdarray_domain_model
- \ref performance
- \ref user_guide_step1
- \ref user_guide_step2
- \ref user_guide_step3
- \ref user_guide_step4
- \ref user_guide_step5
- \ref user_guide_step6
- \ref questions

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

Once Phalanx is integrated into a code by a template savvy developer,
the only operation application users should perform is to extend the
evaluation routines to new equations and new models by adding new
Evaluators.  We realize that application developers and application
users have two distictly different skill sets. Therefore much design
work has been invested in making the Evaluator implementation very
clean and as template free as possible.  Therefore, users of the code
who only write new Evaluators DO NOT need to know templates.

\subsection ug_dummy_2 B. Learn the Phalanx Nomenclature
Users should then learn the nomeclature used in the package defined in
the \ref user_guide_domain_model.

\subsection ug_dummy_3 C. Tutorial

The main concept of Phalanx is to evaluate fields for solving PDEs
(although it is not limited to this).  We demonstrate the integration
process using the simple example found in phalanx/example/EnergyFlux.
Suppose that we want to solve the heat equation over the physical
space \f$ \Omega \f$:

\f[
  \nabla \cdot (-\rho k \nabla T) + s = 0
\f] 

where \f$ T \f$ is the temparature (and our degree of freedom), \f$ \rho \f$ is the density, \f$ k \f$ is the thermal conductivity, and \f$s\f$ is a nonlinear source term.  We pose this in terms of a conservation law system:

\f[
  \nabla \cdot (\mathbf{q}) + s = 0
\f]

where \f$ \mathbf{q} = -\rho k \nabla T \f$ is the heat flux.  The specific discretization technique whether finite element (FE) or finite volume (FV) will ask for \f$\mathbf{q}\f$ and \f$s\f$ at points on the cell.  Phalanx will evaluate \f$\mathbf{q}\f$ and \f$s\f$ at those points and return them to the discretization driver.

Using finite elements, we pose the problem in variational form:  Find \f$ u \in {\mathit{V^h}} \f$ and \f$ \phi \in {\mathit{S^h}} \f$ such that:

\f[
  - \int_{\Omega} \nabla \phi \cdot \mathbf{q} d\Omega 
  + \int_{\Gamma} \phi \mathbf{n} \cdot \mathbf{q} 
  + \int_{\Omega} \phi s d\Omega = 0 
\f]

Phalanx will evaluate \f$\mathbf{q}\f$ and \f$s\f$ at the quadrature
points of the cells and pass them off to the integrator such as <a
href="http://trilinos.sandia.gov/packages/intrepid">Intrepid</a>.

This is a trivial example, but the dependency chains can grow quite
complex if performing something such as a reacting flow calculation
coupled to Navier-Stokes.

Follow the steps below to integrate Phalanx into your application. The example code shown in the steps comes from the energy flux example in the directory "phalanx/example/EnergyFlux".  Note that many classes are named with the word "My" such as MyWorkset, MyTraits, and MyFactory traits.  Any object that starts with the word "My" denotes that this is a user defined class.  The user must implement this class specific to their application.  All Evaluator derived objects are additionally implemented by the user even though they do not follow the convention of starting with the word "My".

- \ref user_guide_step1
- \ref user_guide_step2
- \ref user_guide_step3
- \ref user_guide_step4
- \ref user_guide_step5
- \ref user_guide_step6

\section user_guide_domain_model Phalanx Domain Model

<ul>
<li><b>Cell</b>

Partial differential equations are solved in a domain.  This domain is
discretized into cells (also called elements for the finite element
method).  This library assumes that the block of cells being iterated
over is of the same type!  If different evaluators (i.e. different
material properties) are required in different blocks of cells, a new
FieldMangager must be used for each unique block of elements.  This is
required for efficiency.

<li><b>Parallel and Serial Architectures</b> 

Phalanx can be used on both serial and multi-processor architectures.
The library is designed to perform "local" evalautions on a "local"
set of cells.  The term local means that all cell and field data
required for an evaluation is on the processor that the evaluation is
executed.  So for parallel runs, the cells are distributed over the
processors and a FieldManager is built on each processor to evaluate
only the cells that are assigned to that processor.  If there is any
data distributed to another processor that is required for the
evaluation, the user must handle pulling that information on to the
evaluation processor.  The design of Phalanx will also allow for users
to take advantage of multi-core architectures through a variety of
implementations.

<li><b>%Workset</b>

For efficiency, the evaluation of fields on a set of cells can be
divided into worksets.  The goal of using worksets is to fit all the
required fields into the processor cache so that an evaluation is not
slowed down by paging memory.  Suppose we have 2020 cells to evaluate
on a 4 processor machine.  We might distribute the load so that 505
cells are on each processor.  Now the user must figure out the workset
size.  This is the number of cells to per evaluation call so that the
field memory will fit into cache.  If we have 505 cells on a
processor, suppose we find that only 50 cells at a time will fit into
cache.  Then we will create a FieldManager with a size of 50 cells.
This number is specified in the construction of data layouts.  During the call to postRegistrationSetup(), the FieldManager will allocate workspace storage for all fields relevant to the evaluation.

For our example, there will be 11 worksets.  The first 10 worksets
will have the 50 cell maximum and the final workset will have the 5
remaining cells.  The evaluation routine is called for each workset,
where workset information can be passed in through the evaluate call:

\code
    std::vector<WorksetData> workset_data;
        .
        .
    // Initialize workset data 
        .
        .
    for (std::size_t i = 0; i < workset_data.size(); ++i) {
      field_manager.evaluateFields<MyTraits::Residual>(workset_data[i]);
	.
        .  
      // use evaluated fields 
        .
        .
    }
\endcode

Note that you do not have to use the workset idea.  You could just
pass in workset size equal to the number of local cells on the
processor or you could use a workset size of one cell and wrap the
evaluate call in a loop over the number of cells.  Be aware that this
can result in a possibly large performance hit.

Phalanx, in fact, does not restrict you to cell based iteration.  You can iterate over any entity type such as edge or face structures.    

<li><b>Consistent Evaluation</b>

Phalanx was written to perform consistent evaluations.  By consistent,
we mean that all dependencies of a field evaluation are current with
respect to the current degree of freedom values.  For example, suppose
we need to evaluate the the energy flux.  This has dependencies on the
density, diffusivity, and the temperature gradient.  Each of these
quantities in turn depends on the temperature.  So before the density,
diffusivity and temperature gradient are evaluated, the temperature
must be evaluated.  Before the energy flux can be evaluated, the
density, diffusivity, and temperature gradient must be evaluated.
Phalanx forces an ordered evaluation that updates fields in order to
maintain consistency of the dependency chain.  Without this, one might
end up with lagged values being used from a previous evaluate call.

<li><b>Scalar Type</b>

A scalar type, typically the template argument ScalarT in Phalanx
code, is the type of scalar used in an evaluation.  It is typically a
double or float, but can be special object types for embedded methods such as
sensitivity analysis.  For example, for sensitivity analysis, a double
scalar type is replaced with a foward automatic differentiation object
(FAD) or a reverse automatic differentaion object (RAD) to produce
sensitivity information.  Whatever type is used, the standard
mathematical operators are overloaded for the particular embedded
technology.  For an example of this, see the <a
href="http://trilinos.sandia.gov/packages/sacado">Sacado Automatic
Differentiation Library</a>. Some sample scalar types include:

<ul>
<li> float
<li> double
<li> Sacado::Fad::DFad<double> (for sensitivity analysis)
</ul>

<li><b>Algebraic Type</b>

An algebraic type is the type of objects that hold data.  It is usually a rank n tensor.  Simple examples include a scalar (rank-0 tensor), a vector (rank-1 tensor) or a matrix (rank-2 tensor).  It is not actually restircted to tensors, but can be any struct/class/data type that a user implements.  The algebraic type is a description of how data is stored but does NOT have a corresponding type in the Phalanx code.  It is a notion or idea we use to describe a data type without specifying the actual scalar type (See "Data Type" for more information).  These types are defined by the user.  In the example in the directory "phalanx/example/EnergyFlux", the user selects three algebraic types to represent scalars, vectors and tensors.  The scalar algebraic type is equivalent to the scalar type used in the evaluation.  The vector and tensor objects are objects templated on the scalar type:

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

Note that instead of using an algebraic type, most users can use a multidimensional array to meet their needs.

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

The DataLayout object is used to distinguish fields with the same name, but exist at different locations in the discretization structure (i.e. cell or face).  For example, supposed we have written an evaluator the computes the "Density" field for a set of points in the cell.  Now we want to evaluate the density at a different set of points in the cell.  We might have a "Density" field in a cell associated with a set of integration points (quadrature points in finite elements) and another field associated with the nodes (nodal basis points in finite elements).  We use the same field name (so we can reuse the same Evaluator), "Density", but use two different DataLayouts, one for integration points and one for nodal point.  Now a FieldTag comparison will differentiate the fields due to the different DataLayout.  

Additionally, the DataLayout contains the number of DataT objects associated with the field in the discretized structure.  This size() parameter is not needed to distinguish uniqueness, since the number of objects can be the same for different fields.  It is stored here for convenience when figuring out the size of field arrays.

<li><b>Field Tag</b>

The FieldTag is a description of a field.  It is templated on the data type, DataT.  It is used to identify a field stored in the field manager.  It contains a unique identifier (an stl std::string) and a pointer to a data layout object.  If two FieldTags are equal, the DataLayout, the data type, and the string name are exactly the same.

<li><b>Field</b>

A Field is a set of values for a particular quantity of interest that must be evaluated.  It is templated on the data type, DataT.  It consists of a FieldTag and a reference counted smart pointer (Teuchos::RCP<DataT>) to the data array where the values are stored.  

<li><b>Evaluator</b>

An Evaluator is an object that evaluates a set of Fields.  It contains two vectors of Fields, one set is the fields it will evaluate and one set is the fields it depends on for the evaluation.  For example to evaluate a density field that is a function of temperature and pressure, the field the evaluator will evaluate (the first field set) is a density field, and the set of fields it requires to meet its dependencies are the temperature and pressure fields (the second field set).  The evaluator is templated on the evaluation type and the traits class.  

<li><b>Evaluator Manager</b>

The main object that stores all Fields and Evaluators.  The evaluator manager (EM) sorts the evaluators and determines which evaluators to call and the order necessary to ensure consistency in the fields.  The EM also allocates the memory for storage of all fields so that if we want to force all fields to be in a contiguous block, we have that option available.  Users can write their own allocator for memory management.

<li><b>Multidimensional Array</b>

Instead of using the concept of an algebraic type that the user implements, it may be easier to use a multidimensional array.  For example, suppose we have ten cells, and in each cell there are four points (specifically quadrature points) where we want to store a 3x3 matrix for the stress tensor.  Using the concepts of algebraic types, we would use a user defined matrix in a Phalanx Field:

\code
  PHX::MDA::Layout<Cell,QP> layout(10,4);
  PHX::Field< MyTensor<double> > stress("Stress",layout);
\endcode

However, Phalanx also implements the idea of a multidimensional array with optional compile time checking on rank accessors.

\code
  PHX::MDA::Layout<Cell,QP,Dim,Dim> layout(10,4,3,3);
  PHX::MDField<double,Cell,QP,Dim,Dim> stress("Stress",layout);
\endcode

Here, the "Cell", "QP", and "Dim" objects are small structs that allow users to describe the ordinates associated with the multidimensional array.

The benefits of using the multidimensional array are that (1) checking of the rank accessor at either compile time or runtime (runtime checking is only enabled for debug builds for efficiency) prevent coding errors and (2) the documentation of the ordinals is built into the code - no relying on comments that go out of date.

The EnergyFlux example is reimplemented using the multidimensional array instead of the algebric types in the directory Trilinos/packages/phalanx/example/MultiDimensionalArray/.

  Our recomendation is to use the multidimensional array version as future codes plan to use this object.  The PHX::MDField is fully compatible with Intrepid whereas the PHX::Field is not.  We keep the PHX::Field object around for performance measurements since it directly accesses the Teuchos::ArrayRCP object where the PHX::MDField uses a shards::Array object.  If the compiler optimizes correctly, there should be no difference in performance.  By testing both the Field and MDField, we can test compiler optimization.   

</ul>

\section performance Performance

Some recomendations for efficient code:
<ul>

<li> <b>Use worksets</b> This may eliminate cache misses. 

<li> <b>Enable compiler optimization:</b>  The Field and MDField classes use inlined bracket operators for data access.  That means if you build without optimization some compilers will produce very slow code.  To see if your compiler is optimizing away the bracket operator overhead, run the test found in the directory "phalanx/test/Performance/BracketOperator".  If the timings are the same between a raw pointer array and the Field and MDField classes, your compiler has removed the overhead.   For example, on gnu g++ 4.2.4, compiling with -O0 shows approximately 2x overhead for bracket accessors, while -O2 shows about a 10% overhead, while -O3 completely removes the overhead.

       <li> <b>Algebraic Types:</b> Implementing your own algebraic types, while convenient for users, can introduce overhead as opposed to using raw arrays or the multidimensional array.  The tests found in "phalanx/test/Performance/AlgebraicTypes" demonstrate some of this overhead.  Expression templates should remove this overhead, but in our experience, this seems to be very compiler and implementation dependent.  If you build in tvmet support, you can compare expression templates, our "dumb" implementation for vectors and matrices, and our multidimensional array against raw array access.  Our testing on gnu compilers shows an overhead of about 20-25% when using "dumb" objects and the expresion templates as opposed to raw arrays.  We recommend using raw arrays with the multi-dimensional array (MDField) for fastest runtimes.

       <li> <b>Use Contiguous Allocator:</b> Phalanx has two allocators, one that uses the "new" command for each separate field, and one that allocates a single contiguous array for ALL fields.  If cache performance is the limiting factor, the contiguous allocator could have a big effect on performance.  Additionally, alignment issues can play a part in the allocators depending on how you implement your algrbraic types.  Our ContiguousAllocator allows users to choose the alignment based on a template parameter.  Typically, this is double.

<li> <b>Limit the number of Fields and/or Evaluators:</b>  The more evaluators used in your code, the more the loop sturcutre is broken up.  You go from a single loop to a bunch of small loops.  This can have an effect on the overall performance.  Users should also be judicious on choosing Fields.  Only select Fields as a place to introduce variability in your models or for reuse.  For example, if you require density in a single place in the code and there is only one single model that won't change, do not make it a field.  Just evaluate it once where needed.  But if you need density in multiple providers or you want to swap models at runtime, then it should be a field.  This prevents having to recompute the model or recompile your code to switch models.  This is usually not an issue as long as the amount of work in each evaluator is larger than vtable lookup to make the evaluate call.  In our experience, we have never observed this behaviour.  We point it our here just in case.

<li> <b>Slow compilation times:</b>  As the number of Evaluators in a code grows, the compilation times can become very long.  Making even a minor change can result in the code recompiling all Evaluator code.  Therefore, we recommend using explicit template instantiation.  An example can be found in Trilinos/packages/phalanx/example/FEM_Nonlinear.  All evaluators use explicit template instantiation if phalanx is built with explicit template instation enabled (in the cmake build system, use the flag -D Phalanx__EXPLICIT_TEMPLATE_INSTANTIATION=ON).

</ul>

\section user_guide_mdarray_domain_model Multi-Dimensional Array Domain Model

The multidimensional array was designed for interoperability between a number of Trilinos packages including shards, phalanx and intrepid.  These codes each have their own implementations but follow a basic set of requirements so that functions templated on the array type can use any of the multidimensional array implementations.  The required functions are:

<ul>
<li> ScalarT& operator[] - the bracket operator
<li> ScalarT& operator(1,2,3,...,N) - accessor operators for each rank N array.
<li> size_type rank() - integer number of ordinates
<li> size_type dimension(size_type ordinate) - size of ordinate
<li> size_type size() - total size of array
</ul>

Phalanx implements the multidimensional array in the PHX::MDField class and supports arrays with up to 8 ranks (one more than Fortran arrays support).  More information can be found in the shards library.

\section user_guide_step1 Step 1: Configuring, Building, and installing Phalanx

\subsection ug_step1_general A. General Library Requirements
Phalanx is distributed as a package in the <a href="http://trilinos.sandia.gov">Trilinos Framework</a>.  It can be enabled as part of a trilinos build with the configure option "-D Trilinos_ENABLE_Phalanx=ON".  Phalanx currently has direct dependencies on the following third party libraries:

- <b>Requires</b> the <a href="http://trilinos.sandia.gov/packages/teuchos">Teuchos</a> utilities library, part of the <a href="http://trilinos.sandia.gov/">Trilinos Framework</a>.  This will automatically be enabled when you enable the phalanx library.
 
 - <b>Requires</b> the <a href="http://trilinos.sandia.gov/packages/sacado">Sacado Automatic Differentiation Library</a>, part of the <a href="http://trilinos.sandia.gov/">Trilinos Framework</a>.  This will automatically be enabled when you enable the phalanx library.

 - <b>Requires</b> the <a href="http://www.boost.org">Boost Template Metaprogramming (MPL) Library</a>.  This is a third party library (TPL) that must be installed on your machine.  You must enable the TPL and point to the path to the Boost library during Trilinos configuration.  An example configuration file can be found in Trilinos/packages/phalanx/maintenance/reconfigure.linux.mpi.cmake.

\subsection ug_step1_performance B. Performance Example Requirements

 - <b>Optional:</b> Some performance tests run comparisons against <a href="http://tvmet.sourceforge.net/">TVMET: Tiny Vector Matrix library using Expression Templates</a>.  This is to get a feel for how our "dumb" vector matrix objects perform compared to expression templates.  You must enable the tvmet TPL add the path to the TVMET library during Trilinos configuration.  An example configuration file can be found in Trilinos/packages/phalanx/maintenance/reconfigure.linux.mpi.cmake.

TVMET is optional and hidden behind an ifdef.  The performance tests will be built regardless of whether tvmet is enabled/disabled. 

\subsection ug_step1_fem C. Nonlinear Finite Element Example Requirements

To build the example problem in "phalanx/example/FEM_Nonlinear", distributed vector, distributed sparse matrix, and corresponding linear solvers are required.  The following <b>Trilinos</b> packages need to be enabled to build the FEM_Nonlinear example:

 - <b>Requires:</b> the <a href="http://trilinos.sandia.gov/packages/epetra">Epetra Library</a>, supplies the linear algebra data structures for spares matrices.  Can be enabled during the Trilinos configure with the flag "-D Trilinos_ENABLE_Epetra=ON". 

 - <b>Requires:</b> the <a href="http://trilinos.sandia.gov/packages/ifpack">Ifpack Library</a>, supplies incomplete factorization preconditioners.  Can be enabled during the Trilinos configure with the flag "-D Trilinos_ENABLE_Ifpack=ON". 

 -  <b>Requires:</b> the <a href="http://trilinos.sandia.gov/packages/belos">Belos Library</a>, supplies block GMRES iterative linear solver.  Can be enabled during the Trilinos configure with the flag "-D Trilinos_ENABLE_Belos=ON".  

This example will be disabled if the above packages are not enabled.

\subsection da D. Install Boost

You must have boost installed on your system.  As we only require the MPL library headers, you only need to untar the source code and should not even need to run the configure script.

\subsection db E. Configure Trilinos/Phalanx
The general instructions for building trilinos can be found at <a href="http://trilinos.sandia.gov/documentation.html">Trilinos Documentation Page</a>.  Of particular importance are the Overview, User Guide, and Tutorial documents.  At a minimum you must enable the Teuchos, Sacado, and Phalanx packages.  An example configure script is:

\verbinclude reconfigure.linux

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
      
- \b SetupData - A user defined type to be passed in to the postRegistrationSetup() call.  Allows users to pass in arbitrary data to the evaluators during setup.

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

The typedefs RealType and FadType are done only for convenience.  They are not actually required but cut down on the typing.  Only the EvalTypes typedef is required in the code above.

\subsection td3 C. EvaltoDataMap
Next we need to link the data types to the evaluation type.  Note that one could use the same data type in multiple evaluation types:

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
  Users can pass their own data to the postRegistrationSetup(), evaluateFields(), preEvaluate() and postEvaluate() methods of the PHX::FiledManager class.  In this example, the user passes in a struct that they have written called MyEvalData.  This contains information about the cell workset.  The user is not required to write their own object.  They could just pass in a null pointer if they don't need auxiliary information passed into the routine.  This is demonstrated in the SetupSetup, PreEvalData, and PostEvalData.  A void* is set for the data member. 
\code
     .
     .
     .
    // ******************************************************************
    // *** User Defined Object Passed in for Evaluation Method
    // ******************************************************************
    typedef void* SetupData;
    typedef const MyEvalData& EvalData;
    typedef void* PreEvalData;
    typedef void* PostEvalData;

  };
\endcode

\section user_guide_step4 Step 4: Specialize the PHX::TypeString Object
For debugging information, Phalanx makes a forward declaration of the PHX::TypeString object.  This must be specialized for each evaluation type and each data type so that if there is a run-time error, phalanx can report detailed information on the problem.  We could have used the typeinfo from the stl, but the name() method is not demangled on every platform, so it can make debugging a challenge.  The specialized classes can go into their own file or can be added to the traits file above depending on how you use the Traits class.  During linking, if the compiler complains about multiple defninitions of your specialized traits classes, separate the traits implementation into their own .cpp file.

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
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
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
postRegistrationSetup(typename Traits::SetupData d,
		      PHX::FieldManager<Traits>& vm)
{
  this->utils.setFieldData(density,vm);
  this->utils.setFieldData(temp,vm);

  data_layout_size = density.fieldTag().dataLayout().size();
}

// **********************************************************************
template<typename EvalT, typename Traits>
void Density<EvalT, Traits>::evaluateFields(typename Traits::EvalData d)
{ 
  std::size_t size = d.num_cells * data_layout_size;
  
  for (std::size_t i = 0; i < size; ++i)
    density[i] =  temp[i] * temp[i];
}
\endcode

The constructor pulls out data from the parameter list to set the correct data layout.  Additionally, it tells the FieldManager what fields it will evaluate and what fields it requires/depends on to perform the evaluation.

The postRegistrationSetup method gets pointers from the FieldManager to the array for storing data for each particular field.

  Writing evaluators can be tedious.  We have invested much time in minimizing the amount of code a user writes for a new evaluator.  Our experience is that you can literally have hundreds of evaluators.  So we have added macros to hide the boilerplate code in each evaluator.  Not only does this streamline/condense the code, but it also hides much of the templating.  So if your user base is uncomfortable with C++ templates, the macro definitions could be very helpful.  The definitions are found in the file Phalanx_Evaluator_Macros.hpp.  The same evaluator shown above is now implemented using the macro definitions:

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
PHX_POST_REGISTRATION_SETUP(Density,data,fm)
{
  this->utils.setFieldData(density,fm);
  this->utils.setFieldData(temp,fm);

  data_layout_size = density.fieldTag().dataLayout().size();
}

//**********************************************************************
PHX_EVALUATE_FIELDS(Density,d)
{ 
  std::size_t size = d.num_cells * data_layout_size;
  
  for (std::size_t i = 0; i < size; ++i)
    density[i] =  temp[i] * temp[i];
}
\endcode

The evaluators for the example problem in "phalanx/example/EnergyFlux" have been rewritten using the macro definitions and can be found in the directory "phalanx/test/Utilities/evaluators".  

Finally, since writing even the above code contains much boilerplate, we have written a python script that will generate the above skeleton files for you.  All you need provide is the class name and the filename.  The script is called phalanx_create_evaluator.py and can be found in the "maintenance" directory.  A "make install" will place the script in the "bin" directory.  To generate a skeleton for the above function, you would execute the following command at the prompt:

\code
> ./phalanx_create_evaluator.py -c -n Density Evaluator_Density
\endcode

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

      { // Constant Thermal Conductivity
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_constant;
	p->set<int>("Type", type);
	p->set<string>("Name", "Thermal Conductivity");
	p->set<double>("Value", 2.0);
	p->set< RCP<DataLayout> >("Data Layout", qp);
	evaluators_to_build["Thermal Conductivity"] = p;
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

You are free to request fields to evaluate until the time you call the method postRegistrationSetup() on the field manager.  Once this is called, you can not request any more fields. 

\subsection fmd3 C. Call FieldManager::postRegistrationSetup()

Once the evaluators are registered with the FieldManager and it knows which field it needs to provide, call the postRegistrationSetup() method on the FieldManager.  This method requires that the user specify the maximum number of cells for each evaluation - the workset size.  This number should be selected so that all fields can fit in the cache of the processor if possible.  

\code
      // Assume we have 102 cells on processor and can fit 20 cells in cache
      const std::size_t num_local_cells = 102;
      const std::size_t workset_size = 20;

      fm.postRegistrationSetup(workset_size);
\endcode

The postRegistrationSetup() method causes the following actions to take place in the FieldManager:

<ol>
<li> Based on the requested fields in \ref fmd2, the FieldManager will trace through the evaluators to determine which evaluators to call and the order in which they need to be called to achieve a consistent evaluation. Not all evaluators that are registered will be used.  They will only be called to satisfy dependencies of the required fields.
<li> Once the dependency chain is known, we can pull together a flat list of all fields that will be used.  Now the FieldManager will allocate memory to store all fields.  It will use the Allocator object from the users traits class to accomplish this.  By unifying the allocation into a single object, we can force all fields to be allocated in a single contiguous block of memory if desired.
<li> Once the memory for field data is allocated, we must set the pointer to that memory block inside each Field or MDField object in each evaluator.  The FieldManager does this by calling the postRegistrationSetup() method on each evaluator that is required for the computation. 
</ol>

Note: you can no longer register evaluators or request fields to evaluate once postRegistrationSetup() method is called.

\subsection fmd6 D. Setup Worksets

If the user plans to use worksets, these objects are typically passed in through the member of the evaluate call as defined in your traits class.  In our example, we chose to pass a user defined class called a MyWorkset in to the evaluate call.  Since we plan to use worksets, each object of MyWorkset contains information about the workset.  Here is the MyWorkset implementation: 

\code
#ifndef PHX_EXAMPLE_MY_WORKSET_HPP
#define PHX_EXAMPLE_MY_WORKSET_HPP

#include "Phalanx_ConfigDefs.hpp" // for std::vector
#include "Cell.hpp"

struct MyWorkset {
  
  std::size_t local_offset;

  std::size_t num_cells;
  
  std::vector<MyCell>::iterator begin;

  std::vector<MyCell>::iterator end;

};

#endif
\endcode

The user has written a MyCell class that contains data for each specific local cell.  MyWorkset contains iterators to the beginning and end of the chunk of cells for this particular workset.  The local_offset is the starting index into the local cell array.  The num_cells is the number of cells in the workset.   The MyWorkset objects are created one for each workset via the following code:

\code
      // Create Workset information: Cells and EvalData objects
      std::vector<MyCell> cells(num_local_cells);
      for (std::size_t i = 0; i < cells.size(); ++i)
	cells[i].setLocalIndex(i);
      std::vector<MyWorkset> worksets;

      std::vector<MyCell>::iterator cell_it = cells.begin();
      std::size_t count = 0;
      MyWorkset w;
      w.local_offset = cell_it->localIndex();
      w.begin = cell_it;
      for (; cell_it != cells.end(); ++cell_it) {
	++count;
	std::vector<MyCell>::iterator next = cell_it;
	++next;
	
	if ( count == workset_size || next == cells.end()) {
	  w.end = next;
	  w.num_cells = count;
	  worksets.push_back(w);
	  count = 0;

	  if (next != cells.end()) {
	    w.local_offset = next->localIndex();
	    w.begin = next;
	  }
	}
      }
\endcode

Note that you do not have to use the workset idea.  You could just pass in workset size equal to the number of local cells on the processor or you could use a workset size of one cell and wrap the evaluate call in a loop over the number of cells.  Be aware that this can result in a possible performance hit.

\subsection fmd4 E. Call evaluate()

Finally, users can call the evlauate routines and the pre/post evaluate routines if required.

\code
      fm.preEvaluate<MyTraits::Residual>(NULL);

      // Process all local cells on processor by looping over worksets
      for (std::size_t i = 0; i < worksets.size(); ++i) {

	fm.evaluateFields<MyTraits::Residual>(worksets[i]);
	  
	// Use workset values
                .
                .
                .
      }

      fm.postEvaluate<MyTraits::Residual>(NULL);
\endcode

\subsection fmd5 F. Accessing Data

Accessing field data is achieved as follows:

\code
      Field< MyVector<double> > ef("Energy_Flux", qp);

      fm.getFieldData<MyVector<double>,MyTraits::Residual>(ef);
\endcode

You do not need to use the Field objects to access field data.  You can get the reference counted smart pointer to the data array directly by using the field tag:

\code
      RCP<DataLayout> qp = rcp(new FlatLayout("QP", 4));
      PHX::FieldTag<double> s("Nonlinear Source", qp);
      Teuchos::ArrayRCP<double> source_values;
      fm.getFieldData<double,MyTraits::Residual>(s,source_values);
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

\ref faq5

\ref faq6

\ref faq7

\ref faq8

\ref faq9

\section faq1 1. Why name it Phalanx?  
The phalanx was one of the most dominant military formations of the Greek armies in the classical period.  It was a strictly ordered formation.  The Phalanx software library figures out ordered dependencies of field evaluators.   A second more obscure reference relates to the US Navy.  The Phalanx software package was designed to provide nonlinear functionality for the <a href="http://trilinos.sandia.gov/packages/intrepid">Intrepid</a> itegration library.  Intrepid was the name of an aircraft carrier during in WW II.  Modern US aircraft carriers are protected by a Close-In Weapons System (CIWS) named the Phalanx.  Finally, the PI of this project is an avid strategy warfare gamer and leans towards military references.

\section faq2 2. Why does the evaluation type exist?
This allows for multiple evaluations that used the same Scalar type, but different evaluators specific to the CalculationType.  For example, evaluating the Jacobian (derivative of the residual with respect to the degrees of freedom) and the parameter sensitivities (derivative of the residual with respect to the parameters) both use the same Sacado data type.  They both want derivatives but with respect to a different set of entities.

\section faq3 3. Why do evaluation types have arbitrary data types.  Shouldn't all fields of an evaluation type require the same scalar type?
In 99% of the use cases, the answer is yes!  By changing the scalar type for a specific type of computation, you can break the dependency chain and get incorrect answers.  For example, if computing a Jacobian, you would use a FAD scalar type instead of a double scalar type for the evaluation.  By mixing FAD and double fields in an evaluation type, you will lose derivative components and get incorrect sensitivities.  But there are some use cases where one might want to use different data types.  If the user knows that certain evaluation routines have no derivative components, they can bypass the AD chain by switching to doubles for that part of the computation.  The is an expert only user mode, but we allow for it.

\section faq4 4. Why are the evaluators templated on the evaluation type?
Each evaluator must be templated on the evaluation type to be able to (1) allow for the expert use mode described in FAQ \ref faq3; (2) use the automated factory class (PHX::EvaluatorFactory) to build the evaluator for each evalution type.  

\section faq5 5. Why is the MDALayout object not templated on the ArrayOrder type?
This introduces many complications in the comparison operator (operator==) during a comparison of a reverse and natural data layout with the same (but reversed index types).  Should these objects be the same?  Probably, but this might introduce some confusion.  So in an effort to keep things as clean as possible, there is no differentiation based on forward or reverse (Fortran) array ordering.  If we were to introduce this, then we would have to add another method to the DataLayout base class to allow for checking the ordering.  If mixing natural and reverse ordering becomes common, we may revisit this decision.

\section faq6 6. My code seems to be running slow.  How can I speed it up?
See the section on \ref performance in the Users Guide that gives recomendations on how to maximize the efficiency of your code.

\section faq7 7. Compilation take a long time when minor changes are made to an evaluator.  How can I speed this up?
Explicit template instantiation can be used to speed up evaluator development builds.  See the section on \ref performance in the Users Guide.  

\section faq8 8. Valgrind gives a memory leak error in all evaluators when using the PHX::ContiguousAllocator with Sacado::DFad based fields.  What is the problem?
The contiguous allocator is allocating space for all variables in an evaluation type in a single contiguous chunk of memory.  We allow for multiple scalar types in an evaluation type so the array is allocated using the type char and then based on scalar type sizes (using sizeof method) along with proper alignment, pointers to blocks of memory are handed to the fields using a reinterpret_cast.  DFad uses new to allocate the derivative arrays, so this memory is not contiguous.  The problem is that the field array memory is deleted as a char instead of it's various data types, so all DFad derivative arrays are leaked during the destruction.  There are two ways to alleviate this: (1) write a special Allocator to call the destructor of the field type for each field.  (2) Have Sacado implement a version of DFad that allows the user to allocate the DFad array and pass in the pointer during construction.  This will also require a new allocator.  Choice (2) will be handed in Phalanx in a later release.  Our current suggestion is that one should use DFad only if using the PHX::NewAllocator and use SFad if using the PHX::ContiguousAllocator. 

\section faq9 9. The evaluators contain boilerplate.  Is there an easy way to automate creating the evaluators?
YES!  There is a python script that will automatically generate the skeleton for evaluators.  It is found in the "maintenance" directory of phalanx and is called "phalanx_create_evaluator.py".  It will be install in the "bin" directory when trilinos is installed.  Example use can be found at the end of section \ref user_guide_step5 of the \ref user_guide.

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

  \todo Nothing pending!

*/

/* ************************************************************************ */
/* ************************************************************************ */

#endif
