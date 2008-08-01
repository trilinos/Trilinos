#ifndef PHX_DOXYGEN_DOCUMENTATION_HPP
#define PHX_DOXYGEN_DOCUMENTATION_HPP

/*!

\mainpage

*** WARNING: DOCUMENTATION IS UNDER CONSTRUCTION ***

This will be updated for Trilinos release 9.0

\section contents Contents

  - \ref introduction

  - \ref overview

  - \ref domain_design

  - \ref user_guide

  - \ref developer_guide

  - \ref faq

  - \ref history

  - \ref authors

\section introduction Introduction

Phalanx is a library to handle arbitrary function evaluation with complex nonlinear dependency chains for discretized partial differential equation systems.  It provides a flexible and efficient mechanism for switching dependencies and computing sensitivities either analytically or via automatic differentiation.  It can be used with any cell-based discretization technique including finite element, finite volume and finite difference.

Phalanx is used to break down complex dependencies into manageable algorithmic blocks.  It was written to address a variety of difficulties encountered when writing a general PDE code that allows for flexible switching of dependencies at runtime.  

\section overview Overview

A simple example (found in phalanx/example/energyFlux) is the construction of a Fourier energy flux for the heat equation.  Suppose that we want to solve the heat equation using finite elements:

\f[
  \nabla \cdot (\mathbf{q}) + s = 0
\f] 

where \f$\mathbf{q}\f$ is the heat flux, \f$ \mathbf{q} = -\rho C_p \nabla T \f$  and \f$s\f$ is a nonlinear source term.  Integrating by parts and using the test function \f$ \phi \f$ results in the following weak form equation:

\f[
  - \int_{\Omega} \nabla \phi \cdot \mathbf{q} d\Omega 
  + \int_{\Gamma} \phi \mathbf{n} \cdot \mathbf{q} 
  + \int_{\Omega} \phi s d\Omega = 0 
\f]




\f[
  \mathbf{q} = -\rho C_p \nabla T
\f]

<h3>Advantages of using Phalanx:</h3>

 - Increased flexibility because each simpler piece becomes an extension point that can be swapped out with different implementations.

 - Easier to implement new code because each piece is simpler, more focused and easier to test in isolation.  

 - Ensures that complex dependency chains are evaluated correctly.  

 - Determine the order of variable evaluation so that the correct
 evaluation is performed: For example, to compute the density at the
 quadrature points, it could be a constant, it could be a function of
 temperature or it could be a function of temperature and pressure.  A
 factory creates the density provider and from this, the manager
 automatically figures out the correct evaluation order - whether T
 and/or P will have to be evaluated before the density.  This is an
 easy example, but there are much more complex dependency chains.

<h3>Disadvantages of using Phalanx:</h3>

 - A potential performance loss due to fragmentation of the over-all algorithm (e.g., several small loops instead of one large loop) 

 - A potential loss of visibility of the original, composite problem (since the code is scattered into multiple places).  

Managing these trade-offs can result in application code that both performs well and supports rapid development and extensibility.  

\section user_guide User's Guide

Phalanx is distributed as a package in the <a href="http://trilinos.sandia.gov">Trilinos Framework</a>.  It can be enabled as part of a trilinos build with the configure option "--enable-phalanx".  Phalanx currently has direct dependencies on the following third party libraries:

 - Requires the <a href="http://trilinos.sandia.gov/packages/teuchos">Teuchos</a> utilities library, part of the <a href="http://trilinos.sandia.gov/">Trilinos Framework</a>.
 
 - Requires the <a href="http://trilinos.sandia.gov/packages/sacado">Sacado Automatic Differentiation Library</a>, part of the <a href="http://trilinos.sandia.gov/">Trilinos Framework</a>.

 - Requires the <a href="http://www.boost.org">Boost Template Metaprogramming (MPL) Library</a>

 - Optional: Some performance examples runs comparisons against <a href="http://tvmet.sourceforge.net/">TVMET: Tiny Vector Matrix library using Expression Templates</a>.  This is to get a feel for how our "dumb" vector matrix objects perform compared to expression templates.  Use the configure option "--with-tvmet" to enable the tvmet functionality in the performance tests.

We plan to remove the Boost dependencies in a future release.  

\section domain_design Domain Design Model

<h3>A. Concepts</h3>

<ul>
<li><h3>Cell</h3>
Partial differential equations are solved in a domain.  This domain is discretized into cells (also called elements for the finite element method).  This library assumes that the block of cells being iterated over is of the same type!  If different evaluators (i.e. different material properties) are required in different blocks of cells, a new FieldMangager must be used to switch material properties.  This is required for efficiency.  A contiguous block of memory can be used for all fields allowing for very fast evaluation.

<li><h3>Scalar Type</h3>
A scalar type, typically the template argument ScalarT in Phalanx code, is the default type of scalar used in an evaluation.  It is typically a double, but can be special object types such as a foward automatic differentiation object (FAD) or a reverse automatic differentaion object when used to produce sensitivity information.

<li><h3>Algebraic Type</h3>
An algebraic type is the type of objects that hold data.  It is a rank n tensor.  For example it can be a scalar (rank-0 tensor), a vector (rank 1 tensor) or a matrix (rank 2 tensor).  It is not actually restircted to tensors, but can be any struct/class.  The only requirement is that it be templated on the scalar type.

<li><h3>Data Type</h3>
A data type, typically the template argument DataT in Phalanx code, is an actual type used for storing fields.  It is the combination of a specific scalar type and an algebraic type.  For example it can be a vector object of doubles: Vector<double>, or it could be a tensor object of FAD type: Tensor<Scadao::FAD<double> >, or it could just be a double.

<li><h3>Evaluation Type</h3>
The evaluation type, typically the template argument EvalT in Phalanx code, defines a unique type of evaluation to perform.  An EvaluationContainer is allocated for each conputation type specified in the traits class.  Examples include a residual type, a Jacobian type, and a parameter sensitivity type.  The EvaluationContainer is associated with at least one scalar type and possibly more.  The scalar type usually determines what is being evaluated. For example, to evaluate the equation residuals, the scalar type is usually a double.  To evaluate a Jacobian, the scalar type could be a forward automatic differentiation object, FAD<double>.  By using an evaluation type, the same scalar type can be used for different evaluation types.  For example computing the Jacobian and computing parameter sensitivities both could use the FAD scalar type.

<li><h3>Storage</h3>
A DataContainer object stores all fields of a particular data type.  Each evaluation type holds a vector of containers, one container for each vaid data type.

<li><h3>Data Layout</h3>
The DataLayout object is used to define a unique entity in a cell.  It's main purpose is to distinguish fields that have the same name, but exist on a different entity in the cell.  For example, supposed we have a field named "density".  We might want to store values in a cell associated with a set of integration points (quadrature points in finite elements) and another set of values associated with the nodes (nodal basis points in finite elements).  It additionally contains the number of entities in the cell.  This is not needed to distinguish uniqueness, since the number of entities can be the same for different entities.  It is stored here for convenience.

<li><h3>Field Tag</h3>
The FieldTag is a description of a field.  It is used to describe every field that is stored in the field manager.  It contains a unique identifier and a data layout.  It DOES NOT distinguish the actual data type for the field.  The same field tag can be passed to all evaluation type containers and each container can decide what data type it is associated with.  This leaves flexibility so that each evaluation type can change the data type for efficiency.  For example, one evaluation might want a field to be a FAD type, while another evaluation might want to change a field to a double type to eliminate wasted flops. 

<li><h3>Field</h3>
A Field is a set of values for a particular quantity of interest that must be evaluated.  It is templated on the data type, DataT.  It consists of a FieldTag and a reference counted smart pointer (Teuchos::RCP<DataT>) to the data array where the values are stored.  

<li><h3>Evaluator</h3>
An Evaluator is an object that evaluates a set of Fields.  It contains two sets of Fields, one set are the fields it will evaluate and one set are the fields it depends on for the evaluation.  For example to evaluate a density field that is a function of temperature and pressure, the field the evaluator will evaluate (the first field set) is a density field, and the set of field it requires to meet its dependencies are the temperature and pressure fields (the second field set).  The evaluator is templated on the evaluation type and the traits class.

<li><h3>Evaluator Manager</h3>
The main object that stores all fields and evauation routines. Users can get access to all fields through this object. 

</ul>

<h3>B. Design</h3>

Roger - these need to be updated!

- The evaluate call is fixed for each scalar type.  Each scalar type figures out it's own dependency chain.  However you can trick the evaluate routine into using mixed scalar types based on writing the evaluator to use a hard-coded scalar type.  This is nasty and can have problems if the provider is only registered for one scalar type.  This will need to be addressed in the future.  The problem is that the variable is independent of the ScalarType.

- Previous design built a manager for each scalar type.  This is too restrictive, so an additional layer was added called a computation type.  This allows for multiple evaluations that used the same Scalar type, but different evaluators specific to the CalculationType.

- The fields of a Computation type should be constructed on the same scalar type.  If you are doing a Jacobian evaluation using the scalar type FAD to get sensitivities, changing some fields to double will break the dependency chain.  There are, however, compelling use cases where one might want ot break that chain and drop sensitivites with respect to certain components.  Therefore, a calculation type cal have multiple scalar types.  This is a departure from the original design in Charon that enforced each calculation type to be a unique scalar type.

- Each evaluator must be templated on the calculation type to be able to use teh automated factory class to build the evaluator for each calculation type.  Each evaluator must be constructed with a parameter list as the single argument to be able to use the automated factory class.

\section developer_guide Developer's Guide

\section faq Frequently Asked Questions

  - Why name it Phalanx?  The phalanx was one of the most dominant military formations of the Greek armies in the classical period.  It was a strictly ordered formation.  The Phalanx software library figures out ordered dependencies of field evaluators.   A second more obscure reference relates to the US Navy.  The Phalanx software package was designed to provide nonlinear functionality for the <a href="http://trilinos.sandia.gov/packages/intrepid">Intrepid</a> itegration library.  Intrepid was the name of an aircraft carrier during in WW II.  Modern US aircraft carriers are protected by a Close-In Weapons System (CIWS) named the Phalanx.  Finally, the PI of this project is an avid strategy warfare gamer and leans towards military references.

  - Why is the scalar type not embedded in the FieldTag?  We designed the FieldTag to be scalar type independent.  The reason is that we could then use FieldTags as arguments for Evaluator objects of ANY scalar type.  We have a factory that can automatically build Evaluators for all scalar types.  This automation requires constructor arguments that are not dependent on the scalar type.  

\section history History

Phalanx grew out of the Variable Manger in the Charon code and the Expression Manager in the Aria code.  It is an attempt at merging the two capabilities and is slated to provide nonlinear function evaluation to the Intrepid discretiztion package.

\section authors Authors

  - Roger Pawlowski (PI), SNL 01414

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
