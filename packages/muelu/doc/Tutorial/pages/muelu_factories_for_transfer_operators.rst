======================================
MueLu factories for transfer operators
======================================

For this example,
we reuse the **Recirc2D** example as introduced in :ref:`multigrid_for_non_symmetric_problems/test_example`.
The resulting linear systems are (slightly) non-symmetric
and classical smoothed aggregation methods may be able to solve the problem but are not optimal in sense of convergence.

Multigrid setup phase - algorithmic design
==========================================

Smoothed aggregation based algebraic multigrid methods originally have not been designed for non-symmetric linear systems.
Inappropriately smoothed transfer operators may significantly deteriorate the convergence rate or even break convergence completely.

Unsmoothed transfer operators
-----------------------------

.. warning::

	Insert and link missing figures

Before we introduce smoothed aggregation methods for non-symmetric linear systems,
we first go back one step and demonstrate how to use non-smoothed transfer operators which are eligible for non-symmetric linear systems.
Figure :ref:`muelu_factories_for_transfer_operators/figure_simpledesignnonsmoothed` gives a simplified example how to build the coarse level matrix :math:`A_c`
using the fine level matrix :math:`A` only.
First, we "somehow" build aggregates using the information of the fine level matrix :math:`A`.
The aggregates are then used to build the tentative non-smoothed prolongation operator.
The restrictor is just the transpose of the (tentative) prolongator and finally the coarse level matrix :math:`A_c` is calculated by the triple product :math:`A_c=RAP`.

In Figure :ref:`muelu_factories_for_transfer_operators/figure_simpledesignaamg`,
the **SaPFactory** has been added after the **TentativePFactory**.
Therein the non-smoothed transfer operator from the **TentativePFactory** is smoothed using information of the fine level matrix :math:`A`.
This transfer operator design is used per default when the user does not specify its own transfer operator design.
The default settings are optimal for symmetric positive definite systems.
However, they might be problematic for our non-symmetric problem.

Smoothed transfer operators for non-symmetric systems
-----------------------------------------------------

.. warning::

	Insert missing figures

In case of non-symmetric linear systems it is :math:`A\neq A^T`.
Therefore it is a bad idea just to use the transposed of the smoothed prolongation operator for the restrictor.
Let :math:`\widehat{P}` be the non-smoothed tentative prolongation operator.
Then the smoothed prolongation operator :math:`P` is built using

:math:`P = \bigl(I-\omega A\bigr) \widehat{P}`

with some reasonable smoothing parameter :math:`\omega>0`.
The standard restrictor is

:math:`R = P^T = \widehat{P}^T - \omega \widehat{P}^T A^T = \widehat{P}^T\bigl(I-\omega A^T\bigr).`

That is, the restrictor would be smoothed using the information of :math:`A^T`.
However, for non-symmetric systems we want to use the information of matrix :math:`A` for smoothing the restriction operator, too.
The restriction operator shall we built by the formula

:math:`R = P^T = \widehat{P}^T - \omega \widehat{P}^T A.`

This corresponds to apply the same smoothing strategy to the non-smoothed restriction operator :math:`\widehat{R}=\widehat{P}^T`,
which is applied to the (tentative) prolongation operator with using :math:`A^T` as input instead of matrix :math:`A`.
Figure :ref:`muelu_factories_for_transfer_operators/figure_simpledesignpamg` shows the changed factory design.
The dashed line denotes, that the same smoothing strategy is used than for the prolongation operator.
The concept is known as Petrov-Galerkin smoothed aggregation approach in the literature.
A more advanced transfer operator smoothing strategy for non-symmetric linear systems that is based on the Petrov-Galerkin approach is described in [1]_.
Another approach based on SchurComplement approximations can be found in [2]_.

**Insert missing figures here**

.. _muelu_factories_for_transfer_operators/figure_simpledesign:
.. _muelu_factories_for_transfer_operators/figure_simpledesignpgamg:
.. _muelu_factories_for_transfer_operators/figure_simpledesignaamg:
.. _muelu_factories_for_transfer_operators/figure_simpledesignnonsmoothed:

XML Interface
=============

Unsmoothed transfer operators
-----------------------------

To construct a multigrid hierarchy with unsmoothed transfer operators one can use the following XML file (stored in ../../../test/tutorial/s3a.xml)

.. literalinclude:: ../../../test/tutorial/s3a.xml
  :language: xml
	:caption:

Beside the **TentiativePFactory** which is responsible to generate the unsmoothed transfer operators we also introduce the **UncoupledAggregationFactory** with this example.
In the **Factories** section of the XML file, you find both an entry for the aggregation factory and the prolongation operator factory with its parameters.
In the **Hierarchy** section the defined factories are just put in into the multigrid setup algorithm.
That is, the factory with the name **UncoupledAggregationFact** is used to generate the **Aggregates**
and the **myTentativePFact** is responsible for generating both the (unsmoothed) prolongation operator P and the (coarse) near null space vectors **Nullspace**.

.. note::

	For some more details about the (hidden) **NullspaceFactory**,
	which is internally used to handle the null space information and the dependencies,
	the reader might refer to :ref:`multigrid_for_multiphysics/blocktransfersetup`.

.. note::

	One can also use the **Ptent** variable for registering a **TentativePFactory**.
	This makes the **TentativePFactory** somewhat special in its central role for generating an aggregation based multigrid hierarchy.
	MueLu is smart enough to understand that you want to use the near null space vectors generated by the factory registered as **Ptent** for setting up the transfer operators.

	That is, the following code would explicitly use the **TentativePFactory** object that is created as **myTentativePFact**.
	Since no factory is specified for the prolongation operator **P**,
	MueLu would decide to use a smoothed aggregation prolongation operator (represented by the **SaPFactory**),
	which correctly uses the factory for **Ptent** for the unsmoothed transfers with all its dependencies.

	.. code-block:: xml

		<ParameterList name="MueLu">
		<ParameterList name="Factories">
			<ParameterList name="myTentativePFact">
			<Parameter name="factory" type="string" value="TentativePFactory"/>
			</ParameterList>
		</ParameterList>

		<ParameterList name="Hierarchy">
			<ParameterList name="Levels">
			<Parameter name="Ptent" type="string" value="myTentativePFact"/>
			</ParameterList>
		</ParameterList>
		</ParameterList>

.. admonition:: Exercise 1

	Create a sublist in the **Factories** part of the XML file for the restriction operator factory.
	Use a **TransPFactory** which builds the transposed of **P** to generate the restriction operator **R**.
	Register your restriction factory in the **Hierarchy** section to generate the variable **R**.

Smoothed aggregation for non-symmetric problems
-----------------------------------------------

Next, let's try smoothed transfer operators for the non-symmetric linear system and compare the results of the transfer operator designs.
Take a look at the XML file (in **../../../test/tutorial/s3b.xml**).

.. literalinclude:: ../../../test/tutorial/s3b.xml
  :language: xml
	:caption:

The interesting part is the **Factories** section where several different factories for the restriction operator are defined

.. admonition:: Description

	* **myTentRestrictorFact:** just uses the transposed of the unsmoothed prolongator for restriction.
	* **mySymRestrictorFact:** uses the transposed of the smoothed prolongator for restriction.
	* **myNonsymRestrictorFact:** uses the special non-symmetric smoothing for the restriction operator (based on the **SaPFactory** smoothing factory).


.. note::

	The MueLu framework is very flexible and allows for arbitrary combinations of factories.
	However, be aware that the **TentativePFactory** cannot be used as input for the **GenericRFactory**.
	That is no problem since this combination not really makes sense.
	If you are using the **TentativePFactory** as your final prolongation operator you always have to use the **TransPFactory** for generating the restriction operators.


.. admonition:: Exercise 2

	Run the **Recirc2D** example with the different restriction operator strategies and compare the results for the iterative solver.
	What do you observe? What is the best choice for the transfer operators in the non-symmetric case?

.. admonition:: Exercise 3

	Change the **myProlongatorFact** from type **SaPFactory** to **PgPFactory**,
	which uses automatically calculated local damping factors instead of a global damping factor (with some user parameter **sa: damping factor**).
	Note that the **PgPFactory** might not accept the **sa: damping factor** parameter such that you have to comment it out (using **<!-- ... -->**).

.. admonition:: Exercise 4

	Try to set up a multigrid hierarchy with unsmoothed transfer operators for the transition from the finest level to level 1
	and then use smoothed aggregation for the coarser levels (starting from level 1).


Footnotes
=========
.. [1] Sala, M. and Tuminaro, R. S., A new Petrov-Galerkin Smoothed Aggregation Preconditionerfor nonsymmetric Linear Systems, SIAM J. Sci. Comput., 2008, 31, p. 143–166
.. [2] Wiesner, T. A., Tuminaro, R. S., Wall, W. A. and Gee, M. W., Multigrid transfers fornonsymmetric systems based on Schur complements and Galerkin projections., Numer.Linear Algebra Appl., 2013, doi: 10.1002/nla.1889
