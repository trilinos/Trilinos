=================
Advanced concepts
=================

As already mentioned in the beginning,
MueLu is designed as a multigrid framework, and, even though initiated as an aggregation-based algebraic multigrid method,
it can also be used for other kinds of coarsening methods.
In this chapter we demonstrate the combination of a semi-coarsening method with an aggregation-based coarsening on the coarser levels.
Semi-coarsening is combined with a line-smoothing method which then changes to a point-relaxation smoothing once no further semi-coarsening is possible.
In both cases, the semi-coarsening and the line-smoothing, the key element here is the dynamic switch from one to the another coarsening or smoothing strategy during runtime.

Semi-coarsening
===============

Basic idea
----------

Assuming that you have a 3D problem which is based on an extruded 2D mesh,
semi-coarsening might be an interesting option.
That is, on the finer levels we apply a semi-coarsening transfer operator which basically reduces the problem to a pseudo 2D problem
which then is handled by any other type of (smoothed aggregation based) transfer operator the usual way.

Factory layout without rebalancing
----------------------------------

The semi-coarsening is provided by the **SemiCoarsenPFactory** for generating the semi-coarsening transfer operators
in combination with the **LineDetectionFactory** which performs the line detection,
i.e. it searches for the vertical node lines along the extrusion axis of the 2D mesh.
Once all mesh layers are reduced to one by the **SemiCoarsenPFactory**,
we switch to the aggregation-based standard coarsening process.
There is a **TogglePFactory** which allows to switch back and forth between two different transfer operator strategies
(such as semi-coarsening and standard aggregation-based coarsening).
In principle, any combination of two different transfer operator strategies is allowed.
However, the current implementation only contains decision criteria to switch between semi-coarsening and aggregation-based types of coarsening.

Figure :ref:`advanced_concepts/figure_rebalancedtoggledesign` shows a typical factory layout for the combination of semi-coarsening
with a standard (non-smoothed) aggregation-based coarsening with repartitioning enabled.
First, one can see the **TogglePFactory** which has knowledge about the two different transfer operator branches
and makes a decision which transfer operator is used and provided to the restriction operator factory and the **RAPFactory**.
Depending on the line detection algorithm in **LineDetectionFactory**,
one might also need some information about the fine level mesh (such as the fine level coordinates).
The user only has to provide the mesh on the finest level.
On the coarser levels the **LineDetectionFactory** uses a standard ordering of the degrees of freedom which corresponds to a vertical node ordering.

.. _advanced_concepts/figure_rebalancedtoggledesign:

.. warning::

    Insert missing figure here

Factory layout with rebalancing
-------------------------------
Figure :ref:`advanced_concepts/figure_rebalancedtogglerebalancingdesign` gives the extended factory layout when rebalancing is enabled.
There is a new **ToggleCoordinatesTransferFactory** which is controlled by the **TogglePFactory**
and appropriately generates the coarse coordinates depending on the used transfer operator.
In case of semi-coarsening, the **SemiCoarsenPFactory** provides the coarse coordinates,
which are then piped through the **ToggleCoordinatesTransferFactory**.
In case of standard aggregation,
the **CoordinatesTransferFactory** calculates the coarse coordinates using the aggregation information provided by the **AggregationFactory**.
The coarse coordinate information is finally rebalanced by the **RebalanceTransferFactory** based on the rebalancing information provided by the **RepartitionFactory**.

.. note::

    Note, that the **LineDetectionFactory** algorithm expects all nodes of a vertical line along the extrusion axis of the underlying 2D mesh to be owned by the same processor.
    Do not allow for rebalancing before semi-coarsening is complete!
    Alternatively, you can implement a new interface class to replace the **ZoltanInterface** which makes sure that the nodes are rebalanced appropriately.
    This could be easily done by rebalancing the **VertLineIds** info
    that is provided by the **LineDetectionFactory** and reconstruct the node based **Partition** data.
    The **Chosen P** variable provided by the **TogglePFactory** would tell the new interface class whether we are in semi-coarsening mode or in standard aggregation mode.



.. _advanced_concepts/figure_rebalancedtogglerebalancingdesign:

.. warning::

    Insert missing figure here

The following listing shows exemplary the content of the xml file to set up a factory hierarchy similar to the one shown in Figure :ref:`advanced_concepts/figure_rebalancedtogglerebalancingdesign`.

In the **Factories** section,
first the **SemiCoarsenPFactory** and the line detection algorithm are defined representing the left transfer operator branch.
Next, the smoothed aggregation coarsening branch is defined in Part II
together with an instance of the **TransferCoordinatesTransferFactory** for calculating the corresponding coarse coordinates.
In Part III, the **TogglePFactory** is defined.
It contains a **TransferFactories** sublist where all different coarsening branches (i.e., the semi-coarsening and aggregation-based coarsening) are defined.
The corresponding prolongation operators are listed using the variable name **P** with a number between 1 and 9.
In addition to the prolongation operator factories,
one has also to declare the factories providing the coarse level null space and the *tentative* prolongation operator.
In case of semi-coarsening this is the **SemiCoarsenPFactory** and for the aggregation-based coarsening this is the **TentativePFactory**.
These factories are declared for generating the **Nullspace** variable with a number between 1 and 9 corresponding to the associated transfer operator branch.
Similar for the *tentative* prolongation operators denoted by the variable **Ptent**.

.. note::
    **SemiCoarsenPFactory** provides this information for compatibility reasons,
    even though there is no tentative prolongation operator for the geometrically generated **SemiCoarsenPFactory** operator.
    But the **TogglePFactory** is designed to be more general and allows for combining different kinds of smoothed prolongation operators.
    These need information about the non-smoothed transfer operators in the variable **Ptent**.


Similarly, a **ToggleCoordinatesTransferFactory** is declared with an internal list of all factories providing the coarse level coordinates.
This is the previously defined **CoordinatesTransferFactory** for the standard aggregation-based coarsening branch and the **SemiCoarsenPFactory** for the semi-coarsening branch.

Part IV contains the standard factories for the rebalancing.

Finally, it is important to declare all necessary main factories in the **Hierarchy** section of the xml file.

.. code-block:: xml

    <ParameterList name="MueLu">
    <ParameterList name="Factories">

        <!-- =======================  PART I  ======================= -->
        <ParameterList name="myLineDetectionFact">
        <Parameter name="factory" type="string" value="LineDetectionFactory"/>
        <Parameter name="linedetection: orientation" type="string" value="coordinates"/>
        </ParameterList>

        <ParameterList name="mySemiCoarsenPFact1">
        <Parameter name="factory" type="string" value="SemiCoarsenPFactory"/>
        <Parameter name="semicoarsen: coarsen rate" type="int" value="6"/>
        </ParameterList>

        <!-- =======================  PART II  ======================= -->
        <ParameterList name="UncoupledAggregationFact2">
        <Parameter name="factory" type="string" value="UncoupledAggregationFactory"/>
        <Parameter name="aggregation: ordering" type="string" value="graph"/>
        <Parameter name="aggregation: min agg size" type="int"    value="9"/>
        </ParameterList>

        <ParameterList name="MyCoarseMap2">
        <Parameter name="factory" type="string" value="CoarseMapFactory"/>
        <Parameter name="Aggregates" type="string" value="UncoupledAggregationFact2"/>
        </ParameterList>


        <ParameterList name="myTentativePFact2">
        <Parameter name="factory"     type="string" value="TentativePFactory"/>
        <Parameter name="Aggregates"  type="string" value="UncoupledAggregationFact2"/>
        <Parameter name="CoarseMap"   type="string" value="MyCoarseMap2"/>
        </ParameterList>

        <ParameterList name="mySaPFact2">
        <Parameter name="factory"     type="string" value="SaPFactory"/>
        <Parameter name="P"           type="string" value="myTentativePFact2"/>
        </ParameterList>

        <ParameterList name="myTransferCoordinatesFact">
        <Parameter name="factory"     type="string" value="CoordinatesTransferFactory"/>
        <Parameter name="CoarseMap"   type="string" value="MyCoarseMap2"/>
        <Parameter name="Aggregates"  type="string" value="UncoupledAggregationFact2"/>
        </ParameterList>

        <!-- =======================  PART III  ======================= -->

        <ParameterList name="myTogglePFact">
        <Parameter name="factory"              type="string" value="TogglePFactory"/>
        <Parameter name="semicoarsen: number of levels"       type="int" value="2"/>
        <ParameterList name="TransferFactories">
            <Parameter name="P1"                type="string" value="mySemiCoarsenPFact1"/>
            <Parameter name="P2"                type="string" value="mySaPFact2"/>
            <Parameter name="Ptent1"            type="string" value="mySemiCoarsenPFact1"/>
            <Parameter name="Ptent2"            type="string" value="myTentativePFact2"/>
            <Parameter name="Nullspace1"        type="string" value="mySemiCoarsenPFact1"/>
            <Parameter name="Nullspace2"        type="string" value="myTentativePFact2"/>
        </ParameterList>
        </ParameterList>

        <ParameterList name="myRestrictorFact">
        <Parameter name="factory"   type="string" value="TransPFactory"/>
        <Parameter name="P"         type="string" value="myTogglePFact"/>
        </ParameterList>

        <ParameterList name="myToggleTransferCoordinatesFact">
        <Parameter name="factory"   type="string" value="ToggleCoordinatesTransferFactory"/>
        <Parameter name="Chosen P"  type="string" value="myTogglePFact"/>
        <ParameterList name="TransferFactories">
            <Parameter name="Coordinates1" type="string" value="mySemiCoarsenPFact1"/>
            <Parameter name="Coordinates2" type="string" value="myTransferCoordinatesFact"/>
        </ParameterList>
        </ParameterList>

        <ParameterList name="myRAPFact">
        <Parameter name="factory" type="string" value="RAPFactory"/>
        <Parameter name="P"       type="string" value="myTogglePFact"/>
        <Parameter name="R"       type="string" value="myRestrictorFact"/>
        <ParameterList name="TransferFactories">
            <Parameter name="For Coordinates" type="string" value="myToggleTransferCoordinatesFact"/>
        </ParameterList>
        </ParameterList>

        <!-- =======================  PART IV (Repartitioning)  ======================= -->
        <ParameterList name="myZoltanInterface">
        <Parameter name="factory"      type="string" value="ZoltanInterface"/>
        <Parameter name="A"            type="string" value="myRAPFact"/>
        <Parameter name="Coordinates"  type="string" value="myToggleTransferCoordinatesFact"/>
        </ParameterList>

        <ParameterList name="myRepartitionFact">
        <Parameter name="factory"    type="string" value="RepartitionFactory"/>
        <Parameter name="A"          type="string" value="myRAPFact"/>
        <Parameter name="Partition"  type="string" value="myZoltanInterface"/>
        <Parameter name="repartition: min rows per proc"  type="int"    value="800"/>
        <Parameter name="repartition: max imbalance"      type="double" value="1.1"/>
        <Parameter name="repartition: start level"        type="int"    value="3"/>
        <Parameter name="repartition: remap parts"        type="bool"   value="false"/>
        </ParameterList>

        <ParameterList name="myRebalanceProlongatorFact">
        <Parameter name="factory"      type="string" value="RebalanceTransferFactory"/>
        <Parameter name="type"         type="string" value="Interpolation"/>
        <Parameter name="P"            type="string" value="myTogglePFact"/>
        <Parameter name="Coordinates"  type="string" value="myToggleTransferCoordinatesFact"/>
        <Parameter name="Nullspace"    type="string" value="myTogglePFact"/>
        </ParameterList>

        <ParameterList name="myRebalanceRestrictionFact">
        <Parameter name="factory"      type="string" value="RebalanceTransferFactory"/>
        <Parameter name="type"         type="string" value="Restriction"/>
        <Parameter name="R"            type="string" value="myRestrictorFact"/>
        </ParameterList>

        <ParameterList name="myRebalanceAFact">
        <Parameter name="factory"      type="string" value="RebalanceAcFactory"/>
        <Parameter name="A"            type="string" value="myRAPFact"/>
        </ParameterList>
    </ParameterList>

    <!-- Definition of the multigrid preconditioner -->
    <ParameterList name="Hierarchy">
        <Parameter name="max levels"       type="int"      value="6"/>
        <Parameter name="coarse: max size" type="int"      value="100"/>
        <Parameter name="verbosity"        type="string"   value="High"/>
        <ParameterList name="All">
        <Parameter name="P"              type="string"   value="myRebalanceProlongatorFact"/>
        <Parameter name="Nullspace"      type="string"   value="myRebalanceProlongatorFact"/>
        <Parameter name="CoarseNumZLayers" type="string"   value="myLineDetectionFact"/>
        <Parameter name="LineDetection_Layers" type="string"   value="myLineDetectionFact"/>
        <Parameter name="LineDetection_VertLineIds" type="string"   value="myLineDetectionFact"/>
        <Parameter name="A"              type="string"   value="myRebalanceAFact"/>
        <Parameter name="Coordinates"    type="string"   value="myRebalanceProlongatorFact"/>
        <Parameter name="Importer"       type="string"   value="myRepartitionFact"/>
        <!--<Parameter name="R"          type="string"   value="myRebalanceRestrictionFact"/>-->
    </ParameterList>
    </ParameterList>
    </ParameterList>

.. warning::

    Include the above XML file into testing.


.. _advanced_concepts/line-smoothing:

Line-smoothing
=============

General idea
------------

Semi-coarsening should be combined with line-smoothing as the complementary smoothing operation.
Whereas semi-coarsening coarsens, e.g., along the z-axis trying to produce a 2D representation of a 3D problem,
the line-smoothing operates orthogonal to the coarsening and smoothes in the x- and y-direction and interprets all vertical z-layers technically as nodes in a pseudo 2D problem.

Usage
-----

The following listing shows how to choose a Jacobi line smoother.
The reader might compare the xml code snippets with Section :ref:`level_smoothers` for a detailed description of the different smoothers and parameters.

.. code-block:: xml

    <ParameterList name="MueLu">
    <ParameterList name="Factories">
    <ParameterList name="mySmoother1">
        <Parameter name="factory"   type="string" value="TrilinosSmoother"/>
        <Parameter name="type"      type="string" value="LINESMOOTHING_BANDEDRELAXATION"/>
        <Parameter name="smoother: pre or post"        type="string" value="pre"/>
        <ParameterList name="ParameterList">
        <Parameter name="relaxation: type"           type="string" value="Jacobi"/>
        <Parameter name="relaxation: sweeps"         type="int"    value="2"/>
        <Parameter name="relaxation: damping factor" type="double" value="0.3"/>
        </ParameterList>
    </ParameterList>
    <ParameterList name="Hierarchy">
        <ParameterList name="All">
        <Parameter name="Smoother"    type="string"   value="mySmoother1"/>
        </ParameterList>
    </ParameterList>
    </ParameterList>

.. warning::

    Include the above XML file into testing.


The parameters are standard except of **type**.
The standard choice would be **RELAXATION** for relaxation based smoothers.
To use line-smoothing instead one has the following options:

.. admonition:: Description

   * [LINESMOOTHING\BANDEDRELAXATION] Use banded containers to store the local block associated with one vertical line.
   This is the recommended variant as it saves memory and is usually faster.
   * [LINESMOOTHING\BLOCKEDRELAXATION] Use a dense matrix container to store the local block associated with one vertical line.
   This is the safe fallback variant.
   Use the **LINESMOOTHING\_BANDEDRELAXATION** variant instead.


All the other parameters in the parameter sublist correspond to the usual parameters for relaxation based smoothers such as Jacobi, Gauss-Seidel or Symmetric Gauss-Seidel methods.
Refer to Section :ref:`level_smoothers` or the MueLu user guide [1]_ for an overview of all available parameters.

Footnotes
=========

.. [1] L. Berger-Vergiat, C. A. Glusa, G. Harper, J. J. Hu, M. Mayr, P. Ohm, A. Prokopenko, C. M. Siefert, R. S. Tuminaro, and T. A. Wiesner. MueLu User's Guide. Technical Report SAND2023-12265, Sandia National Laboratories, Albuquerque, NM (USA) 87185, 2023.
