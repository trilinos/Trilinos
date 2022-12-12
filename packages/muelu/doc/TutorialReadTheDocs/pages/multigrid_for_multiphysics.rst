==========================
Multigrid for Multiphysics
==========================

General Concept
===============

Assuming we have a multiphysics application which generates a :math:`n\times n` block operator
(e.g. a Fluid-Structure-Interaction problem,
or a Navier-Stokes problem with the splitting into velocity and pressure degrees of freedom),
the idea is to preserve the block structure on the coarse levels
using block-diagonal segregated transfer operators and block smoothers on all multigrid levels.
The blocks :math:`P_{i,i}` and :math:`R_{i,i}` in the block-diagonal transfer operators :math:`P` and :math:`R`
usually are built using the diagonal blocks :math:`A_{i,i}`
or - for applications where :math:`A_{i,i}` does not contain sufficient information to build aggregates -
any other kind of application-specific method.

For the block smoothers we apply well-established block smoothers
(e.g. block relaxation methods or Schur complement based methods)
as well as application-specific adaptions.

Exemplary setup for a :math:`2\times 2` problem
================================================

.. _multigrid_for_multiphysics/blocktransfersetup:

Setup of block transfer operators
---------------------------------

.. warning:: Insert missing figure.

Figure \ref{fig:transferoperatorsetup} shows a typical layout for the setup of a :math:`2\times 2` blocked operator block transfer operators.
We group factories for building transfer operators for the upper-left and lower-right blocks and build the corresponding blocked transfer operators.
In general, we can distinguish volume-coupled and interface-coupled problems.
In the simpler case of volume-coupled problems,
we usually can use the same aggregates for the second block than we have built for the first block.
That is, we reuse the same aggregates (built, e.g., by *myAggFact1*) for the second block.
For interface-coupled problems we usually need a separate aggregation strategy for the interface DOFs.
This can be a second standard aggregation factory object or an application-specific aggregation routine.

In the following we give more details on the corresponding xml files for setting up the transfer operator layout given in Figure \ref{fig:transferoperatorsetup}.

Factory list
------------
All the following definitions of factories are contained in the **Factories** sublist of the **MueLu** parameter list,
which basically is meant as a factory collection for all used factories in the multigrid setup.
To keep things simple, we only give the factories necessary to define the setup for the most upper-left block
in the :math:`n\times n` block matrix.


.. code-block:: xml

   <!-- sub block factories -->
    <!-- BLOCK 1 (for submatrix A_{00}) -->
    <ParameterList name="mySubBlockAFactory1">
      <Parameter name="factory" type="string" value="SubBlockAFactory"/>
      <Parameter name="block row"                 type="int"     value="0"/>
      <Parameter name="block col"                 type="int"     value="0"/>
      <Parameter name="Range map: Striding info"  type="string"  value="{ 2 }"/>
      <Parameter name="Domain map: Striding info" type="string"  value="{ 2 }"/>
    </ParameterList>

    <ParameterList name="myAggFact1">
      <Parameter name="factory" type="string" value="UncoupledAggregationFactory"/>
      <Parameter name="aggregation: min agg size" type="int" value="5"/>
      <Parameter name="aggregation: max selected neighbors" type="int" value="1"/>
    </ParameterList>

    <!-- tell the tentative prolongator that we have 2 DOFs per node on the coarse levels -->
    <ParameterList name="myCoarseMap1">
      <Parameter name="factory" type="string" value="CoarseMapFactory"/>
      <Parameter name="Striding info" type="string" value="{ 2 }"/>
      <Parameter name="Strided block id" type="int" value="-1"/>
    </ParameterList>

    <ParameterList name="myTentativePFact1">
      <Parameter name="factory" type="string" value="TentativePFactory"/>
      <Parameter name="A" type="string" value="mySubBlockAFactory1"/>
      <Parameter name="Aggregates" type="string" value="myAggFact1"/>
      <Parameter name="CoarseMap" type="string" value="myCoarseMap1"/>
    </ParameterList>

    <!-- We have to use Nullspace1 here. If "Nullspace1" is not set the
         Factory creates the default null space containing of constant
         vectors -->
    <ParameterList name="myNspFact1">
      <Parameter name="factory" type="string" value="NullspaceFactory"/>
      <Parameter name="Fine level nullspace" type="string" value="Nullspace1"/>
      <Parameter name="Nullspace1" type="string" value="myTentativePFact1"/>
    </ParameterList>

.. note::

   Please note, that the ordering of the factories is important in the sense
   that factories have to be defined before they are used as input for other following factories.
   As an example, the **myCoarseMap1** factory has to be defined before the **myTentativePFact1**
   as it is used as input for the **CoarseMap** variable.
   Switching the ordering of the factories in the xml file would result in an error that **myCoarseMap1** is unknown.
   Technically, one could avoid this restriction by a two pass reading process
   which first reads all factories and then resolves the dependencies.
   On the other side, this restriction helps to keep a straightforward linear design of the setup process.

The meaning of the factories is the following:

- **SubBlockAFactory**:
  Given a :math:`n\times n` block operator :math:`A`, the **SubBlockAFactory** extracts the :math:`(i,j)` block
  where :math:`i` is defined by the parameter **block row** and :math:`j` by the parameter **block col** where :math:`0\leq i,j < n`.
  Above example assumes a Thyra-style numbering of the global ids for a simple 2D Navier-Stokes example.
  That is, the matrix block :math:`A_{00}` has two degrees of freedom per node
  (one for the velocity in :math:`x`-direction and one in the :math:`y`-direction).
  The **Range map: Striding info** contains this information (i.e. 2 dofs per node),
  since this information might get lost (or never was stored) when using Thyra block operators.
- **UncoupledAggregationFactory**:
  The aggregation factory is used to build the aggregates.
  In this case, the aggregates shall be built using the graph of :math:`A_{00}` that is returned by the **SubBlockAFactory**.
  In this example, we only give the user parameters for the aggregation.
  Later it is shown how to declare the **FactoryManager**
  which makes sure that this concrete instance of an aggregation factory builds the aggregates for :math:`A_{00}` only.
- **CoarseMapFactory**:
  The **CoarseMapFactory** is used in the **TentativePFactory** and
  basically is responsible to provide a proper domain map for the transfer operator block :math:`P_i`.
  For :math:`P_0`, this is usually a standard map.
  The only information that is important is **Striding info** which means that the coarse domain map has 2 dofs per node again.
  Note: we have 2 dofs per node (for the velocities in :math:`x` and :math:`y`-direction).
  We have 2 null space vectors.
  Therefore, the coarse problem has also 2 dofs per node which means the domain map of :math:`P_i` has to be built for 2 dofs per node.
- **TentativePFactory**:
  Here, the **TentativePFactory** builds the :math:`P_i` block for the blocked transfer operator :math:`P:math:`.
  We explicitly give the names of the factories used to generate :math:`P_i`,
  which include the previously defined factories for the **CoarseMap**, **Aggregates** and **A**.
  This information is not really needed in this place as we later define a **FactoryManager** for the :math:`i`-th block,
  but often makes it easier to understand the dependencies.
  That is, the short version would just be

  .. code-block:: xml

      <ParameterList name="myTentativePFact1">
        <Parameter name="factory" type="string" value="TentativePFactory"/>
      </ParameterList>

- **NullspaceFactory**:
  For defining multigrid methods for multiphysics problems, the **NullspaceFactory** is very important.
  In general, the **TentativePFactory** uses the provided fine level near null space vectors
  to generate the tentative prolongation operator :math:`P_i` together with a coarse representation of the near null space vectors.
  That is, the **TentativePFactory** produces the **Nullspace** information that it needs itself on the next coarser level.
  That is a hidden dependency which usually automatically works without any changes necessary by the user.
  The user is only responsible to provide proper fine level near null space vectors as **Nullspace** variable on the finest level.
  The **NullspaceFactory** is just a helper factory which processes the null space information on the finest level
  and pipes it in into the global setup process.
  For multiphysics problems, the user has to provide :math:`n` partial near null space vectors (one for each mathematical or physical field)
  using the variable names **Nullspace1** to **Nullspace9** on the finest level.
  The **Fine level nullspace** parameter in the **NullspaceFactory** then can be set to the corresponding variable name (e.g. \texttt{Nullspace1}).
  That is, the **NullspaceFactory** checks the fine level variable container for a variable named **Nullspace1**
  and uses the content as fine level null space for input in the **TentativePFactory**.
  It is important, that besides of the **Fine level nullspace** parameter another parameter with the name of the near null space vector
  (in above case **Nullspace1**)
  is declared with the corresponding **TentativePFactory** name as value.
  This closes the circle for the null space generation for block :math:`P_i` on all coarser levels.
  It is important that the **NullspaceFactory** is defined after the corresponding **TentativePFactory** class,
  such that the dependency circle can be closed correctly.

.. note:: Instead of **TentativePFactory** any other factory which generates a coarse null space can be used as well (e.g. **SemiCoarsenPFactory**).

.. note:: Of course, above factory list can be extended for smoothed aggregation. We also skipped the factories for the restriction operators.

Factory manager
---------------

Once the necessary factories for building :math:`P_i` are defined in the **FactoryList** section of the xml file,
we can group them together.
Right after the factories, we can add a **FactoryManager** block in the **FactoryList** section.

.. code-block:: xml

    <!-- Multigrid setup for velocity block (A_{00}) -->
    <ParameterList name="myFirstGroup">
      <Parameter name="group" type="string" value="FactoryManager"/>
      <Parameter name="A" type="string" value="mySubBlockAFactory1"/>
      <Parameter name="P" type="string" value="myTentativePFact1"/>
      <Parameter name="Aggregates" type="string" value="myAggFact1"/>
      <Parameter name="Nullspace" type="string" value="myNspFact1"/>
      <Parameter name="CoarseMap" type="string" value="myCoarseMap1"/>
    </ParameterList>

The name for the group can be chosen freely (e.g. **myFirstGroup**).
Besides the declaration of the **FactoryManager** group using the

.. code-block:: xml

   <Parameter name="group" type="string" value="FactoryManager"/>

parameter, it contains a list of all factories
which are used in context of building the coarsening information for the corresponding block.

The group block defining a **FactoryManager** has a similar role than in the **Hierarchy** section in the xml file later.
It allows to group together factories into subgroups
that can be referred to by the common name (e.g. **myFirstGroup**) later.
These groups help to organize the different factories.
Note, that we basically need one group for each physical/mathematical field in our :math:`n\times n` block operator,
that is we need :math:`n` groups.

.. note::

   The group block for the second row in the block operator could look like the following

   .. code-block:: xml

       <!-- Multigrid setup for pressure block (A_{11}) -->
       <ParameterList name="mySecondGroup">
         <Parameter name="group" type="string" value="FactoryManager"/>
         <Parameter name="A" type="string" value="mySubBlockAFactory2"/>
         <Parameter name="P" type="string" value="myTentativePFact2"/>
         <!-- reuse aggs from velocity block! -->
         <Parameter name="Aggregates" type="string" value="myAggFact1"/>
         <Parameter name="Nullspace" type="string" value="myNspFact2"/>
         <Parameter name="CoarseMap" type="string" value="myCoarseMap2"/>
       </ParameterList>

   This assumes that all the factories have been defined before similar to the first group blocks.
   Note, that in some cases for certain applications,
   it is possible to reuse information from the first block in the second block.
   In this case, we use the abstract aggregation information
   that has been built using the velocity information for the associated pressure degrees of freedom.
   This is possible, since in our example each node has 2 velocity and 1 pressure degree of freedom.

Block prolongation operator
---------------------------

The diagonal block prolongation operator is built using

.. code-block:: xml

   <!-- define block prolongation operator using above blocks -->
   <ParameterList name="myBlockedPFact">
     <Parameter name="factory" type="string" value="BlockedPFactory"/>
     <!-- factory manager for block 1 -->
     <ParameterList name="block1">
       <Parameter name="group" type="string" value="myFirstGroup"/>
     </ParameterList>
     <!-- factory manager for block 2 -->
     <ParameterList name="block2">
       <Parameter name="group" type="string" value="mySecondGroup"/>
     </ParameterList>
   </ParameterList>

in the **FactoryList** section of the xml file after all groups have been defined.
It contains basically sublists with the names **block1** to **blockn**.
Each of these sublists contains a parameter **group** with the group name defined before.

.. note::

   Instead of using groups,
   you could also put all the factory definitions within the corresponding **blockn** parameter list.
   But this would mean that you have to set all inter-factory dependencies in the corresponding block by hand.
   You cannot use the general defaults that are defined in the groups.
   It would also make it impossible to reuse information from factories belonging to a different block
   (e.g., you could not reuse the aggregation information built by the **myAggFact1** for the aggregates in block 2.

Block restriction operator
--------------------------

The following definitions should be the standard for nearly all multiphysics problems.
We use the **GenericRFactory** for building the restriction operator out of the blocked prolongation factory.

.. code-block:: xml

   <!-- define block restriction operator using above blocks -->
   <!-- The block restriction operator is usually always of type
        GenericRFactory since we want to be able to combine, e.g.,
        SmoothedAggregation for block A_{00} with e.g. tentative
        prolongation for block A_{11} (or any other kind of transfer
        strategy for the subblocks -->
   <ParameterList name="myBlockedRFact">
     <Parameter name="factory" type="string" value="GenericRFactory"/>
     <Parameter name="P" type="string" value="myBlockedPFact"/>
   </ParameterList>

It uses the blocks :math:`P_i` to generate the corresponding blocks :math:`R_i` for the diagonal of the restriction operator.
If PG-AMG is used for some or all blocks :math:`P_i` this is automatically considered when generating :math:`R_i`.

.. note::

   Please note, that you cannot use **TransPFactory** as it has no support for block operators.
   However, this is not a problem, since the **TransPFactory** might be used locally for single blocks wherever possible.

Coarse level operator
---------------------

Once, the block diagonal transfer operators :math:`P` and :math:`R` are set up,
the **BlockedRAPFactory** builds the coarse :math:`n\times n` operator:

.. code-block:: xml

    <ParameterList name="myBlockedRAPFact">
      <Parameter name="factory" type="string" value="BlockedRAPFactory"/>
      <Parameter name="P" type="string" value="myBlockedPFact"/>
      <Parameter name="R" type="string" value="myBlockedRFact"/>
    </ParameterList>

Setup of block smoothers
========================

Once the transfer operators are set up properly, it is time to define the block smoothing methods

Xpetra block smoothers
----------------------

The Xpetra package contains a set of general block smoothers,
including a block Gauss-Seidel method for :math:`n\times n` block operators
and a Schur complement based SIMPLE variant for :math:`2\times 2` block operators.
Here we just explain the setup for the Schur complement smoother as an example.

The Schur complement based smoother internally needs two solvers/smoothers.
The first smoother/solver is need for calculating a prediction for the velocities,
the second solver/smoother then has to (approximately) solve the Schur complement equation.

We are still in the **FactoryList** section of the xml file and define the corresponding solver/smoother blocks:

.. code-block:: xml

   <!-- block smoother for block A_{00} -->
   <ParameterList name="mySmooFact1">
     <Parameter name="factory" type="string" value="TrilinosSmoother"/>
     <Parameter name="type" type="string" value="RELAXATION"/>
     <ParameterList name="ParameterList">
       <Parameter name="relaxation: type" type="string" value="Symmetric Gauss-Seidel"/>
       <Parameter name="relaxation: sweeps" type="int"    value="3"/>
       <Parameter name="relaxation: damping factor" type="double" value="0.8"/>
     </ParameterList>
     <Parameter name="A" type="string" value="mySubBlockAFactory1"/>
   </ParameterList>

   <!-- Build Schur Complement factory (for being used in SIMPLE type block smoother)-->
   <ParameterList name="mySchurCompFact">
     <Parameter name="factory" type="string" value="SchurComplementFactory"/>
     <Parameter name="omega" type="double" value="0.8"/>
     <Parameter name="lumping" type="bool" value="false"/>
   </ParameterList>

   <!-- block smoother for block A_{11} respective the Schur complement operator -->
   <ParameterList name="mySchurSmooFact">
     <Parameter name="factory" type="string" value="TrilinosSmoother"/>
     <Parameter name="type" type="string" value="RELAXATION"/>
     <ParameterList name="ParameterList">
       <Parameter name="relaxation: type" type="string" value="Symmetric Gauss-Seidel"/>
       <Parameter name="relaxation: sweeps" type="int"    value="1"/>
       <Parameter name="relaxation: damping factor" type="double" value="0.7"/>
     </ParameterList>
     <!-- You don't have to specify the input matrix A for the smoother -->
     <!-- It is clear from the subblocks in the BlockedSmoother below -->
     <!--<Parameter name="A" type="string" value="mySchurCompFact"/>-->
   </ParameterList>

It is not necessary but helpful to declare variable **A** in **mySmooFact1**
to be the diagonal block :math:`A_{00}` of the blocked operator.
This way it is obvious that this smoother is supposed to generate the prediction within the Schur complement approach.
The second smoother (with the name **mySchurSmooFact** is supposed to solve the Schur complement equation,
that is, the input matrix :math:`A` for this smoother should be the Schur complement operator :math:`A_{11}-A_{10}A_{00}^{-1}A_{01}`
or at least a good approximation of the Schur complement operator.
This operator is provided by the **SchurComplementFactory**.
Be aware, that the **SchurComplementFactory** uses the full :math:`2\times 2` operator :math:`A` as input
to generate the approximation of the Schur complement operator.
It is not defined as input variable **A** since the full :math:`2\times 2` operator is the standard answer for variable **A**.
It would make sense, though, to declare **mySchurCompFact** as variable **A** for **mySchurSmooFact**.

The Schur complement smoother then is defined by the block:

.. code-block:: xml

   <!-- Use SIMPLE: -->
   <!-- User has to define two blocks with each containing a smoother for
        the corresponding sub-block matrix (see above) -->
   <ParameterList name="myBlockSmoother">
     <Parameter name="factory" type="string" value="SimpleSmoother"/>
     <Parameter name="Sweeps" type="int" value="1"/>
     <Parameter name="Damping factor" type="double" value="1.0"/>
     <Parameter name="UseSIMPLEC" type="bool" value="false"/>
     <!-- factory manager for block 1 -->
     <ParameterList name="block1">
       <!-- <Parameter name="group" type="string" value="myFirstGroup"/> -->
       <!-- It's enough to just provide the sub block matrix A_{00} needed
            as input for the smoother "mySmooFact1" as well as the smoother
            itself. Alternatively, one could add below strings to the
            factory group "myFirstGroup" and use it here instead of below lines -->
       <Parameter name="A" type="string" value="mySubBlockAFactory1"/>
       <Parameter name="Smoother" type="string" value="mySmooFact1"/>
     </ParameterList>
     <!-- factory manager for block 2 -->
     <ParameterList name="block2">
       <!-- The second equation in the SIMPLE smoother is the Schur complement
            equation that has to be solved. Therefore, provide the Schur complement
            operator together with the smoother. The smoother object takes the
            Schur complement operator as operator "A" for the internal smoothing process. -->
       <Parameter name="A" type="string" value="mySchurCompFact"/>
       <Parameter name="Smoother" type="string" value="mySchurSmooFact"/>
     </ParameterList>
   </ParameterList>

In this case, we use the **SimpleSmoother** as example.
Besides the typical smoother parameters (number of sweeps, damping, \ldots),
the interesting part are the sublists **block1** and **block2**,
which contain the information about the internal smoothers/solvers.
In above example, we just declare the factories for **A** and **Smoother**.
The variable **A** always gives the internal linear operator that is used within the solver/smoother.
By defining **A** in this place, we do not really have to define it extra in the smoother blocks above.

.. note::

   Please note, that instead of the explicit variable definitions in the **blockn** sublists,
   one could also have given just the group names.
   However, this only works if the **Smoother** variable is also contained in the corresponding groups.
   In above examples from the previous section we skipped the **Smoother** variable.
   This makes sense especially if the aggregation information is built using a different :math:`A` operator
   as the smoother is using.
   In above example we do not build aggregates using the Schur complement operator,
   but want to reuse aggregates from the first block.

As a side note it shall be mentioned,
that you can also directly make all definitions in the **blockn** parameter lists:

.. code-block:: xml

   <ParameterList name="block1">
     <!-- <Parameter name="group" type="string" value="myFirstGroup"/> -->
     <ParameterList name="Smoother">
       <Parameter name="factory" type="string" value="TrilinosSmoother"/>
       <Parameter name="type" type="string" value="RELAXATION"/>
       <ParameterList name="ParameterList">
         <Parameter name="relaxation: type" type="string" value="Symmetric Gauss-Seidel"/>
         <Parameter name="relaxation: sweeps" type="int"    value="3"/>
         <Parameter name="relaxation: damping factor" type="double" value="0.8"/>
       </ParameterList>
     </ParameterList>
     <!-- Note: this is not recommended as it creates a second instance of the sub block AFactory -->
     <ParameterList name="A">
       <Parameter name="factory" type="string" value="SubBlockAFactory"/>
       <Parameter name="block row"                 type="int"     value="0"/>
       <Parameter name="block col"                 type="int"     value="0"/>
       <Parameter name="Range map: Striding info"  type="string"  value="{ 2 }"/>
       <Parameter name="Domain map: Striding info" type="string"  value="{ 2 }"/>
     </ParameterList>
   </ParameterList>

However, this is not really recommended since it prevents the reuse of factories in several places.
E.g., instead of the new **SubBlockAFactory** one should just reuse the **SubBlockAFactory**,
which has been defined and used before for the block transfers.
This drastically simplifies and shortens the factory definitions and reduces the number of potential errors.

.. warning::

   Be aware, that each new block in the xml file means that a new instance of the corresponding factory is instantiated and built.
   In the worst case some expensive information is calculated twice,
   which might heavily impact the overall performance.

Teko block smoothers
--------------------

.. note::
   In Trilinos, the Teko package provides block preconditioners,
   that can be used as alternative to the existing Xpetra block smoothers.
   The Xpetra linear algebra layer also provides support for Thyra block operators,
   which allows us to use the Teko block smoothers within a MueLu multigrid hierarchy.
   In case of the Teko SIMPLE implementation, one again needs to internal solvers/smoothers
   (one for the prediction of the primary variables and one for the solution of the Schur complement equation).
   Teko uses the Stratimikos interface for defining the corresponding smoothers/solvers.
   So, instead of the **SimpleSmoother** object from the previous subsection,
   one can also use the SIMPLE implementation from Teko.

We define a **TekoSmoother** as block smoother using

.. code-block:: xml

   <!-- Use SIMPLE: -->
   <!-- User has to define two blocks with each containing a smoother for
        the corresponding sub-block matrix-->
   <ParameterList name="myTekoSmoother">
     <Parameter name="factory" type="string" value="TekoSmoother"/>
     <Parameter name="Inverse Type" type="string" value="SIMPLE"/> <!-- contains name of sublist within Teko parameters -->
     <ParameterList name="Inverse Factory Library">
       <ParameterList name="SIMPLE">
         <Parameter name="Type" type="string" value="NS SIMPLE"/>
         <Parameter name="Inverse Velocity Type" type="string" value="Amesos"/>
         <Parameter name="Inverse Pressure Type" type="string" value="Amesos"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>

The **TekoSmoother** accepts the full :math:`2\times 2` block operator as input
(not declared above, since it is the default)
and contains a sublist with the name **Inverse Factory Library**.
Within this sublist, all local smoothers/solvers as well as the Teko block smoother
(or several Teko block smoothers) are defined.
In the above example, there is only one Teko block smoother (of type **NS SIMPLE**) declared,
which internally uses direct solvers from the Amesos package for the velocity and pressure (Schur complement) problem.
The **Inverse Type** parameter of the **TekoSmoother** defines the Teko block smoother from the **Inverse Factory Library**.
For the available parameters and block smoothers in Teko,
the reader is referred to the Teko documentation.

Multigrid setup
===============

Last but not least, once both the transfers and the block smoothers are defined,
the multigrid method itself has to be set up.
Note, that all previous definitions and declarations have been made in the **Factories** section of the xml file.
The multigrid setup is now done in the **Hierarchy** section of the xml file and looks like:

.. code-block:: xml

   <ParameterList name="Hierarchy">

     <Parameter name="max levels"          type="int"      value="3"/>
     <Parameter name="coarse: max size"    type="int"      value="10"/>
     <Parameter name="verbosity"           type="string"   value="High"/>

     <ParameterList name="AllLevel">
       <Parameter name="startLevel"        type="int"      value="0"/>
       <Parameter name="Smoother"          type="string"   value="myTekoSmoother"/>
       <Parameter name="CoarseSolver"      type="string"   value="myTekoSmoother"/>
       <Parameter name="P"                 type="string"   value="myBlockedPFact"/>
       <Parameter name="R"                 type="string"   value="myBlockedRFact"/>
       <Parameter name="A"                 type="string"   value="myBlockedRAPFact"/>
     </ParameterList>
   </ParameterList>

The interesting part is the **AllLevel** sublist (you can freely choose the name of this list),
which - in some sense - corresponds to the groups introduced before to setup the block transfers and block smoothers.
In fact, this sublist defines the master **FactoryManager** for the overall multigrid method.
Note, that all variables (**A**, **P**, **R**, **...**) are generated by the block versions
instead of the single block factories.

Exemplary setup for a :math:`2\times 2` problem with rebalancing
================================================================

Transfer operator setup
-----------------------

.. warning:: Include missing figure.

Figure \ref{fig:transferoperatorsetuprebalancing} shows the basic setup for block transfer operators with rebalancing enabled.
Please compare it with the complete XML input deck in Section \ref{sec:xmlinputdeckrebalancing}.

As one can see from the upper part of Figure \ref{fig:transferoperatorsetuprebalancing},
first we build blocked transfer operators and a blocked coarse level operator
using sub-factory manager objects **myFirstGroup** and **mySecondGroup**
in the factories **myBlockedPFact**, **myBlockedRFact** and **myBlockedRAPFact**.
Then, we rebalance the coarse level blocked operator :math:`A` from **myBlockedRAPFact**.

The **myRepartitionHeuristicFact** object will decide whether rebalancing is necessary.
If yes, then it will return the number of required partitions for the coarse level operator.
This input is processed by the repartition interface and repartition factory objects
that finally create **Xpetra::Importer** to do the rebalancing.
The **myRebBlocked{P,R,Ac}Fact** objects use those **Importer** objects to perform the rebalancing.

Please note, that we build additional helper factory manager objects **myRebFirstGroup** and **myRebSecondGroup**
which contain all factories relevant for rebalancing the two blocks.

.. note::

   No changes are necessary when setting up the block smoothers,
   as they use the matrices on the current level as input
   (which may or may not be rebalanced in the previous transfer operator setup process).

Complete XML input deck
-----------------------

.. literalinclude:: ../../../test/tutorial/blocked_rebalancing.xml
  :language: xml
  :caption:

Note, that we are using a coordinate-based rebalancing method from the Zoltan package.
The **myInputCoordsFact** provides the **Coordinates** variable to the **ZoltanInterface**.
That is, **myInputCoordsFact** uses user-provided data on the finest level
and switches to Coordinates provided by the **myTransferCoordinatesFact** on the coarser level.

.. note::

   To deal with potential logical circular dependencies between factories,
   you can use the **dependency for** keyword as demonstrated for the **myTransferCoordinatesFact**
   and **myInputCoordsFact** in this example.
   Note: you can always use the **dependency for** keyword to extend/change factory dependencies in the xml file.

