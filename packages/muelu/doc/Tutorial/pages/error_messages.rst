==============
Error Messages
==============

This section summarizes typical error messages and hints how to address the underlying problem.

Syntax errors
=============

Parser errors
--------------------------

::

    XML parse error at line 27: file ended before closing element 'ParameterList' from line 1

It seems like you forgot to close the **<ParameterList>** section that is opened in line 1 of the xml file.
Close it by adding the line **</ParameterList>** at the end of the (sub-)list.

::

    XML parse error at line 15: start element not well-formed: invalid character

Check line 15 for an invalid xml format.
The reason can be, e.g., a missing closing character **/>** for a parameter.

Parameter list errors
---------------------
::

    All child nodes of a ParameterList must have a name attribute!

You probably forgot to add a name attribute in one or more elements of your xml file, that is you used, e.g.,

* <Parameter  type="string" value="RELAXATION"/>

instead of

* <Parameter name="smoother: type" type="string" value="RELAXATION"/>

::

    Error, the parameter {name="smoother: type",type="int",value="0"} in the parameter (sub)list "ANONYMOUS" exists in the list of valid parameters but has the wrong type.

Choose the coorect type.
The correct type depends on the parameter at hand.
In this example, the correct type is "string".

::

    Use the correct (proposed) value type for the given parameter name, i.e.,

* <Parameter name="smoother: type" type="string" value="RELAXATION"/>

instead of

* <Parameter name="smoother: type" type="int" value="RELAXATION"/>


MueLu errors
============

General errors
--------------

::

    Throw test that evaluated to true: s_.is_null()

    Smoother for Tpetra was not constructed
        during request for data "    PreSmoother" on level 0 by factory NoFactory

Failed to create a level smoother. Check the smoother blocks in your xml file.
The error occurs, e.g., if there is a typing error in the **smoother: type** parameter.
For example

* <Parameter name="smoother: type" type="string" value="REXATION"/>

would trigger above error since the smoother type should be **RELAXATION**.

::

    IFPACK ERROR -2, ifpack/src/Ifpack_PointRelaxation.cpp, line 117

Errors like this indicate that it is a problem within the **smoother: params** section. Most likely a (relaxation) smoother is requested which is not existing (e.g., **Jadobi** instead of **Jacobi**).

.. warning::

    Switch this example to Ifpack2.

::

    The parameter name "smother: type" is not valid. Did you mean "smoother: type"?

There is a typo in your parameter list.
Locate the parameter and fix it (possibly using the suggestions, that come with the error message).

::

    Throw test that evaluated to true: maxNodesPerAggregate < minNodesPerAggregate

Choose the **aggregation: min agg size** parameter to be smaller than the **aggregation: max agg size** parameter for the aggregation routine.

Advanced XML file format
------------------------

::

    Throw test that evaluated to true: bIsZeroNSColumn == true

    MueLu::TentativePFactory::MakeTentative: fine level NS part has a zero column

This error indicates that there is a problem with the provided near null space vectors. There are different reasons which can trigger this problem:

* The near null space vectors are not valid (containing zeros, wrong ordering of internal degrees of freedom).
Please check your near null space vectors.
Maybe there is an empty vector or the ordering of degrees of freedom for the linear operator
does not match with the ordering of the near null space vectors.
* The near null space vectors are correct but used in a wrong way (e.g., a wrong number of degrees of freedom).
Check the screen output for wrong block dimensions (CoalesceDropFactory).
* There is a problem with the aggregates.
Validate the screen output and look for unusual (e.g. very small or empty) aggregates.

::

    Throw test that evaluated to true: factoryManager_ == null

    MueLu::Level(0)::GetFactory(Aggregates, 0): No FactoryManager

This is a typical error when the dependency tree is screwed up.
If aggregates and/or transfer operators are involved,
one usually has forgotten some entries in the **Hierarchy** sublist of the extended XML file format
for the internal factory managers.
These errors can be quite tricky to fix.
In general, it is a good idea to start with a working XML file and extend it step by step if possible.
The following general strategies may help to track down the problem:

* Run the problem with **verbosity=extreme** to get as much screen output as possible.
Check for unusual screen output (such as **Nullspace factory**).
* Try to generate a graphical dependency tree as described in :ref:`useful_commands_and_debugging/dependencytrees`.


For example, above error is caused by the following XML file

.. code-block:: xml

    <ParameterList name="MueLu">
    <ParameterList name="Factories">
        <ParameterList name="myTentativePFact">
        <Parameter name="factory" type="string" value="TentativePFactory"/>
        </ParameterList>
    </ParameterList>

    <ParameterList name="Hierarchy">
        <ParameterList name="Levels">
        <Parameter name="P" type="string" value="myTentativePFact"/>
        <!--<Parameter name="Nullspace" type="string" value="myTentativePFact"/>-->
        </ParameterList>
    </ParameterList>
    </ParameterList>

When looking at the error output, it seems to be a problem with aggregates.
However, no special aggregation factory has been declared in the XML file.
The only factory which has been introduced was a tentative prolongation factory for generating unsmoothed transfer operators.
Therefore, one should start digging into the details of the **TentativePFactory** to find out,
that the unsmoothed transfer operator factory is responsible both for creating the unsmoothed prolongator and the coarse level null space information.
When looking at the screen output one should find that the last called/generated factory is a **NullspaceFactory**,
which can also be a hint that the problem is the null space.

When looking at the XML file,
one can see that the **myTentativePFact** factory has been registered to be responsible for generating  the prolongator :math:`P`,
but the generating factory for the variable **Nullspace** is not declared.
MueLu tries to generate the default null space,
but since it does not know about **myTentativePFact** to be a **TentativePFactory**,
which would already produce the needed information,
the call ordering of the dependent factories (e.g., aggregation) gets mixed up.

Note that the **TentativePFactory** is special.
If you declare an explicit instance of the **TentativePFactory**,
you always have to register it for generating the **Nullspace** variable, too.
Only in very special cases, this would not be necessary.

.. note::

    This is a general rule: if a factory generates more than one output variables,
    always make sure that all these output variables are properly defined in the **FactoryManager** list (or **Hierarchy** sublist in the xml files, respectively).

To solve above problem there are two possibilities:

* Following above comment, just register **myTentativePFact** for generating **Nullspace**.
That is, just comment in the corresponding line in above xml file.
* Alternatively, you can register **myTentativePFact** for generating **Ptent** (and **P**).
This way you mark the **myTentativePFact** object to be used for generating the unsmoothed transfer operators
(and state that they shall be used for the final prolongation operators).
MueLu is smart enough to understand that the factory responsible for generating **Ptent** is also supposed to generate the null space vectors.
