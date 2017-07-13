========================================================================
	README: Zoltan2 Test Driver XML input files
========================================================================

This document details general formatting and requirements for Zoltan2 .xml input files used by the Zoltan2 test driver (/packages/zoltan2/test/driver/test_driver.exe).  A template for creating you own input files is included in this directory (input_template.xml).

========================================================================


Every Zoltan2 test driver input file must contain 2 required sections sections and may contain a 3rd optional section:

	1. Input source definition
	2. Test/Problem parameters and components definitions
	3. Comparison definitions (optional)

Each section is detailed below.  Please note that the sections described must be contained within a “Main” xml block whose name is arbitrary (ex. 0). Please not that the general use of the word “block” in the following refers to an xml parameter list block unless otherwise stated.


========================================================================
Section 1: Input source definition (REQUIRED)
========================================================================

In the first section of your input file you should define the input data source and type.  The test driver is designed such that all problems defined in section 2 share a common data source, therefore only 1 XML input source definition block per input file source is supported.  Currently there are 2 different flavors of input source definitions: one defining input from some supported file type, and the second defining Galari generated data.  Input parameter blocks should be named “InputParameters” and be the first defined block with in the main XML block.

 A block defining input from an existing file source must contain the following parameters (ex. 1):

	* input path: path to the directory containing the input file (relative or absolute)
	* input file: the name of the source file — EXCLUDING extension.
	* file type: Matrix Market, Pamgen, Chaco, or Geometric Generator


A Galari generated block should contain the following parameters (ex. 2):
	
	* x: number of grid points in the x-direction
	* y (optional): number of grid points in the y-direction
	* z (optional): number of grid points in the z-direction
	* equation type: Galari specific equation name, i.e. Laplace3D

An input source block defining a Galari generated data source may define a problem in 1, 2, or 3D.  The dimension of the problem is inferred by which coordinate parameters (x,y,z) have been defined, e.g., if only the x parameter has been defined then the problem is assumed to be 1D.

Both types of input blocks also support optional boolean parameters: distribute input and debug. Both parameters default to “True” when undefined.  Distribute input is only applicable for Chaco and MatrixMarket input formats.  Debug determines the verbosity of the UserInputForTests output stream.
 
======================================================================== Section 2: Zoltan2 problem definition (REQUIRED)
========================================================================

This section contains all of the blocks and associated sub-blocks that define a Zoltan2 problem.  This section must contain at least one problem definition and may contain as many as the user likes for the sake of testing and comparison. Each problem definition block should be uniquely named and must contain a ‘kind’ parameter specifying the kind of Zoltan2 problem: partitioning, coloring, ordering etc.  Each problem definition block must contain the following 2 sub-blocks (ex. 3):
	
	* InputAdapterParameters
	* Zoltan2Parameters

InputAdapterParameters:  This block defines the type of input adapter to be passed to the Zoltan2 problem as well as which data structure it should be constructed with.  Therefore this block is required to contain the following 2 parameters:
	
	* data type:
		- coordinates
		- (x,t,e)petra_vector
		- (x,t,e)petra_multivector
		- (x,t,e)petra_crs_matrix
		- (x,t,e)petra_crs_graph
	* input adapter
		- BasicIdentifier
		- BasicVector
		- XpetraMultivector
		- XpetraCrsMatrix
		- XpetraCrsGraph
		- PamgenMesh

Please note that if you choose to use a multi-vector data type with an adapter then you must additionally define an integer typed “vector_dimension” parameter specifying the dimension of the (x,t,e)petra multi-vector data type.

Zoltan2Parameters:  This block defines all of the parameters applicable to a given problem.  It supports all Zoltan2 parameters as well as those defined by supported TPLs.  Please consult the Zoltan2 documentation for a complete list of parameters.

======================================================================== Section 3: Metric definitions (OPTIONAL)
========================================================================
In addition to the aforementioned required blocks, a problem definition block may contain an optional 3rd “Metrics” block (ex. 4 and 5).
Sublists titled "metriccheck1", "metriccheck2", etc define each check to be conducted.

For example:
     <ParameterList name="Metrics">
       <ParameterList name="metriccheck1">
         <Parameter name="check" type="string" value="imbalance"/>
         <Parameter name="weight" type="int" value="0"/>
         <Parameter name="lower" type="double" value="0.99"/>
         <Parameter name="upper" type="double" value="1.4"/>
       </ParameterList>
    </ParameterList>
    
The 'check' value is type string and can accept the following key names which correspond to API calls in EvaluatePartition
	* ’imbalance’ is for ImbalanceMetrics
	* ’total edge cuts’ is for GraphMetrics
	* ’max edge cuts’ is for GraphMetrics
	
Additional values which can be specified:
	* ’weight’ is optional int type 0 or greater
	* 'normed' is optional bool type for ImbalanceMetrics only - not compatible with the weight option
	* 'lower' is double type and at least one of lower or upper should be specified
	* 'upper' is double type and at least one of lower or upper should be specified
 ======================================================================== Section 4: Comparison definitions (OPTIONAL)
========================================================================
 This section is optional, and may be defined to compare of solutions, metrics, or timers for different algorithms/adapters defined in section 2.  Like section 2 this section may include multiple “Comparison” blocks each specifying two problems/tests to compare.  For solution comparisons we compare problem “A” and problem “B” .  For metric or timer comparisons we compare a “Problem” vs. a “Reference”.  

A  solution “Comparison” block must contain the following 2 parameters (ex. 6):
	
	* A: the name of problem A
	* B: the name of problem B

The value of parameter A and B must be identical to the names of the blocks defining problems A and B in section 2.

A  metric or timer “Comparison” block must contain the following 2 parameters, followed by 1 or more sub-blocks defining metrics or timer lower and/or upper bounds (ex. 7):
	
	* Problem: the name of the problem to compare against the reference
	* Reference: the problem to be used as a reference

Supported timers are as follows: 

	* adapter construction time
	* problem construction time
	* solve time

The timers refer to the the time spent to construct an adapter, the total time spent to construct a Zoltan2 problem, and the time spent to solve the problem.  Pass/fail of a metric comparison is determined as follows:

reference_metric = (Reference’s Max Metric of a Part) / AVG
problem_metric = (Problem’s Max Metric of a Part) / AVG

PASS criteria is:
reference_metric * LOWER <= problem_metric <= reference_metric * UPPER

Pass/fail of a timer comparison is determined as follows:

PASS criteria is:
reference_time * LOWER <= problem_time <= reference_time * UPPER

Where lower and upper refer to the defined lower and/or upper bounds for tolerance.  As before, the values defining the parameters “Problem” and “Reference” must be the same as the values used to define the problems/tests in section 2.  

======================================================================== EXAMPLES:
========================================================================


////////////////////////////////////////////////////////
Example 0: General xml input structure
////////////////////////////////////////////////////////

<ParameterList name=“ARBITRAY NAME”>

SECTION 1 (REQUIRED)
	<ParameterList name=“InputParameters”>
		…
	</ParameterList>

SECTION 2 (REQUIRED)
	<ParameterList name=“Problem 1”>
		…
	<ParameterList name=“InputAdapterParameters”>
		…
	</ParameterList>

	<ParameterList name=“Zoltan2Parameters”>
		…
	</ParameterList>

	<ParameterList name=“Metrics”>
			<ParameterList name=“metic name”>
				…
			</ParameterList>
	</ParameterList>
	</ParameterList>

SECTION 3 (OPTIONAL)
	<ParameterList name=“Comparison”>
		…
	</ParameterList>

	<ParameterList name=“Comparison”>
		…
	</ParameterList>

</ParameterList>


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Example 1: input from file
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  <ParameterList name="InputParameters">
    <Parameter name="input path" type="string" value="PATH/TO/INPUT/DIRECTORY"/>
    <Parameter name="input file" type="string" value="FILE NAME (NO EXTENSION)"/>
    <Parameter name="file type" type="string" value="INPUT FILE TYPE"/>
  </ParameterList>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Example 2: Galari generated input
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

<ParameterList name="InputParameters">
    <Parameter name="x" type="int" value="##"/>
    <Parameter name="y" type="int" value="##"/>
    <Parameter name="z" type="int" value="##"/>
    <Parameter name="equation type" type="string" value="GALERI EQUATION"/>
  </ParameterList>


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Example 3a: A Partitioning problem definition block.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  <ParameterList name="TEST/PROBLEM TITLE #1">
    
    <!--####################################################
     Specify the problem type
     #####################################################-->
    
    <Parameter name="kind" type="string" value="PROBLEM TYPE"/>
    
    
    <!--####################################################
     Define block for the input adapter
     * must define a data type
     ** multivector data types require you to
     define a 'vector_dimension' parameter,
     which is an int value corresponding to
     the multivector dimension
     * must define an adapter type
     #####################################################-->
    
    <ParameterList name="InputAdapterParameters">
      <Parameter name="data type" type="string" value="INPUT DATA TYPE"/>
      <Parameter name="input adapter" type="string" value="INPUT ADAPTER TYPE"/>
    </ParameterList>
    
    
    <!--####################################################
     Define block of Zoltan2 problem parameters
     * all Zoltan2 parameters are valid
     * tell Zoltan to compute metrics if you are
     going to run a pass fail test defined
     by the following metrics block
     #####################################################-->
    
    <ParameterList name="Zoltan2Parameters">
      <Parameter name="algorithm" type="string" value="SPECIFY ALGORITHM" />
      <Parameter name="bisection_num_test_cuts" type="int" value="###" />
      <Parameter name="rectilinear" type="bool" value=“true/false”/>
      <Parameter name="compute_metrics" type=“bool” value="true/false”/>
    </ParameterList>
    
    
    <!--####################################################
     (OPTIONAL) Define block of metric tolerances
     #####################################################-->
    
    <ParameterList name="Metrics">
      <ParameterList name="metriccheck1">
        <Parameter name="check" type="string" value="imbalance"/>
        <Parameter name="lower" type="double" value="0.99"/>
        <Parameter name="upper" type="double" value="1.4"/>
      </ParameterList>
    </ParameterList>
    
  </ParameterList>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Example 4. Metric Definition for ‘object count’ for ImbalanceMetrics
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    <ParameterList name="Metrics">
      <ParameterList name="metriccheck1">
        <Parameter name="check" type="string" value="imbalance"/>
        <Parameter name="lower" type="double" value="0.99"/>
        <Parameter name="upper" type="double" value="1.4"/>
      </ParameterList>
    </ParameterList>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Example 5. Metric Definition for ‘total edge cuts’ for GraphMetrics
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    <ParameterList name="Metrics">
      <ParameterList name="metriccheck1">
        <Parameter name="check" type="string" value=“total edge cuts”/>
        <Parameter name="lower" type="double" value="0.99"/>
        <Parameter name="upper" type="double" value="1.4"/>
      </ParameterList>
    </ParameterList>
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Example 6.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  <ParameterList name="Comparison">
    <Parameter name="A" type="string" value="TEST/PROBLEM TITLE #1"/>
    <Parameter name="B" type="string" value="TEST/PROBLEM TITLE #2"/>
  </ParameterList>
  
  <ParameterList name="Comparison">
    <Parameter name="A" type="string" value="TEST/PROBLEM TITLE #1"/>
    <Parameter name="B" type="string" value="TEST/PROBLEM TITLE #3"/>
  </ParameterList>


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Example 7.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  <ParameterList name=“Comparison">
    <Parameter name="Problem" type="string" value="TEST/PROBLEM TITLE #1"/>
    <Parameter name="Reference" type="string" value="TEST/PROBLEM TITLE #2"/>
    
    <ParameterList name="Metrics">
      <ParameterList name="metriccheck1">
        <Parameter name="check" type="string" value=“total edge cuts”/>
        <Parameter name="lower" type="double" value=“0.5”/>
        <Parameter name="upper" type="double" value=“1.5”/>
      </ParameterList>
    </ParameterList>
    
    <ParameterList name="TIMER NAME">
      <Parameter name="lower" type="double" value="##.####"/>
      <Parameter name="upper" type="double" value="##.####"/>
    </ParameterList>
    
  </ParameterList>
