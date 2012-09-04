<!--
   XML definition for a collection of tests to be run by the test driver.

   This format is very basic.  We could add many attributes to almost
   every section.

   "Tests" is the outer most level.  A collection of tests can have a "name".

   Within "Tests" are are a list of individual tests with tag "Test".  A
  "Test" can have a "name".

   Within a "Test" we can have "TestParameters" and "Zoltan2Parameters".
   Eventually we could add "TPLParameters".

   Within "TestParameters" we specify three things so far: the input
   data, the type of InputAdapter to use and the passing criteria.  
   For now we will just use weights and coordinates if they are found
   in "inputFile".  We should add options to generate weights or
   coordinates in patterns such that we know what imbalance would be
   a good imbalance.

   The "Zoltan2Parameters" section can contain a subsection for
   any Zoltan2 parameter.

-->

<Tests   name="rcbParallelTests" >

  <Test  name="test1" >

     <TestParameters >
        <inputFile="simple.mtx"  />

        <inputMesh  xdim=""  ydim="" zdim=""  matrixType="" />

        <inputGeometry  type=""  size=""  />

        <weights  includeWeightsFromFile="true "
                   objectWeightDimension="1 "
                   edgeWeightDimension="0 " />

        <coordinates  includeCoordinatesFromFile="false"
                      coordinateDimension="3" />


        <inputAdapter name="BasicCoordinateInput" />

        <passingCriteria >
           <imbalance  equals=""  lessThan=""  greaterThan="" />

        </passingCriteria>
        
     </TestParameters>

     <Zoltan2Parameters >
       <parameterName1   value="" />
       <parameterName2   value="" />
     </Zoltan2Parameters>

   </Test> 

</Tests>
