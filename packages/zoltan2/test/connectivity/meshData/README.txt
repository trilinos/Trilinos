Status:  1/22/20

File formats:

*.mtx
MatrixMarket files that represent the input ice sheet mesh as
a matrix. We read these files with UserInputForTests and convert them
to CrsGraphs.

*-basal-friction.ascii
line 1:  #meshvertices
For real test cases:
  remaining lines: nonzero floating point number if vertex is grounded, 0 otherwise
For manually created test cases:
  remaining lines:  1 if vertex is grounded; 0 otherwise

*-borders.quad.msh
line 1:  two meaningless words
line 2:  #meshvertices (<== always zero) #elements (<== always zero) #borderEdges 
line 3:  must be blank
remaining lines:  two mesh vertices contained within each borderEdge and
         meaningless flag (<== always zero)

*-s?-answers.txt
file of mesh vertex IDs for vertices that should be removed by the algorithm
"s?" indicates sensitivity:

   Given sensitivy s, 
   for each vertex v flagged as grounded in *-basal-friction.ascii, 
   if any element incident to v has at least s grounded vertices, v is grounded.

The current test driver cannot filter the grounding in this manner, so all answer
files being used are sensitivity 1, meaning the initial grounding information is untouched.
Any answer file of the form *-answers.txt is sensitivity 1, and equivalent to *-s1-answers.txt.

