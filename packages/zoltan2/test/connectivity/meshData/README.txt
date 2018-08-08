Status:  8/8/18

File formats:

*.quad.msh -- mesh files with quad elements
line 1:  two meaningless words
line 2:  #meshvertices #elements #floatingBorderEdges (<== always zero)
lines 3-#meshvertices+2:  mesh vertex numbers (one based)
remaining lines:  mesh vertices contained within each element
          four vertices and meaningless flag (<== always zero)

*-basal-friction.ascii
line 1:  #meshvertices
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

NOTE:  Currently, these files are zero-based; it would be better for them to 
match the mesh files (one-based) so that we could more easily hand-check them.

***NOTE***:  
Ian is not sure he is checking case where answer file says something
should be removed by distributed does not flag it
Need to check removing same and keeping same.  (See below.)

***NOTE***:  
tictactoe answers need to be checked (particularly for s1-all -- we expected
everything to be kept in this case)

