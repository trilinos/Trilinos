Status:  1/22/20

File formats:

*.mtx
MatrixMarket files that represent the input ice sheet mesh as
a matrix. We read these files with UserInputForTests and convert them
to CrsGraphs.

*-basal-friction.ascii
line 1:  #meshvertices
remaining lines: nonzero floating point number if vertex is grounded, 0 otherwise

*-borders.quad.msh
line 1:  two meaningless words
line 2:  #borderEdges 
remaining lines:  two mesh vertices contained within each borderEdge

*-answers.txt
A list of vertices that are expected to be removed from the input.
The format for this file is one vertex per line, with only removed 
vertices appearing in this file. 
No count of the vertices should appear in the file.
Empty answer files indicate nothing should be removed.
