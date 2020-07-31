After a look and run through the debugger; I'm not sure how the nodal
variable interpolation *ever* worked in the past.  It looks totally
wrong and I can understand why it isn’t working now; I just can't
understand how it ever worked correctly for multi-element block
models.

As I see it, it does the following for the interpolation:

 * For all blocks

   * For all time steps

     * For all nodal variables
       * Iterate all nodes in this block; map from A->B
       * Write values for *all* B nodes at this step for this variable

   * For all time steps

     * For all element variables
       * Iterate all elements in this block; map from A->B
       * Write values for all elements in this block at this step for this variable

This works for element variables since the exodus API can output
elements a block and variable at a time.  For nodes, it doesn't work
and you will end up with the values at the last step for all
nodes/variables except for nodes which are only in the last element
block which seems to be what you are seeing.

Fixing this would be a major undertaking and I'm not sure it would get
prioritized (although you are welcome to try).

This *should* work OK if you only do a single timestep or if you only
have a single element block.  With a single timestep and multiple
element blocks, there is an issue of what happens if the node is
shared between multiple element blocks -- it will only get the
interpolated value from the last block.

Now, what to do...
 * I think that the Percept code can do some mapping from mesh to mesh...
 * Klugy, but can do a timestep at a time and then rejoin all timesteps using `conjoin`
 * Klugy, but can subset down to one block / mesh and then run mapvar on each submesh and then join using `ejoin`

Sorry for the bearer of bad news, but hopefully there is a path to get
done what you need...
