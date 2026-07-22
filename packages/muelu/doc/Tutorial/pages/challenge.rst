=============================
Challenge: elasticity example
=============================

Practical example
=================

For the second challenge, we consider an 2D elasticity example with :math:`7020` degrees of freedom. No further information is provided (geometry, discretization technique, ...).

User-interface
==============

Run the **hands-on.sh** script and choose the option 5 for the elasticity. The script automatically generates a XML file with reference multigrid parameters which are far from being optimal. The resulting problem matrix is symmetric. Therefore, we can use a CG method as outer linear solver.

.. admonition:: Exercise 1

    Open the **Elasticity2D_parameters.xml** file by pressing option 3. Try to find optimized multigrid settings using your knowledge from the previous tutorials. We have 2 (displacement) degrees of freedom per node and 3 vectors describing the near null space components (rigid body modes). All this information is automatically set correctly by the **hands-on.py** script.
    Run the example. Check the screen output (using option 1) and verify **blockdim=2** on level 1 and **blockdim=3** on level 2.
    

In the screen output of the **CoalesceDropFactory** the **blockdim** denotes the number of degrees of freedom per node (or super node on the coarser levels). Since the number of near null space vectors differs from the number of PDE equations, the number of degrees of freedom per node changes on the different multigrid levels.

.. admonition:: Exercise 2

    Open the XML parameter file (choose option 3) and try to find optimized settings. Use the advanced XML file format. Save the file, rerun the example (option 0) and compare the output with the reference results.
  

.. note::

    Use :ref:`cd_example/generalhints` for a general step-by-step procedure to optimize the multigrid parameters.

.. admonition:: Exercise 3

    How do the reference settings and your XML parameter settings perform when increasing the number of processors?


.. admonition:: Exercise 3

    Compare the results of the reference method and your preconditioner parameters when changing to a GMRES solver (instead of CG). What is changing? What about the solver timings?
    
