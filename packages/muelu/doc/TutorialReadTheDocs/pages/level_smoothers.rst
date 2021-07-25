===============
Level smoothers
===============

From the last tutorial we have learned that the used multigrid algorithm may have a significant influence in the convergence speed. When comparing the error plots for the standalone multigrid smoothers with unsmoothed and smoothed aggregation multigrid one finds also a notable difference in the ``smoothness`` of the error.

Background on multigrid methods
===============================

Obviously there are cases where some highly oscillatory error modes are left and overlaying some low frequency modes. In other cases there are only low frequency error modes left. Theses are basically the two typical cases one might find in practice.

Multigrid methods are based on the fact, that (cheap) level smoothing method often are able to smooth out high oscillatory error components whereas they cannot reduce low frequency error components very well. These low frequency error components then are transferred to a coarse level where they can be seen as high frequency error component for a level smoother on the coarse level.

One should not forget that for an efficient multigrid method both the so-called coarse level correction method and the level smoothers have to work together. That is, one has to choose the right multigrid method (e.g., ``unsmoothed`` or ``sa``) in combination with an appropriate level smoothing strategy.

In context of multigrid level smoothers we have to define both the level smoothers and the coarse solver. Usually, a direct solver is used as coarse solver that is applied to the coarsest multigrid levels. However, it is also possible to apply any other kind of iterative smoothing method or even no solver at all (even though this would be non-standard). The following XML file shows how to use a Jacobi smoother both for level smoothing and as coarse solver.

.. literalinclude:: ../../../test/tutorial/s1_easy_jacobi.xml
  :language: xml

The corresponding multigrid hierarchy is

.. warning::

    Insert missing output and ref to figures

Figures \ref{fig:1vcycles} and \ref{fig:5vcycles} show the multigrid effect of different number of Jacobi smoothers on all multigrid levels.

One has even more fine-grained control over pre- and post-smoothing.

.. literalinclude:: ../../../test/tutorial/s1_easy_jacobi2.xml
  :language: xml

This produces the following multigrid hierarchy

.. warning::

    Insert missing output

.. note::

    Note that the relaxation based methods provided by the Ifpack package are embedded in an outer additive Schwarz method.


Of course, there exist other smoother methods such as polynomial smoothers (Chebyshev) and ILU based methods.
A detailed overview of the different available smoothers can be found in the Muelu users guide ([1]_).

.. warning::

    Figure out how to insert following figures correctly

``\begin{figure}
\subfigure[1 level with 1 Jacobi sweep ($\omega=0.9$)\label{fig:2dlap111}]{\includegraphics[width=0.3\textwidth]{pics/1level_1jac09.png}}\hspace{0.03\textwidth}
\subfigure[1 level with 10 Jacobi sweeps ($\omega=0.9$)\label{fig:2dlap1110}]{\includegraphics[width=0.3\textwidth]{pics/1level_10jac09.png}}\hspace{0.03\textwidth}
\subfigure[1 level with 100 Jacobi sweeps ($\omega=0.9$)\label{fig:2dlap11100}]{\includegraphics[width=0.3\textwidth]{pics/1level_100jac09.png}} \\
\subfigure[2 level with 1 Jacobi sweep ($\omega=0.9$)\label{fig:2dlap121}]{\includegraphics[width=0.3\textwidth]{pics/2level_1jac09.png}}\hspace{0.03\textwidth}
\subfigure[2 level with 10 Jacobi sweeps ($\omega=0.9$)\label{fig:2dlap1210}]{\includegraphics[width=0.3\textwidth]{pics/2level_10jac09.png}}\hspace{0.03\textwidth}
\subfigure[2 level with 100 Jacobi sweeps ($\omega=0.9$)\label{fig:2dlap12100}]{\includegraphics[width=0.3\textwidth]{pics/2level_100jac09.png}} \\
\subfigure[3 level with 1 Jacobi sweep ($\omega=0.9$)\label{fig:2dlap131}]{\includegraphics[width=0.3\textwidth]{pics/3level_1jac09.png}}\hspace{0.03\textwidth}
\subfigure[3 level with 10 Jacobi sweeps ($\omega=0.9$)\label{fig:2dlap1310}]{\includegraphics[width=0.3\textwidth]{pics/3level_10jac09.png}}\hspace{0.03\textwidth}
\subfigure[3 level with 100 Jacobi sweeps ($\omega=0.9$)\label{fig:2dlap13100}]{\includegraphics[width=0.3\textwidth]{pics/3level_100jac09.png}} \\
\caption{2D Laplace equation on $50\times 50$ mesh after 1 V-cycle with an AMG multigrid solver and Jacobi smoothers on all multigrid levels. (2 processors)}
\label{fig:1vcycles}
\end{figure}``


``\begin{figure}
\subfigure[1 level with 1 Jacobi sweep ($\omega=0.9$)\label{fig:2dlap511}]{\includegraphics[width=0.3\textwidth]{pics/5sweeps_1level_1jac09.png}}\hspace{0.03\textwidth}
\subfigure[1 level with 10 Jacobi sweeps ($\omega=0.9$)\label{fig:2dlap5110}]{\includegraphics[width=0.3\textwidth]{pics/5sweeps_1level_10jac09.png}}\hspace{0.03\textwidth}
\subfigure[1 level with 100 Jacobi sweeps ($\omega=0.9$)\label{fig:2dlap51100}]{\includegraphics[width=0.3\textwidth]{pics/5sweeps_1level_100jac09.png}} \\
\subfigure[2 level with 1 Jacobi sweep ($\omega=0.9$)\label{fig:2dlap521}]{\includegraphics[width=0.3\textwidth]{pics/5sweeps_2level_1jac09.png}}\hspace{0.03\textwidth}
\subfigure[2 level with 10 Jacobi sweeps ($\omega=0.9$)\label{fig:2dlap5210}]{\includegraphics[width=0.3\textwidth]{pics/5sweeps_2level_10jac09.png}}\hspace{0.03\textwidth}
\subfigure[2 level with 100 Jacobi sweeps ($\omega=0.9$)\label{fig:2dlap52100}]{\includegraphics[width=0.3\textwidth]{pics/5sweeps_2level_100jac09.png}} \\
\subfigure[3 level with 1 Jacobi sweep ($\omega=0.9$)\label{fig:2dlap531}]{\includegraphics[width=0.3\textwidth]{pics/5sweeps_3level_1jac09.png}}\hspace{0.03\textwidth}
\subfigure[3 level with 10 Jacobi sweeps ($\omega=0.9$)\label{fig:2dlap5310}]{\includegraphics[width=0.3\textwidth]{pics/5sweeps_3level_10jac09.png}}\hspace{0.03\textwidth}
\subfigure[3 level with 100 Jacobi sweeps ($\omega=0.9$)\label{fig:2dlap53100}]{\includegraphics[width=0.3\textwidth]{pics/5sweeps_3level_100jac09.png}} \\
\caption{2D Laplace equation on $50\times 50$ mesh after 5 V-cycle with an AMG multigrid solver and Jacobi smoothers on all multigrid levels. (2 processors)}
\label{fig:5vcycles}
\end{figure}``

.. admonition:: Exercise 1

    Play around with the smoother parameters and study their effect on the error plot and the convergence of the preconditioned cg method. For all available smoothing options and parameters refer to the Muelu user guide ([1]_). Hint: use ``unsmoothed`` transfer operator basis functions (i.e., ``multigrid algorithm = unsmoothed``) to highlight the effect of the level smoothers.

.. admonition:: Exercise 2

    Use the following parameters to solve the :math:`50\times 50` Laplace 2D problem on 2 processors ``../../../test/tutorial/s1_easy_exercise.xml``. That is, we change to smoothed aggregation AMG. You can find the xml file also in ``../../../test/tutorial/s1_easy_exercise.xml``. Run the example on 2 processors and check the number of linear iterations and the solver timings in the screen output. Can you find smoother parameters which reduce the number of iterations? Can you find smoother parameters which reduce the iteration timings?


Footnotes
=========
.. [1] A. Prokopenko, T.A. Wiesner and J.J. Hu, MueLu User’s Guide., SAND????-????, 2014

.. warning::

    Complete source [1]!