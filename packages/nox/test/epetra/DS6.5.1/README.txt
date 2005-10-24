The following output is from a serial run on linux.  I'm keeping this
around for comparisons during refactors.


minasmorgul 12:04 > ./DS6_5_1.exe -v
Starting epetra/DS6.5.1/DS_6_5_1.exe
Serial Run

************************************************************************

-- Parameters Passed to Nonlinear Solver --

     Direction ->
       Method = "Newton"
       Newton ->
         Forcing Term Method = "Constant"
         Linear Solver ->
           Aztec Solver = "GMRES"
           Compute Scaling Manually = true   [default]
           Max Age Of Prec = 1   [default]
           Max Iterations = 20   [unused]
           Output Frequency = 1
           Output Solver Details = true   [default]
           Preconditioner = "None"
           Preconditioner Operator = "Use Jacobian"   [default]
           Tolerance = 0.0001
           Zero Initial Guess = false   [default]
         Rescue Bad Newton Solve = true   [default]
     Line Search ->
       Method = "Polynomial"
       Polynomial ->
         Allowed Relative Increase = 100   [default]
         Alpha Factor = 0.0001   [default]
         Default Step = 1   [default]
         Force Interpolation = false   [default]
         Interpolation Type = "Cubic"   [default]
         Max Bounds Factor = 0.5   [default]
         Max Iters = 100   [default]
         Maximum Iteration for Increase = 0   [default]
         Min Bounds Factor = 0.1   [default]
         Minimum Step = 1e-12   [default]
         Recovery Step = 1   [default]
         Recovery Step Type = "Constant"   [default]
         Sufficient Decrease Condition = "Armijo-Goldstein"   [default]
         Use Counters = true   [default]
     Nonlinear Solver = "Line Search Based"
     Printing ->
       MyPID = 0
       Output Information = 191
       Output Precision = 5
       Output Processor = 0
     Solver Options ->
       Status Test Check Type = 1   [default]
************************************************************************
-- Status Test Results --
**...........OR Combination ->
  **...........F-Norm = 1.699e+00 < 1.000e-06
               (Length-Scaled Two-Norm, Absolute Tolerance)
  **...........Number of Iterations = 0 < 25
************************************************************************

************************************************************************
-- Nonlinear Solver Step 0 --
f = 2.40284e+00  step = 0.00000e+00  dx = 0.00000e+00
************************************************************************

       CALCULATING FORCING TERM
       Method: Constant
       Forcing Term: 0.0001

                *******************************************************
                ***** Problem: Epetra::CrsMatrix
                ***** Preconditioned GMRES solution
                ***** No preconditioning
                ***** No scaling
                *******************************************************

                iter:    0           residual = 1.000000e+00
                iter:    1           residual = 2.400435e-01
                iter:    2           residual = 2.720359e-32


                Solution time: 0.000000 (sec.)
                total iterations: 2

************************************************************************
-- Polynomial Line Search --
       Merit Function = Sum Of Squares (default): 0.5 * ||F|| * ||F||
  1: step = 1.00000e+00 oldf = 2.40284e+00 newf = 1.07586e+03
  2: step = 1.00000e-01 oldf = 2.40284e+00 newf = 4.44026e+00
  3: step = 5.00000e-02 oldf = 2.40284e+00 newf = 2.72730e+00
  4: step = 1.16098e-02 oldf = 2.40284e+00 newf = 2.39590e+00 (STEP ACCEPTED!)
************************************************************************
************************************************************************
-- Status Test Results --
**...........OR Combination ->
  **...........F-Norm = 1.694e+00 < 1.000e-06
               (Length-Scaled Two-Norm, Absolute Tolerance)
  **...........Number of Iterations = 1 < 25
************************************************************************

************************************************************************
-- Nonlinear Solver Step 1 --
f = 2.39590e+00  step = 1.16098e-02  dx = 1.01874e+01
************************************************************************

       CALCULATING FORCING TERM
       Method: Constant
       Forcing Term: 0.0001

                *******************************************************
                ***** Problem: Epetra::CrsMatrix
                ***** Preconditioned GMRES solution
                ***** No preconditioning
                ***** No scaling
                *******************************************************

                iter:    0           residual = 1.000000e+00
                iter:    1           residual = 2.396208e-01
                iter:    2           residual = 1.132705e-32


                Solution time: 0.000000 (sec.)
                total iterations: 2

************************************************************************
-- Polynomial Line Search --
       Merit Function = Sum Of Squares (default): 0.5 * ||F|| * ||F||
  1: step = 1.00000e+00 oldf = 2.39590e+00 newf = 1.90226e+01
  2: step = 1.00000e-01 oldf = 2.39590e+00 newf = 2.24959e+00 (STEP ACCEPTED!)
************************************************************************
************************************************************************
-- Status Test Results --
**...........OR Combination ->
  **...........F-Norm = 1.591e+00 < 1.000e-06
               (Length-Scaled Two-Norm, Absolute Tolerance)
  **...........Number of Iterations = 2 < 25
************************************************************************

************************************************************************
-- Nonlinear Solver Step 2 --
f = 2.24959e+00  step = 1.00000e-01  dx = 2.40188e+00
************************************************************************

       CALCULATING FORCING TERM
       Method: Constant
       Forcing Term: 0.0001

                *******************************************************
                ***** Problem: Epetra::CrsMatrix
                ***** Preconditioned GMRES solution
                ***** No preconditioning
                ***** No scaling
                *******************************************************

                iter:    0           residual = 1.000000e+00
                iter:    1           residual = 2.254165e-01
                iter:    2           residual = 5.722776e-33


                Solution time: 0.000000 (sec.)
                total iterations: 2

************************************************************************
-- Polynomial Line Search --
       Merit Function = Sum Of Squares (default): 0.5 * ||F|| * ||F||
  1: step = 1.00000e+00 oldf = 2.24959e+00 newf = 1.31952e+00 (STEP ACCEPTED!)
************************************************************************
************************************************************************
-- Status Test Results --
**...........OR Combination ->
  **...........F-Norm = 9.330e-01 < 1.000e-06
               (Length-Scaled Two-Norm, Absolute Tolerance)
  **...........Number of Iterations = 3 < 25
************************************************************************

************************************************************************
-- Nonlinear Solver Step 3 --
f = 1.31952e+00  step = 1.00000e+00  dx = 8.73166e-01
************************************************************************

       CALCULATING FORCING TERM
       Method: Constant
       Forcing Term: 0.0001

                *******************************************************
                ***** Problem: Epetra::CrsMatrix
                ***** Preconditioned GMRES solution
                ***** No preconditioning
                ***** No scaling
                *******************************************************

                iter:    0           residual = 1.000000e+00
                iter:    1           residual = 1.807388e-02
                iter:    2           residual = 1.185009e-34


                Solution time: 0.000000 (sec.)
                total iterations: 2

************************************************************************
-- Polynomial Line Search --
       Merit Function = Sum Of Squares (default): 0.5 * ||F|| * ||F||
  1: step = 1.00000e+00 oldf = 1.31952e+00 newf = 1.59385e-01 (STEP ACCEPTED!)
************************************************************************
************************************************************************
-- Status Test Results --
**...........OR Combination ->
  **...........F-Norm = 1.127e-01 < 1.000e-06
               (Length-Scaled Two-Norm, Absolute Tolerance)
  **...........Number of Iterations = 4 < 25
************************************************************************

************************************************************************
-- Nonlinear Solver Step 4 --
f = 1.59385e-01  step = 1.00000e+00  dx = 2.32834e-01
************************************************************************

       CALCULATING FORCING TERM
       Method: Constant
       Forcing Term: 0.0001

                *******************************************************
                ***** Problem: Epetra::CrsMatrix
                ***** Preconditioned GMRES solution
                ***** No preconditioning
                ***** No scaling
                *******************************************************

                iter:    0           residual = 1.000000e+00
                iter:    1           residual = 2.982063e-01
                iter:    2           residual = 2.103038e-34


                Solution time: 0.000000 (sec.)
                total iterations: 2

************************************************************************
-- Polynomial Line Search --
       Merit Function = Sum Of Squares (default): 0.5 * ||F|| * ||F||
  1: step = 1.00000e+00 oldf = 1.59385e-01 newf = 1.01081e-02 (STEP ACCEPTED!)
************************************************************************
************************************************************************
-- Status Test Results --
**...........OR Combination ->
  **...........F-Norm = 7.147e-03 < 1.000e-06
               (Length-Scaled Two-Norm, Absolute Tolerance)
  **...........Number of Iterations = 5 < 25
************************************************************************

************************************************************************
-- Nonlinear Solver Step 5 --
f = 1.01081e-02  step = 1.00000e+00  dx = 6.15300e-02
************************************************************************

       CALCULATING FORCING TERM
       Method: Constant
       Forcing Term: 0.0001

                *******************************************************
                ***** Problem: Epetra::CrsMatrix
                ***** Preconditioned GMRES solution
                ***** No preconditioning
                ***** No scaling
                *******************************************************

                iter:    0           residual = 1.000000e+00
                iter:    1           residual = 2.991114e-01
                iter:    2           residual = 4.063169e-33


                Solution time: 0.000000 (sec.)
                total iterations: 2

************************************************************************
-- Polynomial Line Search --
       Merit Function = Sum Of Squares (default): 0.5 * ||F|| * ||F||
  1: step = 1.00000e+00 oldf = 1.01081e-02 newf = 4.62256e-05 (STEP ACCEPTED!)
************************************************************************
************************************************************************
-- Status Test Results --
**...........OR Combination ->
  **...........F-Norm = 3.269e-05 < 1.000e-06
               (Length-Scaled Two-Norm, Absolute Tolerance)
  **...........Number of Iterations = 6 < 25
************************************************************************

************************************************************************
-- Nonlinear Solver Step 6 --
f = 4.62256e-05  step = 1.00000e+00  dx = 4.13207e-03
************************************************************************

       CALCULATING FORCING TERM
       Method: Constant
       Forcing Term: 0.0001

                *******************************************************
                ***** Problem: Epetra::CrsMatrix
                ***** Preconditioned GMRES solution
                ***** No preconditioning
                ***** No scaling
                *******************************************************

                iter:    0           residual = 1.000000e+00
                iter:    1           residual = 3.051246e-01
                iter:    2           residual = 7.398888e-34


                Solution time: 0.000000 (sec.)
                total iterations: 2

************************************************************************
-- Polynomial Line Search --
       Merit Function = Sum Of Squares (default): 0.5 * ||F|| * ||F||
  1: step = 1.00000e+00 oldf = 4.62256e-05 newf = 9.97767e-10 (STEP ACCEPTED!)
************************************************************************

************************************************************************
-- Nonlinear Solver Step 7 --
f = 9.97767e-10  step = 1.00000e+00  dx = 1.92709e-05 (Converged!)
************************************************************************

************************************************************************
-- Final Status Test Results --
Converged....OR Combination ->
  Converged....F-Norm = 7.055e-10 < 1.000e-06
               (Length-Scaled Two-Norm, Absolute Tolerance)
  ??...........Number of Iterations = -1 < 25
************************************************************************

Final Parameters
****************
Direction ->
  Method = "Newton"
  Newton ->
    Forcing Term Method = "Constant"
    Linear Solver ->
      Aztec Solver = "GMRES"
      Compute Scaling Manually = true   [default]
      Max Age Of Prec = 1   [default]
      Max Iterations = 20
      Output ->
        Achieved Tolerance = 1.64e-16   [unused]
        Number of Linear Iterations = 2   [unused]
        Total Number of Linear Iterations = 14   [unused]
      Output Frequency = 1
      Output Solver Details = true   [default]
      Preconditioner = "None"
      Preconditioner Operator = "Use Jacobian"   [default]
      Tolerance = 0.0001
      Zero Initial Guess = false   [default]
    Rescue Bad Newton Solve = true   [default]
Line Search ->
  Adjusted Tolerance = 0   [unused]
  Method = "Polynomial"
  Output ->
    Total Number of Failed Line Searches = 0   [unused]
    Total Number of Line Search Calls = 7   [unused]
    Total Number of Line Search Inner Iterations = 4   [unused]
    Total Number of Non-trivial Line Searches = 2   [unused]
  Polynomial ->
    Allowed Relative Increase = 100   [default]
    Alpha Factor = 0.0001   [default]
    Default Step = 1   [default]
    Force Interpolation = false   [default]
    Interpolation Type = "Cubic"   [default]
    Max Bounds Factor = 0.5   [default]
    Max Iters = 100   [default]
    Maximum Iteration for Increase = 0   [default]
    Min Bounds Factor = 0.1   [default]
    Minimum Step = 1e-12   [default]
    Recovery Step = 1   [default]
    Recovery Step Type = "Constant"   [default]
    Sufficient Decrease Condition = "Armijo-Goldstein"   [default]
    Use Counters = true   [default]
Nonlinear Solver = "Line Search Based"
Output ->
  2-Norm of Residual = 9.98e-10   [unused]
  Nonlinear Iterations = 7   [unused]
Printing ->
  MyPID = 0
  Output Information = 191
  Output Precision = 5
  Output Processor = 0
Solver Options ->
  Status Test Check Type = 1   [default]

Test passed!
minasmorgul 12:04 >
