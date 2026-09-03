# Navier–Stokes Preconditioners {#teko_navier_stokes}

[TOC]

Teko's physics-based preconditioners target the (linearized) incompressible Navier–Stokes
saddle-point system

\f[
\begin{bmatrix} F & B^T \\ B & -C \end{bmatrix}
\begin{bmatrix} u \\ p \end{bmatrix}
=
\begin{bmatrix} f \\ g \end{bmatrix},
\f]

where \f$F\f$ is the (convection–diffusion) velocity operator, \f$B\f$/\f$B^T\f$ the
divergence/gradient coupling, and \f$C\f$ an optional pressure stabilization. The methods
here — **SIMPLE**, **LSC**, and **PCD** — differ in how they approximate the pressure Schur
complement \f$S = -(C + B F^{-1} B^T)\f$. For the theory see
[Cyr, Shadid & Tuminaro (JCP 2011)](http://dx.doi.org/10.1016/j.jcp.2011.09.001).

Some of these methods need operators the algebra cannot produce from the assembled system
alone (a velocity mass matrix, a pressure Laplacian). Teko obtains these through the
[RequestHandler](07-advanced.md) — the application must register a callback that answers
the corresponding request message. Each method below lists the messages it may issue.

---

## SIMPLE — `"NS SIMPLE"`

SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) approximates \f$F^{-1}\f$ by the
inverse of a diagonal of \f$F\f$ when forming the Schur complement.
(`src/NS/Teko_SIMPLEPreconditionerFactory.cpp`)

| Parameter | Type | Default | Meaning |
|-----------|------|---------|---------|
| `"Inverse Type"` | inverse name | Amesos2 if enabled, else empty | Default inverse for both sub-solves. |
| `"Inverse Velocity Type"` | inverse name | — | Override for the velocity (\f$F\f$) solve. |
| `"Preconditioner Velocity Type"` | inverse name | — | Optional preconditioner for the velocity solve. |
| `"Inverse Pressure Type"` | inverse name | — | Override for the pressure (Schur) solve. |
| `"Preconditioner Pressure Type"` | inverse name | — | Optional preconditioner for the pressure solve. |
| `"Alpha"` | double | `1.0` | SIMPLE relaxation factor. |
| `"Explicit Velocity Inverse Type"` | enum / inverse name | — | How \f$F\f$ is approximated in the Schur complement. One of the [`Diagonal Type` values](04-block-preconditioners.md) (`Diagonal`/`Lumped`/`AbsRowSum`/`BlkDiag`), or an inverse label. |
| `"Use Mass Scaling"` | bool | `false` | Scale by the velocity mass matrix (issues the `"Velocity Mass Matrix"` request). |
| `"H options"` | sublist | — | Read only when `"Explicit Velocity Inverse Type" = "BlkDiag"`. |

**RequestHandler messages:** `"Velocity Mass Matrix"` (when `"Use Mass Scaling"` is true).

**`"NS SIMPLE-Timed"`** (`src/NS/Teko_TimingsSIMPLEPreconditionerFactory.cpp`) is identical in
parameters but adds fine-grained timers.

---

## LSC — `"NS LSC"`

The Least-Squares Commutator method. The factory selects one of three Schur-complement
*strategies*. (`src/NS/Teko_LSCPreconditionerFactory.cpp`)

### Factory-level parameters

| Parameter | Type | Default | Meaning |
|-----------|------|---------|---------|
| `"Is Symmetric"` | bool | — | Whether the system is symmetric. |
| `"Strategy Name"` | string | `"Basic Inverse"` | `"Basic Inverse"`, `"Pressure Laplace"`, or `"SIMPLEC"`. |
| `"Strategy Settings"` | sublist | — | Parameters for the chosen strategy (required unless `"Basic Inverse"`). |

### Basic Inverse strategy (`InvLSCStrategy`)

The default and most general strategy. (`src/NS/Teko_InvLSCStrategy.cpp`)

| Parameter | Type | Default | Meaning |
|-----------|------|---------|---------|
| `"Inverse Type"` | inverse name | — | Default inverse for the sub-solves. |
| `"Inverse Velocity Type"` | inverse name | — | Velocity (\f$F\f$) solve. |
| `"Inverse Pressure Type"` | inverse name | — | Pressure (Schur) solve. |
| `"Ignore Boundary Rows"` | bool | — | Ignore boundary rows when forming the commutator. |
| `"Use LDU"` | bool | — | Use the full LDU form. |
| `"Use Mass Scaling"` | bool | — | Scale by the velocity mass matrix. |
| `"Use W-Scaling"` | bool | — | Apply W-scaling. |
| `"Eigen Solver Iterations"` | int | — | Iterations for the internal eigenvalue estimate. |
| `"Scaling Type"` | enum | — | Diagonal used for scaling (`Diagonal`/`Lumped`/`AbsRowSum`/`BlkDiag`). |
| `"Assume Stable Discretization"` | bool | — | Skip the pressure-stabilization term. |

**RequestHandler messages:** `"Velocity Mass Matrix"` (when `"Use Mass Scaling"`),
`"W-Scaling Vector"` (when `"Use W-Scaling"`).

### SIMPLEC strategy (`LSCSIMPLECStrategy`)

(`src/NS/Teko_LSCSIMPLECStrategy.cpp`) Accepts: `"Inverse Type"`, `"Inverse Velocity Type"`,
`"Inverse Pressure Type"`, `"Use LDU"`, `"Use Mass Scaling"`, `"Eigen Solver Iterations"`,
`"Scaling Type"`.

### Pressure Laplace strategy (`PresLaplaceLSCStrategy`)

(`src/NS/Teko_PresLaplaceLSCStrategy.cpp`) Same parameter set as SIMPLEC.
**RequestHandler messages:** `"Pressure Laplace Operator"`, `"Velocity Mass Operator"`.

---

## PCD — `"Block LU2x2"` with `"Strategy Name" = "NS PCD Strategy"`

The Pressure Convection–Diffusion method is delivered as a Schur-complement *strategy* of the
[`Block LU2x2`](04-block-preconditioners.md) factory (rather than its own top-level `"Type"`).
It approximates \f$S^{-1} \approx M_p^{-1} F_p A_p^{-1}\f$ using a pressure mass matrix, a
pressure convection–diffusion operator, and a pressure Laplacian.
(`src/NS/Teko_PCDStrategy.cpp`)

Set `"Type" = "Block LU2x2"`, `"Strategy Name" = "NS PCD Strategy"`, and put these in
`"Strategy Settings"`:

| Parameter | Type | Default | Meaning |
|-----------|------|---------|---------|
| `"Inverse Type"` | inverse name | — | Default inverse for the sub-solves. |
| `"Inverse F Type"` | inverse name | — | Inverse of the velocity operator \f$F\f$. |
| `"Inverse Laplace Type"` | inverse name | — | Inverse of the pressure Laplacian \f$A_p\f$. |
| `"Inverse Mass Type"` | enum / inverse name | — | Approximation of the pressure mass matrix \f$M_p\f$ (a [`Diagonal Type`](04-block-preconditioners.md) or an inverse label). |
| `"Flip Schur Complement Ordering"` | bool | `false` | Reverse the order of the Schur-complement factors. |
| `"Pressure Laplace Parameters"` | sublist | — | Settings forwarded to the Laplace sub-problem. |
| `"Pressure Convection Diffusion Parameters"` | sublist | — | Settings for the convection–diffusion sub-problem. |

**RequestHandler messages:** `"PCD Operator"`, `"Pressure Laplace Operator"`.

---

## Worked SIMPLE example (XML)

A SIMPLE preconditioner on a two-field (`"2 1"`) velocity–pressure system, solving the
velocity block with MueLu and the pressure block with an ILU:

```xml
<ParameterList name="Teko">
  <Parameter name="Strided Blocking" type="string" value="2 1"/>
  <Parameter name="Inverse Type"     type="string" value="SIMPLE"/>
  <ParameterList name="Inverse Factory Library">

    <ParameterList name="SIMPLE">
      <Parameter name="Type"                  type="string" value="NS SIMPLE"/>
      <Parameter name="Alpha"                 type="double" value="0.9"/>
      <Parameter name="Explicit Velocity Inverse Type" type="string" value="AbsRowSum"/>
      <Parameter name="Inverse Velocity Type" type="string" value="VelocityMG"/>
      <Parameter name="Inverse Pressure Type" type="string" value="PressureILU"/>
    </ParameterList>

    <ParameterList name="VelocityMG">
      <Parameter name="Type" type="string" value="MueLu"/>
      <Parameter name="number of equations" type="int" value="2"/>
    </ParameterList>

    <ParameterList name="PressureILU">
      <Parameter name="Type"      type="string" value="Ifpack2"/>
      <Parameter name="Prec Type" type="string" value="ILUT"/>
    </ParameterList>

  </ParameterList>
</ParameterList>
```

## Supplying operators via the RequestHandler

When a method above lists a RequestHandler message, the application must answer it. In
outline:

```cpp
class MassMatrixCallback : public Teko::RequestCallback<Teko::LinearOp> {
public:
  bool handlesRequest(const Teko::RequestMesg& m) override {
    return m.getName() == "Velocity Mass Matrix";
  }
  Teko::LinearOp request(const Teko::RequestMesg& /*m*/) override { return massOp_; }
  void preRequest(const Teko::RequestMesg&) override {}
private:
  Teko::LinearOp massOp_;
};

RCP<Teko::RequestHandler> rh = rcp(new Teko::RequestHandler);
rh->addRequestCallback(rcp(new MassMatrixCallback(/* ... */)));
// hand rh to the StratimikosFactory (setRequestHandler) or the InverseLibrary
```

See [Advanced Topics](07-advanced.md) for the full RequestHandler API.
