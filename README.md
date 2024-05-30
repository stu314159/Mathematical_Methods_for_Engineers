# Mathematical_Methods_for_Engineers
A (draft) textbook based on two undergraduate courses; one in analytical methods and the other in numerical methods.  

## Analytical Methods
The first course is focused on analytical methods for linear ordinary and partial differential equations.  All students come into the course having taken a three-semester sequence of calculus along with a course in ordinary differential equations.  The analytical methods portion briefly reviews methods for constant coefficient linear equations and proceeds to methods for non-constant coefficient ODEs including Cauchy-Euler equations, power series methods, and the method of Frobenieus.  After a review of Fourier Series methods and an introduction to Fourier-Legendre and Fourier-Bessel expansions we thoroughly explore solutions to second-order, linear, partial differential equations.  Since many students are also studying nuclear engineering, there is a heavy focus on addressing boundary value problems in cylindrical and spherical coordinate systems that are applicable to other topics of interest such as reactor physics.  There is also heavy emphasis on heat transfer applications that students will see later on in their undergraduate curriculum.

The materials presented are based heavily on Professor Dennis Zill's excellent book, Advanced Engineering Mathematics, 7th Edition. We lightly select from chapters 1-3 for review; chapter 5 for series solution methods; and chapters 12-14 for Fourier Series and solutions to linear boundary value problems.  

What distinguishes this course from Prof Zill's work is the explicit incorporation of computational tools in the solution process.  These ``semi-analytical methods'' are presented here in MATLABowing to the students preparation with that tool.  Other open-source tools like Octave and Python, of course, could be used.

## Numerical Methods
The second course is focused on numerical methods for a variety of applications relevant to the physical sciences.  Prerequisites for this course include a three-semester sequence of calculus followed by a course in ordinary differential equations. Students also are required to have taken an introductory programming course where they have gained experience using the MATLAB programming environment.  A majority of the students who have taken this course have also taken the Analytical Methods course previously, although it is not a requirement.  

The course starts with an introduction to number representation on a computer and reviews the basics of linear algebra.  This is followed by treatment of algorithms for solving nonlinear equations and systems of equations.  For linear systems of equations, both direct and iterative methods are described; this includes treatment of sparse linear systems and preconditioning for iterative algorithms.  Least squares curve fitting is described in some detail as is interpolation with Lagrange polynomials.  Numeric differentiation is treated with both finite difference formulas as well as differentiation with Lagrange polynomials.  Several algorithms for numeric integration are described along with a fairly complete derivation of Gauss quadrature.  Several Runge-Kutta methods are illustrated for solving initial value problems, this includes a description of embedded RK methods with adaptive time-stepping.  Boundary value problems are solved using the shooting method, finite difference methods, and finite element methods.  A moderately detailed derivation of the Galerkin formulation of the method of weighted residuals is used to introduce Finite element methods.  The course is finished with several computational workshops using COMSOL, the details of which are not treated in this manuscript. 

The materials presented are based heavily on the excellent textbook by Prof Amos Gilat and Prof Vish Subramaniam, Numerical Methods for Scientists and Engineers, 3rd Edition.  The review is taken from selected portions of chapter 1 and 2; methods for solving non-linear and linear equations from chapter 3 and 4; curve fitting, interpolation, differentiation and integration are derived, in part, from chapters 6, 8, and 9; and most of the initial and boundary value problem examples are taken from chapters 10 and 11.

Several topics have been added to the treatment of Gilat and Subramaniam in an effort to selectively add depth and, in a few cases, provide updates to highlight tools available with more recent releases of MATLAB.  Depth is added with the treatment of sparse matrices and preconditioned iterative solution schemes for linear systems of equations; a more detailed derivation of Gauss quadrature is provided; several more Runge-Kutta methods are shown for solving initial value problems including implicit and embedded RK methods.  The treatment of finite element methods is based, in part, on the exposition from Finite Element Methods Using MATLAB (2nd Edition), by Prof Young Kwon.

This course is intended to serve scientists and engineers who will use numerical methods *as a tool* in the context of computational problems that they may face.  As a consequence, a fairly pragmatic approach is taken with many of the derivation details summarized or, in some cases, skipped altogether.  Many methods are described not derived and theorems are used and applied, not proven.  I hope that any mathematician with the patience to read through these pages will graciously forgive me and, if needed, look for a different reference.
