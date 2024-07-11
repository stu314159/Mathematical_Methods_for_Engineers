# Parameters
C = 2 # radius of the disc

[Mesh]
  [disc]
    type = AnnularMeshGenerator
    nr = 40
    nt = 40
    rmin = 0
    rmax = ${C}
  []
[]

[Problem]
  type = FEProblem
[]

[Variables]
  [u]
    family = LAGRANGE
    order = FIRST
  []
[]

[Functions]
  [bdy]
    type = Lec30aBCFunction
    #type = ParsedFunction
    #expression = 'x'
    #expression = 'atan(y/x)'  # does not seem to support atan2
  []
[]

[Kernels]
  [diff]
    type = ADDiffusion
    variable = u
  []
[]

[BCs]
  [outside]
    type = ADFunctionDirichletBC
    variable = u
    boundary = rmax
    function = bdy 
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
[]
  
