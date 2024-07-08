[Mesh]
  [gmg]
    type=GeneratedMeshGenerator
    dim = 2
    nx = 100
    ny = 100
    xmin = 0.0
    xmax = 5.0
    ymin = 0.0
    ymax = 3.0
  []
[]

[Problem]
  type = FEProblem
[]

[Variables]
  [u]
    family = LAGRANGE # default shape function type
    order = FIRST  # default order
  []
[]

[Kernels]
  [diff]
    type = ADDiffusion
    variable = u
  []
[]

[Functions]
  [top_parta]
    type = ParsedFunction
    expression = 'x^2'
  []
  [top_partb]
    type = ConstantFunction
    value = 6.25
  []
  [top_fun]
    type = PiecewiseFunction
    axis = x
    axis_coordinates = '2.5'
    functions = 'top_parta top_partb'
  []
[]

[BCs]
  [./top]
    type = FunctionDirichletBC
    variable = u
    boundary = top
    function = top_fun
  [../]
  [./bottom]
    type = DirichletBC
    variable = u
    boundary = bottom
    value = 0
  [../]
  [./sides]
    type = NeumannBC
    variable = u
    boundary = 'left right'
    value = 0
  [../]
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

