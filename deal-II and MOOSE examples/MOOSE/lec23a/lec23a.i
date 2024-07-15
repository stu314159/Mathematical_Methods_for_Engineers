k = 0.1 # W/m-K, diffusion coefficient
L = 1.0 # m, length of the 

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 100
    xmin = 0.0
    xmax = 1.0
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
  [ic_fun]
    type = ParsedFunction
    expression = 'x*(${L}-x)'
  []
[]

[Kernels]
  [diff]
    type = CoefDiffusion
    variable = u
    coef = ${k}
  []
  
  [time]
    type = ADTimeDerivative
    variable = u
  []  
[]

[BCs]
  [left]
    type = ADDirichletBC
    variable = u
    boundary = left
    value = 0
  []
  [right]
    type = ADDirichletBC
    variable = u
    boundary = right
    value = 0
  []
[]

[ICs]
  [u_ic]
    type = FunctionIC
    variable = u
    function = ic_fun
  []
[]

[Executioner]
  type = Transient
  end_time = 0.5
  num_steps = 40
  scheme = implicit-euler
  verbose = true
[]

[Outputs]
  exodus = true
[]


