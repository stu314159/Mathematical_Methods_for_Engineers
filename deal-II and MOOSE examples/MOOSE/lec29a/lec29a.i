A = 5.0
B = 3.0
k = 2.0


[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 100
    ny = 100
    xmin = 0.0
    xmax = ${A}
    ymin = 0.0
    ymax = ${B}  
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
    expression = 'x*(${A}-x)*y*(${B}-y)'
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
  [lec29a_bcs]
    type = DirichletBC
    variable = u
    boundary = 'left right top bottom'
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
  end_time = 1.0
  num_steps = 50
  scheme = implicit-euler
  verbose = true
[]

[Outputs]
  exodus = true
[]

