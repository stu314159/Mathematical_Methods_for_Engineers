#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_faces.h>
#include <deal.II/numerics/data_postprocessor.h>

#include <iostream>
#include <fstream>
#include <cmath>

using namespace dealii;

const double _A = 4.0; // m, extent in the x-direction
const double _B = 3.0; // m, extent in the y-direction
const double _K = 1.0; // W/m-K, thermal conductivity

template <int dim>
class LaplaceCoeff : public Function<dim>
{
public:
  virtual double value (const Point<dim> &p,
                        const unsigned int component = 0) const override;
};

template<int dim>
double LaplaceCoeff<dim>::value(const Point<dim> & /*p*/,
                      const unsigned int /*component*/) const
{
  
  return _K;
}

template <int dim>
class InitialValues : public Function<dim>
{
public:
  virtual double value(const Point<dim> &p,
                 const unsigned int component = 0) const override;

};

template <int dim>
double InitialValues<dim>::value(const Point<dim> &p,
                           const unsigned int /*component*/) const
{
  double _x = p[0];
  double _y = p[1];
  
 /* if ((_x < _A/2.) && (_y < _B/2.))
  {
    return -1.0;
  } else
  {
    return 1.0;
  }*/
  
  return _x*(_A-_x)*_y*(_B-_y);

}

class Lec29a
{
public:
  Lec29a ();
  void run ();
  
private:
  void make_grid ();
  void setup_system ();
  void solve_time_step ();
  void output_results () const;
  
  Triangulation<2> 		triangulation;
  FE_Q<2>			fe;
  DoFHandler<2>			dof_handler;
  AffineConstraints<double>	constraints;
  
  SparsityPattern		sparsity_pattern;
  SparseMatrix<double>		mass_matrix;
  SparseMatrix<double>		laplace_matrix;
  SparseMatrix<double>		system_matrix;
  
  Vector<double>		solution;
  Vector<double>		old_solution;
  Vector<double>		system_rhs;
  
  double			time;
  double			time_step; // actually, step length
  unsigned int			timestep_number;
  
  const double			theta;
  const double			max_time;

};

Lec29a::Lec29a()
        : fe(2) // use 2nd-order polynomial shape functions
        , dof_handler(triangulation)
        , time_step(1./50.)
        , theta(0.5) // theta = 0.5 corresponds to the Newmark method
	, max_time(2.0)
{}

void
Lec29a::make_grid()
{
  const Point<2> lower_left(0,0);
  const Point<2> upper_right(_A,_B);
  GridGenerator::hyper_rectangle(triangulation, lower_left,
                            upper_right, false);
  
  triangulation.refine_global(6);
  
  std::cout << "Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl;
                              
}

void 
Lec29a::setup_system()
{
  dof_handler.distribute_dofs(fe);
  constraints.clear();
  DoFTools::make_hanging_node_constraints(dof_handler,constraints);
  // I don't think the above line is needed since I will not be
  // using AMR and thus there should not be any hanging nodes.
  // Intentionally leaving anyway until I really figure the whole
  // thing out.
  constraints.close();
  
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler,dsp);
  sparsity_pattern.copy_from(dsp);
  
  mass_matrix.reinit(sparsity_pattern);
  laplace_matrix.reinit(sparsity_pattern);
  system_matrix.reinit(sparsity_pattern);
  
  MatrixCreator::create_mass_matrix(dof_handler,
                                    QGauss<2>(fe.degree + 1),
                                    mass_matrix);
  
  LaplaceCoeff<2> my_coeff;
  MatrixCreator::create_laplace_matrix(dof_handler,
                                       QGauss<2>(fe.degree + 1),
                                       laplace_matrix,
                                       &my_coeff);
  
  solution.reinit(dof_handler.n_dofs());
  old_solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());                                    
}

void
Lec29a::solve_time_step()
{
SolverControl			solver_control(1000,1e-8*system_rhs.l2_norm());
SolverCG<Vector<double>> 	cg(solver_control);
PreconditionSSOR<SparseMatrix<double>> preconditioner;

preconditioner.initialize(system_matrix,1.0);
cg.solve(system_matrix, solution, system_rhs, preconditioner);

constraints.distribute(solution);

std::cout << "  " << solver_control.last_step() << " CG iterations."
          << std::endl;
}

void
Lec29a::output_results() const
{
  DataOut<2> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution,"U");
  data_out.build_patches();
  data_out.set_flags(DataOutBase::VtkFlags(time,timestep_number));
  
  const std::string filename = 
     "solution-" + Utilities::int_to_string(timestep_number,3) + ".vtk";
  std::ofstream output(filename);
  data_out.write_vtk(output);
}

void
Lec29a::run()
{
  time = 0.0;
  timestep_number = 0;
  
  make_grid();
  setup_system();
  
  std::cout << "Number of active cells: " 
            << triangulation.n_active_cells()
            << std::endl;
  
  std::cout << "Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << std::endl;
            
  VectorTools::interpolate(dof_handler,
                           InitialValues<2>(),
                           old_solution);
                           
  Vector<double> tmp;
  
  solution = old_solution;
  output_results();
  system_rhs = 0.;
  
  tmp.reinit(solution.size());
  
  while (time <= max_time)
  {
    time += time_step;
    ++timestep_number;
    std::cout << "Time step " << timestep_number
              << " at t=" << time
              << std::endl;
              
    // MU^{n-1}-(1-theta)k_n A U^{n-1} is the system rhs
    // construct system rhs below
    mass_matrix.vmult(system_rhs,old_solution); // M*U^{n-1} --> system_rhs
    laplace_matrix.vmult(tmp,old_solution); // A*U^{n-1} -> tmp
    system_rhs.add(-(1.0-theta)*time_step,tmp); // system_rhs = -(1-theta)k_n (tmp) -> system_rhs
    system_matrix.copy_from(mass_matrix); // SM = M
    system_matrix.add(theta*time_step,laplace_matrix); // SM = M+k_n*theta*A
    
    constraints.condense(system_matrix, system_rhs);
    
    std::map<types::global_dof_index,double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<2>(),
                                             boundary_values);
    MatrixTools::apply_boundary_values(boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);
                                       
    solve_time_step();
    output_results();
    old_solution = solution;
  
  }

}

int
main()
{
  Lec29a heat_equation_solver;
  heat_equation_solver.run();
  
  return 0;

}






