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

// this is an incredibly bad way of doing this but...
const double L = 1.0; // m, length
const double k = 0.1; // W/m-K, thermal conductivity

template <int dim>
class LaplaceCoeff : public Function<dim>
{
public:
  virtual double value(const Point<dim> &p,
                  const unsigned int component = 0) const override;
};

template <int dim>
double LaplaceCoeff<dim>::value(const Point<dim> & /*p*/,
                 const unsigned int /*component*/) const
{
  return k;
}



template <int dim>
class InitialValues : public Function<dim>
{
public:
  virtual double value(const Point<dim> &p,
             const unsigned int component = 0) const override;
};

template <int dim>
double InitialValues<dim>::value(const Point<dim> & p,
                  const unsigned int /*component*/) const
{
  double _x = p[0];
  
  return _x*(L-_x);
}

class Lec23a
{
public:
  Lec23a ();
  void run ();

private: 
  void make_grid ();
  void setup_system ();
  //void assemble_system(); // keep this since I will not use AMR 
  void solve_time_step ();
  void output_results () const;
  
  Triangulation<1> 	triangulation;
  FE_Q<1>		fe;
  DoFHandler<1>		dof_handler;
  
  AffineConstraints<double> constraints;
  
  SparsityPattern 	sparsity_pattern;
  SparseMatrix<double>	mass_matrix;
  SparseMatrix<double>	laplace_matrix;
  SparseMatrix<double>	system_matrix;
  
  Vector<double>	solution;
  Vector<double> 	old_solution;
  Vector<double>	system_rhs;
  
  double		time;
  double		time_step;
  unsigned int		timestep_number;
  
  const double 		theta;
  const double		max_time;

};

Lec23a::Lec23a ()
         : fe(1) 
         , dof_handler(triangulation)
         , time_step(1./100.)
         , theta(0.5)
         , max_time(0.5)
{}

void
Lec23a::make_grid()
{
GridGenerator::hyper_cube(triangulation,
				0., L, true); // left is 0 and right is L

triangulation.refine_global(5);
std::cout << "Number of active cells: " <<
             triangulation.n_active_cells()
          << std::endl;
          		
}

void Lec23a::setup_system()
{
  dof_handler.distribute_dofs(fe);
  
  constraints.clear();
  DoFTools::make_hanging_node_constraints(dof_handler,constraints);
  
  // include Dirichlet BCs in the constraints object
 // VectorTools::interpolate_boundary_values(dof_handler,
 //                                          0, // left boundary
 //                                          Functions::ZeroFunction<1>(),
 //                                          constraints);
 // VectorTools::interpolate_boundary_navles(dof_handler,
 //                                          1, // right boundary
 //                                          Functions::ZeroFunction<1>(),
 //                                          constraints);
                                           
  
  constraints.close();
  
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler,dsp);
  sparsity_pattern.copy_from(dsp);
  
  mass_matrix.reinit(sparsity_pattern);
  laplace_matrix.reinit(sparsity_pattern);
  system_matrix.reinit(sparsity_pattern);
  
  // cool functions to create and assemble system
  MatrixCreator::create_mass_matrix(dof_handler,
  				    QGauss<1>(fe.degree + 1),
  				    mass_matrix);
  				    
  LaplaceCoeff<1> my_coeff;
  // how to make the function bit work...				    
  MatrixCreator::create_laplace_matrix(dof_handler,
                                       QGauss<1>(fe.degree + 1),
                                       laplace_matrix,
                                       &my_coeff);
  
  solution.reinit(dof_handler.n_dofs());
  old_solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());

}

void Lec23a::solve_time_step()
{
  SolverControl		solver_control(1000,1e-8*system_rhs.l2_norm());
  SolverCG<Vector<double>> cg(solver_control);
  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  
  preconditioner.initialize(system_matrix,1.0);
  cg.solve(system_matrix,solution,system_rhs,preconditioner);
  constraints.distribute(solution);
  
  std::cout << "  " << solver_control.last_step() << " CG iterations."
            << std::endl;
  
}

void Lec23a::output_results() const
{
  DataOut<1> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution,"U");
  data_out.build_patches();
  data_out.set_flags(DataOutBase::VtkFlags(time,timestep_number));
  
  const std::string filename = 
     "solution-" + Utilities::int_to_string(timestep_number,3) + ".vtk";
  std::ofstream output(filename);
  data_out.write_vtk(output);
}

void Lec23a::run()
{
  time = 0.0;
  timestep_number = 0;
  
  std::cout << "in run function..." << std::endl;
  
  std::cout << "making the grid..." << std::endl;
  make_grid();
  
  std::cout << "setting up the system..." << std::endl;
  setup_system();
  
  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;
            
  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;
  
  // apply initial condition (problem with this function...)
  VectorTools::interpolate(dof_handler,
                           InitialValues<1>(),
                           old_solution);
                           
  std::cout << "initial condition applied... " << std::endl;
                           
  Vector<double> tmp; // used for a variety of purposes
  //Vector<double> forcing_terms; // leaving this off for this exercise.
  
  solution = old_solution;
  output_results(); // output initial solution
  
  system_rhs = 0.; // try this now...
  
  tmp.reinit(solution.size());
    
  while (time <= max_time)
  {
    time += time_step;
    ++timestep_number;
    std::cout << "Time step " << timestep_number << " at t=" << time
              << std::endl;
              
    // MU^{n-1}-(1-theta)k_n A U^{n-1} (system rhs
    
    mass_matrix.vmult(system_rhs,old_solution); // M*U^{n-1}, puts result into system_rhs 
    laplace_matrix.vmult(tmp,old_solution); // A*U^{n-1}, puts resut into tmp
    system_rhs.add(-(1-theta)*time_step,tmp); // I think it actually multiplies the two arguments and adds to system_rhs. (yes, per the documentation)
    
    system_matrix.copy_from(mass_matrix);// SM = M
    system_matrix.add(theta*time_step,laplace_matrix);// SM = M + k_n*theta*A
    
    constraints.condense(system_matrix, system_rhs);
    
    std::map<types::global_dof_index,double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<1>(),
                                             boundary_values);
    VectorTools::interpolate_boundary_values(dof_handler,
                                              1,
                                              Functions::ZeroFunction<1>(),
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

int main ()
{
  Lec23a heat_equation_solver;
  
  std::cout << "heat equation solver object created..." << std::endl;
  heat_equation_solver.run();
  
  return 0;

}






