// Deal.II libraries that I will (may) use
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

// object to contain the final sparse system matrix
#include <deal.II/lac/sparse_matrix.h>

// object that will pre-compute the sparsity pattern in a 
// memory-efficient way
#include <deal.II/lac/dynamic_sparsity_pattern.h>

// solver and pre-conditioner for the sparse linear system of equations
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

// the fundamental class for creating meshes (grids)
#include <deal.II/grid/tria.h>

// some standard grid generation tools for simple geometries (like this one)
#include <deal.II/grid/grid_generator.h>

// output a graphic representation of the grid
#include <deal.II/grid/grid_out.h> 
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

// the class FE_Q handles finite elements using Lagrange polynomials 
// as shapre functions
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

//C++ libraries that I need to use
#include <iostream>
#include <fstream>
#include <cmath>

using namespace dealii;


// declare the class

class Lec26a
{
  public:
    Lec26a ();
    void run ();
    
  private:
    void make_grid ();
    void setup_system ();
    void assemble_system ();
    void solve ();
    void output_results () const;
    
    Triangulation<2> 	triangulation;
    FE_Q<2> 		fe;
    DoFHandler<2>	dof_handler;
    
    SparsityPattern	sparsity_pattern;
    SparseMatrix<double> system_matrix;
    Vector<double>	solution;
    Vector<double>	system_rhs;
    Point<2> lower_left;
    Point<2> upper_right;
};

// constructor
Lec26a::Lec26a()
  : fe(2), dof_handler(triangulation), 
    lower_left(0,0), upper_right(5,3)
{}

// make_grid member function
void Lec26a::make_grid()
{
  GridGenerator::hyper_rectangle(triangulation, lower_left, upper_right);
  triangulation.refine_global(6);
  
  std::cout << "  Number of active cells: " << triangulation.n_active_cells()
            << std::endl
            << "  Total number of cells: " << triangulation.n_cells()
            << std::endl;
            
  const bool OUTPUT_GRID = true;
  if(OUTPUT_GRID)
  {
    std::ofstream 	out("grid.vtk");
    GridOut  		grid_out;
    grid_out.write_vtk(triangulation,out);
    std::cout << "Grid written to grid.vtk" << std::endl;
  }
}

// allocate memory for the linear system of equations
void Lec26a::setup_system()
{
  dof_handler.distribute_dofs(fe);
  std::cout << "  Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;
  
  // these are some really cool objects.  It would be nice
  // to spend some time to figure out how they work.
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler,dsp);
  sparsity_pattern.copy_from(dsp);
  
  system_matrix.reinit(sparsity_pattern);
  
  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}

// assemble system
void Lec26a::assemble_system()
{
  QGauss<2> quadrature_formula(fe.degree + 1); // generally this should work
  FEValues<2> fe_values(fe,
                        quadrature_formula,
                        update_values | update_gradients | update_JxW_values);
  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  std::cout << "  Number of dofs per cell: " << dofs_per_cell
            << std::endl;
  
  // create containers for element contributions to the system matrix/ rhs
  FullMatrix<double> cell_matrix(dofs_per_cell,dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  
  // std::vector must be included in one of the deal.ii headers...
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  
  // iterate over all of the finite elements
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      cell_matrix = 0; // nice syntax
      cell_rhs = 0; // ditto
      
      // iterate over all quadrature points
      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
          for (const unsigned int i : fe_values.dof_indices())
            for(const unsigned int j : fe_values.dof_indices())
              cell_matrix(i,j) += 
                (fe_values.shape_grad(i,q_index) * // grad phi_i (x_q)
                 fe_values.shape_grad(j,q_index) * // grad phi_j (x_q)
                 fe_values.JxW(q_index));          // dx
                 
             // I do not have a right hand side for this problem
             // but for illustration purposes, I'll just go ahead and do this
             for (const unsigned int i : fe_values.dof_indices())
               cell_rhs(i) += (fe_values.shape_value(i,q_index) *
                               1.0 * // change this later...
                               fe_values.JxW(q_index));
        
        }
      // do this so I will know where to put the cell matrix values
      // into the system matrices
      cell->get_dof_indices(local_dof_indices);
      
      for (const unsigned int i : fe_values.dof_indices())
        for (const unsigned int j : fe_values.dof_indices())
          system_matrix.add(local_dof_indices[i],
                            local_dof_indices[j],
                            cell_matrix(i,j));
      
      for (const unsigned int i : fe_values.dof_indices())
        system_rhs(local_dof_indices[i]) += cell_rhs(i); 
    
    }
    
    // now I fiddle with the system matrix and rhs to 
    // handle boundary conditions
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           Functions::ZeroFunction<2>(),
                                           boundary_values);
  MatrixTools::apply_boundary_values(boundary_values,
                                     system_matrix,
                                     solution,
                                     system_rhs);

}

void Lec26a::solve()
{
  std::cout << "  Solving system. " << std::endl;
  
  SolverControl		solver_control(1000,1e-6*system_rhs.l2_norm());
  SolverCG<Vector<double>> solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
}

void Lec26a::output_results() const
{
  
  std::cout << "  Writing output data. " << std::endl;
  
  DataOut<2> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution,"solution");
  data_out.build_patches();
  
  std::ofstream output("solution.vtk");
  data_out.write_vtk(output);

}

// run the simulation
void Lec26a::run()
{
  make_grid();
  setup_system();
  assemble_system();
  solve();
  output_results();  
}

int main()
{
  Lec26a example_problem;
  example_problem.run();

  return 0;
}
