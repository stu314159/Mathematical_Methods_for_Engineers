// Deal.II libraries that I will (may) use
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
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
  //  void setup_system ();
  //  void assemble_system ();
 //   void solve ();
  //  void output_results () const;
    
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
  : fe(1), dof_handler(triangulation), 
    lower_left(0,0), upper_right(5,3)
{}

// make_grid member function
void Lec26a::make_grid()
{
  GridGenerator::hyper_rectangle(triangulation, lower_left, upper_right);
  triangulation.refine_global(4);
  
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

// run the simulation
void Lec26a::run()
{
  make_grid();
  
}

int main()
{
  Lec26a example_problem;
  example_problem.run();

  return 0;
}
