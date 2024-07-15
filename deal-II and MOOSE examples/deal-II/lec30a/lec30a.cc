#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

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

template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  virtual double value(const Point<dim> & p,
                       const unsigned int component = 0) const override;
};

template <int dim>
double BoundaryValues<dim>::value(const Point<dim> & p,
                                  const unsigned int /*component*/) const
{
  double _x = p[0];
  double _y = p[1];
  
  return std::atan2(_y,_x);

}


class Lec30a
{
public:
  Lec30a ();
  void run ();
  
private:
  void make_grid ();
  void setup_system ();
  void assemble_system ();
  void solve ();
  void output_results () const;
  
  Triangulation<2> 	triangulation;
  FE_Q<2>		fe;
  DoFHandler<2>		dof_handler;
  SparsityPattern	sparsity_pattern;
  SparseMatrix<double>	system_matrix;
  Vector<double>	solution;
  Vector<double>	system_rhs;
  
  double 		_R;
  
};

// constructor
Lec30a::Lec30a()
    : fe(2), dof_handler(triangulation), _R(2)
 {}
 
 void 
 Lec30a::make_grid()
 {
   const Point<2> center(0,0);
   GridGenerator::hyper_ball(triangulation, 
                              center, _R);
   triangulation.refine_global(5);
   std::cout << "Number of active cells: " << triangulation.n_active_cells()
             << std::endl;
 }
 
 void
 Lec30a::setup_system()
 {
   dof_handler.distribute_dofs(fe);
   DynamicSparsityPattern dsp(dof_handler.n_dofs());
   DoFTools::make_sparsity_pattern(dof_handler,dsp);
   sparsity_pattern.copy_from(dsp);
   
   system_matrix.reinit(sparsity_pattern);
   
   solution.reinit(dof_handler.n_dofs());
   system_rhs.reinit(dof_handler.n_dofs());
 }
 
 void
 Lec30a::assemble_system()
 {
   QGauss<2> quadrature_formula(fe.degree + 1);
   FEValues<2> fe_values(fe,
                         quadrature_formula,
                         update_values | update_gradients | update_JxW_values);
   
   // create containers for element contributions to the system matrix and rhs
   const unsigned int 	dofs_per_cell = fe.n_dofs_per_cell();
   FullMatrix<double> 	cell_matrix(dofs_per_cell,dofs_per_cell);
   Vector<double> 	cell_rhs(dofs_per_cell);
   
   std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
   
   // iterate over all of the finite elements
   for (const auto &cell : dof_handler.active_cell_iterators())
   {
     fe_values.reinit(cell);
     cell_matrix = 0.;
     cell_rhs = 0.;
     
     // iterate over all quadrature points
     for (const unsigned int q_index : fe_values.quadrature_point_indices())
     {
       for (const unsigned int i : fe_values.dof_indices())
       {
         for(const unsigned int j : fe_values.dof_indices())
           cell_matrix(i,j) +=
               (fe_values.shape_grad(i,q_index)* //grad phi_i(x_q)
                fe_values.shape_grad(j,q_index)* //grad phi_j(x_q)
                fe_values.JxW(q_index));
       
      
       // I do not have a right hand side for this problem
       // but for illustration purposes, I'll just go ahead and do this
       
         cell_rhs(i) += (fe_values.shape_value(i,q_index) *
                         0.0 * // can now delete, I guess...
                         fe_values.JxW(q_index));    
                         
       }     
     
     } // q_index
     
     // put the cell matrix and cell rhs into the system data structures
     cell->get_dof_indices(local_dof_indices);
     for (const unsigned int i : fe_values.dof_indices())
     {
       for (const unsigned int j : fe_values.dof_indices())
         system_matrix.add(local_dof_indices[i],
                           local_dof_indices[j],
                           cell_matrix(i,j));
       
       system_rhs(local_dof_indices[i]) += cell_rhs(i);
       	  
     }  
     
     
   } // cell
   
   // now I fiddle with the system matrix and rhs to 
   // handle boundary conditions
   std::map<types::global_dof_index,double> boundary_values;
   VectorTools::interpolate_boundary_values(dof_handler,
                                            0, // all boundaries
                                            BoundaryValues<2>(),
                                            boundary_values);
                                            
   MatrixTools::apply_boundary_values(boundary_values,
                                      system_matrix,
                                      solution,
                                      system_rhs);
   
 } // assemble_system
 
 void
 Lec30a::solve()
 {
   std::cout << "  Solving system. " << std::endl;
   SolverControl		solver_control(1000,1e-6*system_rhs.l2_norm());
   SolverCG<Vector<double>>	solver(solver_control);
   solver.solve(system_matrix,solution,system_rhs,PreconditionIdentity()); 
 }
 
 void
 Lec30a::output_results() const
 {
   std::cout << "  Writing output data. " << std::endl;
   
   DataOut<2> data_out;
   data_out.attach_dof_handler(dof_handler);
   data_out.add_data_vector(solution,"solution");
   data_out.build_patches();
   
   std::ofstream output("solution.vtu");
   data_out.write_vtu(output); 
 }
 
 void
 Lec30a::run()
 {
   make_grid();
   setup_system();
   assemble_system();
   solve();
   output_results(); 
 }
 
 int main()
 {
   Lec30a example_problem;
   example_problem.run();
   
   return 0;
 }

 
 
 
 
 
 
 
 
 

