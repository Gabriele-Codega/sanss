#ifndef navier_stokes_h
#define navier_stokes_h


#include "deal.II/base/data_out_base.h"
#include "deal.II/base/mpi.h"
#include "deal.II/base/index_set.h"
#include "deal.II/base/quadrature_lib.h"
#include "deal.II/base/conditional_ostream.h"
#include "deal.II/base/types.h"
#include "deal.II/base/parameter_handler.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/distributed/tria.h"
#include "deal.II/grid/tria.h"
#include "deal.II/fe/fe_system.h"
#include "deal.II/fe/fe_q.h"
#include "deal.II/lac/affine_constraints.h"
#include "deal.II/lac/sparsity_pattern.h"
#include "deal.II/lac/sparse_matrix.h"
#include "deal.II/lac/trilinos_vector.h"
#include "deal.II/lac/trilinos_sparse_matrix.h"

#include "time_control.h"

#include <memory>
#include <string>


using namespace dealii;


enum class LinearSolverType
{
    direct,
    gmres
};

enum class StabilisationType
{
    none,
    gls
};

enum class TestCase
{
    cavity,
    cylinder,
    channel
};

template <int dim>
class NavierStokes
{
public:
	NavierStokes(const std::string);
	void solve();
private:
	ParameterHandler prm;
	void declare_parameters();
	void read_parameters(const std::string fname);

	unsigned int	    degree_u;
	unsigned int	    degree_p;
	std::shared_ptr<FESystem<dim>>  fe;
	DoFHandler<dim>	    dof_handler;
	parallel::distributed::Triangulation<dim>	triangulation;
	unsigned int	degree_q;
	QGauss<dim>	quadrature_cells;
	std::map<types::boundary_id, std::shared_ptr<Function<dim>>>	boundary_conditions;
	AffineConstraints<double>	zero_constraints;
	AffineConstraints<double>	nonzero_constraints;
	bool	pressure_has_zero_mean;
	void	make_constraints();
	void	pressure_mean_to_zero();

	TestCase    test_case;
	unsigned int n_refinements;

	TimeControl	time_control;
	double			nu;
	StabilisationType	stabilisation;


	void	setup_system();
	void	assemble(const bool, const bool);
	void	assemble_system(const bool);
	void	assemble_rhs(const bool);

	void 	solve_newton(bool);
	void	solve_linear_system();

	unsigned int max_nonlinear_iterations;
	unsigned int max_linear_iterations;
	double nonlinear_tolerance;
	double linear_tolerance;
	LinearSolverType linear_method;

	TrilinosWrappers::MPI::Vector	solution_nm1;
	TrilinosWrappers::MPI::Vector	solution;
	TrilinosWrappers::MPI::Vector	newton_update;

	SparsityPattern			sparsity_pattern;
	TrilinosWrappers::SparseMatrix	system_matrix;
	TrilinosWrappers::MPI::Vector	system_rhs;

	MPI_Comm	mpi_communicator;
	int		mpi_rank;

	IndexSet	locally_owned_dofs;
	IndexSet	locally_relevant_dofs;

	ConditionalOStream	pcout;
	std::string output_directory;
	unsigned int write_interval;
	bool   compute_strain_rate;
	bool   compute_vorticity;

	void	output_results(TimeControl &);

};

#endif
