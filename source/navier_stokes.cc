// sanss includes
#include "navier_stokes.h"
#include "mesh_utils.h"
#include "postprocess.h"

// deal.ii includes
#include "deal.II/base/mpi.h"
#include "deal.II/base/exceptions.h"
#include "deal.II/base/patterns.h"
#include "deal.II/base/data_out_base.h"
#include "deal.II/base/index_set.h"
#include "deal.II/base/tensor.h"
#include "deal.II/base/types.h"

#include "deal.II/fe/fe_q.h"
#include "deal.II/fe/fe_system.h"
#include "deal.II/fe/component_mask.h"
#include "deal.II/fe/fe_values_extractors.h"

#include "deal.II/dofs/dof_tools.h"
#include "deal.II/dofs/dof_renumbering.h"

#include "deal.II/lac/affine_constraints.h"
#include "deal.II/lac/dynamic_sparsity_pattern.h"
#include "deal.II/lac/sparsity_tools.h"
#include "deal.II/lac/full_matrix.h"
#include "deal.II/lac/vector_operation.h"
#include "deal.II/lac/solver_control.h"
#include "deal.II/lac/trilinos_solver.h"
#include "deal.II/lac/trilinos_precondition.h"

#include "deal.II/numerics/vector_tools.h"
#include "deal.II/numerics/data_component_interpretation.h"
#include "deal.II/numerics/data_out.h"
#include "deal.II/numerics/vector_tools_mean_value.h"

// c++ includes
#include <filesystem>
#include <iostream>
#include <memory>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <vector>

template <int dim>
NavierStokes<dim>::NavierStokes(const std::string fname)
    : degree_u(1)
    , degree_p(1)
    , triangulation(MPI_COMM_WORLD)
    , degree_q(degree_u + 2)
    , quadrature_cells(degree_u + 2)
    , max_nonlinear_iterations(10)
    , max_linear_iterations(100)
    , nonlinear_tolerance(1e-6)
    , linear_tolerance(1e-6)
    , mpi_communicator(MPI_COMM_WORLD)
    , mpi_rank(Utilities::MPI::this_mpi_process(mpi_communicator))
    , pcout(std::cout, (mpi_rank == 0))
    , output_directory("simulation_output")
    , compute_strain_rate(false)
    , compute_vorticity(false)
{
    declare_parameters();
    read_parameters(fname);

    fe = std::make_shared<FESystem<dim>>(FE_Q<dim>(degree_u), dim, FE_Q<dim>(degree_p), 1);

    if (test_case == TestCase::cavity)
	makeCavity(triangulation, boundary_conditions,6);
    else if (test_case == TestCase::cylinder)
	makeCylinder(triangulation, boundary_conditions, n_refinements);
    else if (test_case == TestCase::channel)
	makeChannel(triangulation, boundary_conditions);

    dof_handler.clear();
    dof_handler.reinit(triangulation);

}

template <int dim>
void NavierStokes<dim>::declare_parameters()
{
    time_control.declare_parameters(prm);

    prm.enter_subsection("Geometry");
    prm.declare_entry("Test case",
		      "cavity",
		      Patterns::Selection("cavity|cylinder|channel"),
		      "Predefined test case to run. Either `cavity` for the lid-driven cavity flow, `cylinder` for the 2D flow past a cylinder, or `channel` for a 2D plane channel flow.");
    prm.declare_entry("Number of refinements",
		      "1",
		      Patterns::Integer(1),
		      "How many times should the mesh be refined.");
    prm.leave_subsection();

    prm.enter_subsection("Governing equations");
    prm.declare_entry("Kinematic viscosity",
		      "0.01",
		      Patterns::Double(0.0),
		      "Kinematic viscosity.");
    prm.declare_entry("Stabilisation",
		      "gls",
		      Patterns::Selection("none|gls"),
		      "What kind of stabilisation should be used. Either `none` or `gls`.");
    prm.leave_subsection();

    prm.enter_subsection("Finite Element");
    prm.declare_entry("u degree",
		      "1",
		      Patterns::Integer(1),
		      "Polynomial degree for velocity.");
    prm.declare_entry("p degree",
		      "1",
		      Patterns::Integer(1),
		      "Polynomial degree for pressure.");
    prm.declare_entry("Q degree",
		      "3",
		      Patterns::Integer(1),
		      "Number of quadrature points. ");
    prm.declare_entry("Pressure has zero mean",
		      "true",
		      Patterns::Bool(),
		      "Whether the mean pressure should be zero.");
    prm.leave_subsection();

    prm.enter_subsection("Output control");
    prm.declare_entry("Output directory",
		      "simulation_output",
		      Patterns::DirectoryName(),
		      "Name of the directory where output should be placed.");
    prm.declare_entry("Write interval",
		      "1",
		      Patterns::Integer(1),
		      "Number of time steps after which output is generated.");
    prm.declare_entry("Compute strain rate",
		      "false",
		      Patterns::Bool(),
		      "Whether to compute the strain rate for output.");
    prm.declare_entry("Compute vorticity",
		      "false",
		      Patterns::Bool(),
		      "Whether to compute the vorticity for output.");
    prm.leave_subsection();

    prm.enter_subsection("Newton method");
    prm.declare_entry("Max nonlinear iterations",
		      "10",
		      Patterns::Integer(1),
		      "Maximum number of iterations in the nonlinear solver.");
    prm.declare_entry("Nonlinear tolerance",
		      "1e-6",
		      Patterns::Double(0),
		      "Tolerance of the nonlinear solver");
    prm.enter_subsection("Linear solver");
	prm.declare_entry("Method",
			  "gmres",
			  Patterns::Selection("direct|gmres"),
			  "Algorithm for the linear solver. Either `direct` or `gmres`.");
	prm.declare_entry("Max linear iterations",
			  "100",
			  Patterns::Integer(1),
			  "Maximum number of iterations in the linear solver.");
	prm.declare_entry("Linear tolerance",
			  "1e-6",
			  Patterns::Double(0),
			  "Tolerance of the linear solver");
    prm.leave_subsection(); // linear
    prm.leave_subsection(); // nonlinear
}

template <int dim>
void NavierStokes<dim>::read_parameters(const std::string fname)
{
    std::ifstream   file(fname);
    AssertThrow(file, ExcFileNotOpen(fname));

    prm.parse_input(file);

    time_control.read_parameters(prm);

    prm.enter_subsection("Geometry");
    auto val = prm.get("Test case");
    if (val == "cavity")
	test_case = TestCase::cavity;
    else if (val == "cylinder")
	test_case = TestCase::cylinder;
    else if (val == "channel")
	test_case = TestCase::channel;
    else
	std::runtime_error("Unkown test case. Please select between `cavity`, `cylinder` or `channel`.");
    n_refinements = prm.get_integer("Number of refinements");
    prm.leave_subsection();

    prm.enter_subsection("Governing equations");
    nu = prm.get_double("Kinematic viscosity");
    if (prm.get("Stabilisation") == "none")
	stabilisation = StabilisationType::none;
    else if (prm.get("Stabilisation") == "gls")
	stabilisation = StabilisationType::gls;
    else
	std::runtime_error("Unknown stabilisation type. Possible values are `none` or `gls`");
    prm.leave_subsection();

    prm.enter_subsection("Finite Element");
    degree_u = prm.get_integer("u degree");
    degree_p = prm.get_integer("p degree");
    degree_q = prm.get_integer("Q degree");
    pressure_has_zero_mean = prm.get_bool("Pressure has zero mean");
    prm.leave_subsection();

    prm.enter_subsection("Output control");
    output_directory = prm.get("Output directory");
    {
	std::filesystem::path dir(output_directory);
	if (!std::filesystem::exists(dir))
	    std::filesystem::create_directories(dir);
    }
    write_interval = prm.get_integer("Write interval");
    compute_strain_rate = prm.get_bool("Compute strain rate");
    compute_vorticity = prm.get_bool("Compute vorticity");
    prm.leave_subsection();

    prm.enter_subsection("Newton method");
    max_nonlinear_iterations = prm.get_integer("Max nonlinear iterations");
    nonlinear_tolerance = prm.get_double("Nonlinear tolerance");
    prm.enter_subsection("Linear solver");
	if (prm.get("Method") == "direct")
	    linear_method = LinearSolverType::direct;
	else if (prm.get("Method") == "gmres")
	    linear_method = LinearSolverType::gmres;
	else
	    std::runtime_error("Unknown linear solver type. Available options are `direct` and `gmres`.");
	max_linear_iterations = prm.get_integer("Max linear iterations");
	linear_tolerance = prm.get_double("Linear tolerance");
    prm.leave_subsection(); // linear
    prm.leave_subsection(); // nonlinear
}



template <int dim>
void NavierStokes<dim>::setup_system()
{
    dof_handler.distribute_dofs(*fe);
    DoFRenumbering::Cuthill_McKee(dof_handler);

    locally_owned_dofs	    = dof_handler.locally_owned_dofs();
    locally_relevant_dofs   = DoFTools::extract_locally_relevant_dofs(dof_handler);


    // current and previous solutions are ghosted since they are needed in assembly
    solution_nm1.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    // newton update and system rhs anre non-ghosted since not needed in assembly
    // (though they need to be written to)
    newton_update.reinit(locally_owned_dofs, mpi_communicator);
    system_rhs.reinit(locally_owned_dofs, mpi_communicator);

    make_constraints();

    DynamicSparsityPattern dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, zero_constraints);
    SparsityTools::distribute_sparsity_pattern(dsp, dof_handler.locally_owned_dofs(), mpi_communicator, locally_relevant_dofs);

    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(locally_owned_dofs, locally_owned_dofs, sparsity_pattern, mpi_communicator);

}


template <int dim>
void NavierStokes<dim>::make_constraints()
{
    const FEValuesExtractors::Vector	velocities(0);
    const FEValuesExtractors::Scalar	pressure(dim);

    zero_constraints.clear();
    nonzero_constraints.clear();

    zero_constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler,
					      zero_constraints);
    for (const auto &[id, type] : boundary_conditions)
    {
	DoFTools::make_zero_boundary_constraints(dof_handler,
						  id,
						  zero_constraints,
						  fe->component_mask(velocities));
    }

    nonzero_constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler,
					    nonzero_constraints);
    for (const auto &[id, function] : boundary_conditions)
    {
	VectorTools::interpolate_boundary_values(dof_handler,
						 id,
						 *function,
						 nonzero_constraints,
						 fe->component_mask(velocities));
    }
    zero_constraints.close();
    nonzero_constraints.close();
}


template <int dim>
void NavierStokes<dim>::pressure_mean_to_zero()
{
    TrilinosWrappers::MPI::Vector   non_ghosted_solution(locally_owned_dofs, mpi_communicator);
    non_ghosted_solution = solution;
    double mean_pressure = VectorTools::compute_mean_value(dof_handler, quadrature_cells, solution, dim);
    const FEValuesExtractors::Scalar pressure(dim);
    const ComponentMask	pressure_mask = fe->component_mask(pressure);
    const IndexSet pressure_dofs = DoFTools::extract_dofs(dof_handler, pressure_mask);
    for (const auto i : pressure_dofs)
    {
    non_ghosted_solution[i] -= mean_pressure;
    }
    solution = non_ghosted_solution;
}


template <int dim>
void NavierStokes<dim>::assemble(const bool initial_step,
				 const bool assemble_matrix)
{
    system_matrix = 0;
    system_rhs = 0;

    FEValues<dim> fe_values(*fe, quadrature_cells,
			    update_values | update_gradients | update_quadrature_points | update_hessians | update_JxW_values);
    const unsigned int	n_dofs = fe->n_dofs_per_cell();
    const unsigned int  n_q_points = fe_values.n_quadrature_points;

    FullMatrix<double>	local_matrix(n_dofs, n_dofs);
    Vector<double>	local_rhs(n_dofs);

    std::vector<types::global_dof_index>    local_dof_indices(n_dofs);

    double  dt = time_control.get_dt();
    double  dt_inv = 1./dt;

    // current solution and derivatives
    std::vector<Tensor<1, dim>>	current_velocity(n_q_points);
    std::vector<Tensor<2, dim>>	current_velocity_grad(n_q_points);
    std::vector<double>		current_velocity_div(n_q_points);
    std::vector<Tensor<1, dim>>	current_velocity_lapl(n_q_points);
    std::vector<double>		current_pressure(n_q_points);
    std::vector<Tensor<1, dim>>	current_pressure_grad(n_q_points);

    // velocity at the previous time
    std::vector<Tensor<1, dim>>	velocity_nm1_vec(n_q_points);

    // residual and jacobian of the strong form
    std::vector<Tensor<1, dim>>	strong_residual(n_q_points);
    std::vector<std::vector<Tensor<1, dim>>>	strong_jacobian(n_q_points, std::vector<Tensor<1, dim>>(n_dofs));

    // shape functions
    std::vector<std::vector<Tensor<1, dim>>>	phi_u(n_q_points, std::vector<Tensor<1, dim>>(n_dofs));
    std::vector<std::vector<Tensor<2, dim>>>	grad_phi_u(n_q_points, std::vector<Tensor<2, dim>>(n_dofs));
    std::vector<std::vector<double>>		div_phi_u(n_q_points, std::vector<double>(n_dofs));
    std::vector<std::vector<Tensor<3, dim>>>	hess_phi_u(n_q_points, std::vector<Tensor<3, dim>>(n_dofs));
    std::vector<std::vector<Tensor<1, dim>>>	lapl_phi_u(n_q_points, std::vector<Tensor<1, dim>>(n_dofs));
    std::vector<std::vector<double>>		phi_p(n_q_points, std::vector<double>(n_dofs));
    std::vector<std::vector<Tensor<1, dim>>>	grad_phi_p(n_q_points, std::vector<Tensor<1, dim>>(n_dofs));

    const FEValuesExtractors::Vector	velocities(0);
    const FEValuesExtractors::Scalar	pressure(dim);
    
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
	if (!cell->is_locally_owned())
	    continue;

	fe_values.reinit(cell);
	if (assemble_matrix)
	    local_matrix = 0;
	local_rhs = 0;


	// get solution values
	fe_values[velocities].get_function_values(solution, current_velocity);
	fe_values[velocities].get_function_gradients(solution, current_velocity_grad);
	fe_values[velocities].get_function_divergences(solution, current_velocity_div);
	fe_values[velocities].get_function_laplacians(solution, current_velocity_lapl);
	
	fe_values[pressure].get_function_values(solution, current_pressure);
	fe_values[pressure].get_function_gradients(solution, current_pressure_grad);

	fe_values[velocities].get_function_values(solution_nm1, velocity_nm1_vec);

	// strong residual at current solution
	for (const unsigned int q : fe_values.quadrature_point_indices())
	{
	    strong_residual[q] = current_velocity_grad[q] * current_velocity[q] 
				+ current_pressure_grad[q] 
				- nu * current_velocity_lapl[q];
	    if (!time_control.is_steady)
		strong_residual[q] += dt_inv * (current_velocity[q] - velocity_nm1_vec[q]);
	}

	// get shape values
	for (const unsigned int q : fe_values.quadrature_point_indices())
	    for (unsigned int k = 0; k < n_dofs; ++k)
	    {
		phi_u[q][k] = fe_values[velocities].value(k,q);
		grad_phi_u[q][k] = fe_values[velocities].gradient(k,q);
		div_phi_u[q][k] = fe_values[velocities].divergence(k,q);
		hess_phi_u[q][k] = fe_values[velocities].hessian(k,q);
		for (unsigned int d = 0; d < dim; ++d)
		    lapl_phi_u[q][k][d] = trace(hess_phi_u[q][k][d]);
		phi_p[q][k] = fe_values[pressure].value(k,q);
		grad_phi_p[q][k] = fe_values[pressure].gradient(k,q);
	    }

	// precompute strong jacobian at current solution in shape directions
	for (const unsigned int q : fe_values.quadrature_point_indices())
	    for (unsigned int k = 0; k < n_dofs; ++k)
	    {
		strong_jacobian[q][k] = - nu * lapl_phi_u[q][k] 
					+ (grad_phi_u[q][k] * current_velocity[q]) 
					+ (current_velocity_grad[q] * phi_u[q][k]) 
					+ grad_phi_p[q][k];
		if (!time_control.is_steady)
		    strong_jacobian[q][k] += dt_inv * phi_u[q][k] ;
	    }

	// define characteristic size h as the diameter of
	// a sphere with the same volume as current cell.
	double h = 0;
	double cell_volume = 0;
	for (const unsigned int q : fe_values.quadrature_point_indices())
	    cell_volume += fe_values.JxW(q);
	if (dim == 2)
	    h = std::sqrt(4 * cell_volume/numbers::PI)/degree_u ;
	else if (dim == 3)
	    h = std::sqrt(6 * cell_volume/numbers::PI)/degree_u ;
	else
	    std::runtime_error("Error in assembly: can only perform simulations in 2D or 3D.");

	for (const unsigned int q : fe_values.quadrature_point_indices())
	{
	    // GaLS and LSIC constants
	    double tau	    = 0;
	    double tau_lsic = 0;

	    double u_mag = current_velocity[q].norm();

	    tau = time_control.is_steady ? 1./std::sqrt(4*u_mag*u_mag/(h*h) + 9 * 16 * nu*nu/(h*h*h*h)) : 1./std::sqrt(dt_inv*dt_inv +  4*u_mag*u_mag/(h*h) + 9 * 16 * nu*nu/(h*h*h*h));
	    tau_lsic = u_mag * h * 0.5;

	    for (unsigned int i = 0; i < n_dofs; ++i)
	    {
		if (assemble_matrix)
		{
		    for (unsigned int j = 0; j < n_dofs; ++j)
		    {

			if (!time_control.is_steady)
			    local_matrix(i,j) += ( dt_inv * phi_u[q][i] * phi_u[q][j] ) * fe_values.JxW(q);

			local_matrix(i,j) += ( nu * scalar_product(grad_phi_u[q][i], grad_phi_u[q][j]) 
					    + phi_u[q][i] * (grad_phi_u[q][j] * current_velocity[q])
					    + phi_u[q][i] * (current_velocity_grad[q] * phi_u[q][j]) 
					    - div_phi_u[q][i] * phi_p[q][j] 
					    + phi_p[q][i] * div_phi_u[q][j]) * fe_values.JxW(q);

			if (stabilisation == StabilisationType::gls)
			{
			    local_matrix(i,j) += tau * (strong_residual[q] * (grad_phi_u[q][i] * phi_u[q][j]) + strong_jacobian[q][j] * (grad_phi_u[q][i] * current_velocity[q] + grad_phi_p[q][i] - nu * lapl_phi_u[q][i])) * fe_values.JxW(q);

			    local_matrix(i,j) += tau_lsic * (div_phi_u[q][i] * div_phi_u[q][j]) * fe_values.JxW(q);
			}

		    }
		}
		
		if (!time_control.is_steady)
		    local_rhs(i) +=	- (  dt_inv * (current_velocity[q] - velocity_nm1_vec[q]) * phi_u[q][i]) * fe_values.JxW(q);

		local_rhs(i) += - (nu * scalar_product(grad_phi_u[q][i], current_velocity_grad[q]) 
				    + phi_u[q][i] * (current_velocity_grad[q] * current_velocity[q]) 
				    - div_phi_u[q][i] * current_pressure[q] 
				    + phi_p[q][i] * current_velocity_div[q]) * fe_values.JxW(q);

		if (stabilisation == StabilisationType::gls)
		{
		    local_rhs(i) += -tau * (strong_residual[q] * ( grad_phi_u[q][i] * current_velocity[q] + grad_phi_p[q][i] - nu * lapl_phi_u[q][i])) * fe_values.JxW(q);

		    local_rhs(i) += -tau_lsic * (current_velocity_div[q] * div_phi_u[q][i]) * fe_values.JxW(q);
		}
	    }

	}
	cell->get_dof_indices(local_dof_indices);

	const AffineConstraints<double> &constraints_used = initial_step ? nonzero_constraints : zero_constraints;

	if (assemble_matrix)
	    constraints_used.distribute_local_to_global(local_matrix, 
							local_rhs, 
							local_dof_indices, 
							system_matrix, 
							system_rhs);
	else
	    constraints_used.distribute_local_to_global(local_rhs, 
							local_dof_indices,
							system_rhs);

    }
    if (assemble_matrix)
	system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
}

template <int dim>
void NavierStokes<dim>::assemble_system(const bool initial_step)
{
    assemble(initial_step, true);
}

template <int dim>
void NavierStokes<dim>::assemble_rhs(const bool initial_step)
{
    assemble(initial_step, false);
}

template <int dim>
void NavierStokes<dim>::solve_newton(bool initial_step)
{
    double  update_magnitude = nonlinear_tolerance + 1;
    double  current_residual = nonlinear_tolerance + 1;
    double  old_residual     = nonlinear_tolerance + 1;
    unsigned int it = 0;

    // auxiliary vectors for applying boundary conditions and line search
    TrilinosWrappers::MPI::Vector   non_ghosted_solution(locally_owned_dofs, mpi_communicator);
    TrilinosWrappers::MPI::Vector   backup_solution(locally_owned_dofs, mpi_communicator);

    // first assembly and residual for the initial guess
    assemble_rhs(initial_step);
    nonzero_constraints.set_zero(system_rhs);
    current_residual = system_rhs.l2_norm();
    old_residual = current_residual;
    pcout << "\n  Newton solver -- initial residual = " << current_residual << std::endl;


    while ((current_residual > nonlinear_tolerance && ++it <= max_nonlinear_iterations) || it == 0)
    {
	pcout << "    Iteration " << it << std::endl;
	
	pcout << "      Assemblying..."; 
	std::flush(std::cout);
	assemble_system(initial_step);
	pcout << " done." << std::endl;

	
	const AffineConstraints<double> &constraints_used = initial_step ? nonzero_constraints : zero_constraints;

	// setup and solve the linear system
	constraints_used.set_zero(newton_update);
	solve_linear_system();
	constraints_used.distribute(newton_update);

	update_magnitude = newton_update.l2_norm();
	backup_solution = solution;

	unsigned int alpha_iter = 0;
	double latest_alpha_residual = current_residual;
	for (double alpha = 1.; alpha > 1e-3; alpha*=0.5)
	{
	    non_ghosted_solution = backup_solution;
	    non_ghosted_solution.add(alpha, newton_update);
	    solution = non_ghosted_solution;

	    if (pressure_has_zero_mean)
		pressure_mean_to_zero();

	    assemble_rhs(initial_step);
	    nonzero_constraints.set_zero(system_rhs);
	    current_residual = system_rhs.l2_norm();

	    pcout << "      Line search iteration: " << alpha_iter 
		  << "  --  residual = " << current_residual << std::endl;

	    if (initial_step)
		break;

	    if (current_residual > latest_alpha_residual && alpha_iter != 0)
	    {
		alpha *= 2;
		non_ghosted_solution = backup_solution;
		non_ghosted_solution.add(alpha, newton_update);
		solution = non_ghosted_solution;

		if (pressure_has_zero_mean)
		    pressure_mean_to_zero();

		current_residual = latest_alpha_residual;
		break;
	    }

	    if (current_residual < old_residual)
		break;

	    latest_alpha_residual = current_residual;
	    alpha_iter++;
	} // line search

	initial_step = false;

	old_residual = current_residual;
	pcout << "\n      increment norm = " << update_magnitude
	      << " -- residual = " << current_residual
	      << std::endl;
    }

}

template <int dim>
void NavierStokes<dim>::solve_linear_system()
{
    SolverControl   solver_control(max_linear_iterations, linear_tolerance*system_rhs.l2_norm());
    if (linear_method == LinearSolverType::gmres)
    {
	TrilinosWrappers::SolverGMRES	    solver(solver_control);
	TrilinosWrappers::PreconditionAMG	    preconditioner;
	std::vector<std::vector<bool>>	    constant_modes;
	ComponentMask			    components(dim+1, true);
	DoFTools::extract_constant_modes(dof_handler, components, constant_modes);
	TrilinosWrappers::PreconditionAMG::AdditionalData   additional_data(false,
									    true,
									    1,
									    false,
									    1e-4,
									    constant_modes,
									    2,
									    0,
									    false,
									    "ILU",
									    "ILU");
	preconditioner.initialize(system_matrix, additional_data);
	solver.solve(system_matrix, newton_update, system_rhs, preconditioner);

	pcout << "      Number of GMRES iterations = " << solver_control.last_step() << " -- final residual = " << solver_control.last_value()
	  << std::endl;
    }
    else if (linear_method == LinearSolverType::direct)
    {
	TrilinosWrappers::SolverDirect::AdditionalData additional_data(false, "Amesos_Mumps");
	TrilinosWrappers::SolverDirect solver(solver_control, additional_data);
	solver.initialize(system_matrix); 
	solver.solve(newton_update, system_rhs);
	pcout << "      Final linear residual = " << solver_control.last_value() << std::endl;
    }
}


template <int dim>
void NavierStokes<dim>::solve()
{
    setup_system();

    pcout << "Running on " << Utilities::MPI::n_mpi_processes(mpi_communicator) << " MPI processes." << std::endl;
    pcout << "Using Q" << degree_u << "-Q" << degree_p << " finite elements -- " << dof_handler.n_dofs() << " DoFs." << std::endl;

    // apply boundary constraints before simulation starts
    {
	TrilinosWrappers::MPI::Vector   non_ghosted_solution(locally_owned_dofs, mpi_communicator);
	nonzero_constraints.distribute(non_ghosted_solution);
	solution_nm1 = non_ghosted_solution;

	// write initial condition
	output_results(time_control);
    }

    if (!time_control.is_steady)
    {
	pcout << "Running unsteady simulation:" << std::endl
	      << "  Initial time = " << time_control.get_initial_time() << ", final time = " << time_control.get_final_time() << ", dt = " << time_control.get_dt() << std::endl;
	while (time_control.get_current_time() < time_control.get_final_time())
	{
	    time_control.advance_time();

	    bool initial_step = (time_control.get_current_time_step() == 1);

	    pcout << "\n+"+std::string(30,'-')+"+" <<std::endl
		  << "|" << std::setw(15) << " Time = " << std::setw(10) << time_control.get_current_time() << std::setw(6) << "|" << std::endl
		  << "+"+std::string(30,'-')+"+" <<std::endl;

	    solve_newton(initial_step);
	    solution_nm1 = solution;

	    if (time_control.get_current_time_step() % write_interval == 0)
	    {
		pcout << "  Writing to file..." << std::endl;
		output_results(time_control);
	    }
	}
    }
    else
    {
	pcout << " Running steady simulation" << std::endl;
	bool initial_step = true;
	solve_newton(initial_step);
	pcout << "  Writing to file..." << std::endl;
	output_results(time_control);
    }
}

template <int dim>
void NavierStokes<dim>::output_results(TimeControl &time_control)
{
    unsigned int time_step = time_control.get_current_time_step();
    double t = time_control.get_current_time();

    std::string	file_name = "solution_"+Utilities::int_to_string(time_step/write_interval,6)+".vtu";

    std::vector<std::string> solution_names(dim, "velocity");
    solution_names.push_back("pressure");
 
    std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back (DataComponentInterpretation::component_is_scalar);

    VorticityPostprocessor<dim>	    vorticity_pp;
    StrainRatePostprocessor<dim>    strain_rate_pp; 
    DataOut<dim>    data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.attach_triangulation(triangulation);

    data_out.add_data_vector(solution, solution_names, DataOut<dim>::type_dof_data, data_component_interpretation);

    if (compute_vorticity)
    {
	data_out.add_data_vector(solution, vorticity_pp);
    }

    if (compute_strain_rate)
    {
	data_out.add_data_vector(solution, strain_rate_pp);
    }

    DataOutBase::VtkFlags flags;
    flags.time			    = t;
    flags.cycle			    = time_step;
    flags.write_higher_order_cells  = true;
    data_out.set_flags(flags);

    unsigned int n_patches = std::min((unsigned int)4, degree_u);
    data_out.build_patches(n_patches);

    std::filesystem::path full_path(output_directory);
    full_path = full_path / file_name;
    data_out.write_vtu_in_parallel(full_path, mpi_communicator);
}



template class NavierStokes<2>;
template class NavierStokes<3>;
