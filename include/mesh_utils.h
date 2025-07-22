#ifndef mesh_utils_h
#define mesh_utils_h

#include "deal.II/base/exceptions.h"
#include "deal.II/base/function.h"
#include "deal.II/base/types.h"
#include "deal.II/distributed/tria.h"
#include "deal.II/grid/grid_generator.h"
#include <memory>

#include "boundary_value.h"

using namespace dealii;

template <int dim>
void makeCavity(parallel::distributed::Triangulation<dim> &tria, 
	    	std::map<types::boundary_id, std::shared_ptr<Function<dim>> > &boundary_conditions,
		const unsigned int n_refinements = 6)
{
    Assert((dim == 2), ExcNotImplemented());
    GridGenerator::hyper_cube(tria, 0., 1., true);
    tria.refine_global(n_refinements);

    boundary_conditions.insert({0, std::make_shared<Functions::ZeroFunction<dim>>(dim+1)});
    boundary_conditions.insert({1, std::make_shared<Functions::ZeroFunction<dim>>(dim+1)});
    boundary_conditions.insert({2, std::make_shared<Functions::ZeroFunction<dim>>(dim+1)});
    boundary_conditions.insert({3, std::make_shared<CavityInlet<dim>>()});

}

template <int dim>
void makeChannel(parallel::distributed::Triangulation<dim> &tria, 
	    	 std::map<types::boundary_id, std::shared_ptr<Function<dim>> > &boundary_conditions,
		 const unsigned int n_refinements = 6)
{
    Assert((dim == 2), ExcNotImplemented());
    Point<dim> a(0,0);
    Point<dim> b(2,0.5);
    GridGenerator::hyper_rectangle(tria, a, b, true);
    tria.refine_global(n_refinements);

    boundary_conditions.insert({0, std::make_shared<ChannelInlet<dim>>(0.5)});
    boundary_conditions.insert({2, std::make_shared<Functions::ZeroFunction<dim>>(dim+1)});
    boundary_conditions.insert({3, std::make_shared<Functions::ZeroFunction<dim>>(dim+1)});
}

template <int dim>
void makeCylinder(parallel::distributed::Triangulation<dim> &tria, 
		  std::map<types::boundary_id, std::shared_ptr<Function<dim>> > &boundary_conditions,
		  const unsigned int n_refinements = 2)
{
    Assert((dim == 2), ExcNotImplemented());
    GridGenerator::channel_with_cylinder(tria, 0.03, 2, 2.0, true);
    tria.refine_global(n_refinements);

    boundary_conditions.insert({0, std::make_shared<ChannelInlet<dim>>(0.41)});
    boundary_conditions.insert({2, std::make_shared<Functions::ZeroFunction<dim>>(dim+1)});
    boundary_conditions.insert({3, std::make_shared<Functions::ZeroFunction<dim>>(dim+1)});
    boundary_conditions.insert({4, std::make_shared<Functions::ZeroFunction<dim>>(dim+1)});
}



#endif
