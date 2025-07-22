#ifndef boundary_value_h
#define boundary_value_h

#include "deal.II/base/point.h"
#include "deal.II/base/function.h"
#include <cmath>
#include <cstdlib>
using namespace dealii;

template <int dim>
class ChannelInlet: public Function<dim>
{
public: 
    ChannelInlet(const double channel_height) : Function<dim>(dim+1), h{channel_height} {}

    virtual double value(const Point<dim> &p, const unsigned int component) const override;

private:
    const double h;
};

template <int dim>
double ChannelInlet<dim>::value(const Point<dim> &p, const unsigned int component) const
{
    if (component == 0 && std::abs(p[0]) < 1e-10)
	return 4*p[1] * (this->h - p[1])/(this->h*this->h);
	// return 1.;
    return 0.;
}

template <int dim>
class ChannelPressureBoundary: public Function<dim>
{
public: 
    ChannelPressureBoundary() : Function<dim>(dim+1) {}

    virtual double value(const Point<dim> &p, const unsigned int component) const override;
};

template <int dim>
double ChannelPressureBoundary<dim>::value(const Point<dim> &p, const unsigned int component) const
{
    if (component == dim && std::abs(p[0]) < 1e-10)
	return 8.;
    return 0.;
}

template <int dim>
class CavityInlet: public Function<dim>
{
public: 
    CavityInlet() : Function<dim>(dim+1) {}

    virtual double value(const Point<dim> &p, const unsigned int component) const override;
};

template <int dim>
double CavityInlet<dim>::value(const Point<dim> &p, const unsigned int component) const
{
    if (component == 0 && std::abs(p[1]-1) < 1e-10)
	return 1;
    return 0.;
}

template <int dim>
class ExactSolution : public Function<dim>
{
public:
    ExactSolution() : Function<dim>(dim+1) {}

    virtual double value(const Point<dim> &p, const unsigned int component) const override;
};

template <int dim>
double ExactSolution<dim>::value(const Point<dim> &p, const unsigned int component) const
{
    if (component == 0)
    {
	return std::sin(p[0]) * std::sin(p[1] + this->get_time());
    } 
    else if (component == 1)
    {
	return std::cos(p[0]) * std::cos(p[1] + this->get_time());
    }
    else if (component == dim)
    {
	return std::cos(p[0]) * std::sin(p[1] + this->get_time());
    }
    else
    {
	return 0;
    }
}

#endif
