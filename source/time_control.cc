#include "time_control.h"

TimeControl::TimeControl()
	: is_steady(false)
	, initial_time(0.)
	, final_time(1.)
	, t(0.)
	, dt(0.01)
	, time_step(0)
{}

void TimeControl::advance_time()
{
	this->t += this->dt;
	this->time_step++;
}

void TimeControl::set_dt(const double new_dt)
{
	this->dt = new_dt;
}

double TimeControl::get_dt()
{
	return this->dt;
}

double TimeControl::get_current_time()
{
	return this->t;
}

double TimeControl::get_initial_time()
{
	return this->initial_time;
}

double TimeControl::get_final_time()
{
	return this->final_time;
}

unsigned int TimeControl::get_current_time_step()
{
	return this->time_step;
}

void TimeControl::declare_parameters(ParameterHandler &prm)
{
    prm.enter_subsection("Time parameters");
    prm.declare_entry("Initial time",
                      "0.",
                      Patterns::Double(0.),
                      "Initial time for the simulation.");
    prm.declare_entry("Final time",
                      "1.",
                      Patterns::Double(0.),
                      "Final time for the simulation.");
    prm.declare_entry("Delta t",
                      "0.01",
                      Patterns::Double(0.),
                      "Size of the time step.");
    prm.declare_entry("Is steady",
		      "false",
		      Patterns::Bool(),
                      "Whether the simulation is steady or not.");
    prm.leave_subsection();
}

void TimeControl::read_parameters(ParameterHandler &prm)
{
    prm.enter_subsection("Time parameters");
    initial_time = prm.get_double("Initial time");
    t = initial_time;
    final_time = prm.get_double("Final time");
    dt = prm.get_double("Delta t");
    is_steady = prm.get_bool("Is steady");
    prm.leave_subsection();
}
