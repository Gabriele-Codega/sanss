#ifndef time_control_h
#define time_control_h

#include "deal.II/base/parameter_handler.h"
#include "deal.II/base/patterns.h"
#include <string>

using namespace dealii;

class TimeControl
{
public:
	TimeControl();
	void advance_time();
	void set_dt(const double);
	double get_dt();
	double get_current_time();
	double get_initial_time();
	double get_final_time();
	unsigned int get_current_time_step();

	bool is_steady;

	void declare_parameters(ParameterHandler &);
	void read_parameters(ParameterHandler &);
private:
	double initial_time;
	double final_time;
	double t;
	double dt;
	unsigned int time_step;
};


#endif
