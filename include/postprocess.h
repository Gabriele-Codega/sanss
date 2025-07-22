#ifndef postprocess_h
#define postprocess_h

#include "deal.II/base/exceptions.h"
#include "deal.II/lac/vector.h"
#include "deal.II/numerics/data_postprocessor.h"

using namespace dealii;

template <int dim>
class VorticityPostprocessor : public DataPostprocessorVector<dim>
{
public:
    VorticityPostprocessor() : DataPostprocessorVector<dim>("vorticity", update_gradients){}

    virtual
    void
    evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
			  std::vector<Vector<double>> &computed_quantities) const override
    {
	AssertDimension(input_data.solution_gradients.size(), computed_quantities.size());

	for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
	{
	    AssertDimension(computed_quantities[p].size(), dim);
	    if constexpr (dim == 3)
	    {
		computed_quantities[p](0) = input_data.solution_gradients[p][2][1] - input_data.solution_gradients[p][1][2];
		computed_quantities[p](1) = input_data.solution_gradients[p][0][2] - input_data.solution_gradients[p][2][0];
		computed_quantities[p](2) = input_data.solution_gradients[p][1][0] - input_data.solution_gradients[p][0][1];
	    }
	    if constexpr (dim == 2)
	    {
		computed_quantities[p](0) = input_data.solution_gradients[p][1][0] - input_data.solution_gradients[p][0][1];
	    }
	}

    }
};

template <int dim>
class StrainRatePostprocessor : public DataPostprocessorTensor<dim>
{
public:
    StrainRatePostprocessor() : DataPostprocessorTensor<dim>("strain_rate", update_gradients){}

    virtual
    void
    evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
			  std::vector<Vector<double>> &computed_quantities) const override
    {
	AssertDimension(input_data.solution_gradients.size(), computed_quantities.size());

	for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
	{
	    AssertDimension (computed_quantities[p].size(),
                         (Tensor<2,dim>::n_independent_components));

	    for (unsigned int i=0; i<dim; ++i)
		for (unsigned int j=0; j<dim; ++j)
		    computed_quantities[p][Tensor<2,dim>::component_to_unrolled_index(TableIndices<2>(i,j))]
			= 0.5 * (input_data.solution_gradients[p][i][j] + input_data.solution_gradients[p][j][i]);
	}

    }
};

#endif
