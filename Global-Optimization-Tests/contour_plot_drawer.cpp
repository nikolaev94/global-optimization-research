
#include "contour_plot_drawer.h"


void ContourPlotDataDumper::dump(const std::string& filename)
{
	std::ofstream ofstream(filename);

	OptProblem::ContourData data(DEFAULT_GRID_SIZE);

	target_problem->initContourData(data);

	for (auto i = 0; i < data.grid_size; i++)
	{
		for (auto j = 0; j < data.grid_size; j++)
		{
			ofstream << data.x_values[i] << ' ' << data.y_values[j] << ' '
				<< data.z_values[i][j] << std::endl;
		}
	}

	ofstream.close();
}

