#include "deal.II/base/conditional_ostream.h"
#include "deal.II/base/mpi.h"
#include "navier_stokes.h"
#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv)
{

    try
    {
        Utilities::MPI::MPI_InitFinalize    mpi_initialisation(argc,argv,1);
        ConditionalOStream pcout(std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0));

        if (argc < 3)
        {
            pcout << "Usage: " << argv[0] << " <space_dim>" << " parameter_file" << std::endl;
            return 1;
        }

        pcout << "                  _   _    _____        " << std::endl
              << "                 | \\ | |  / ____|       " << std::endl
              << "    ___    __ _  |  \\| | | (___    ___  " << std::endl
              << "   / __|  / _` | | . ` |  \\___ \\  / __| " << std::endl
              << "   \\__ \\ | (_| | | |\\  |  ____) | \\__ \\ " << std::endl
              << "   |___/  \\__,_| |_| \\_| |_____/  |___/ " << std::endl
              << "                                        " << std::endl
              << "   still another Navier - Stokes solver " << std::endl
              << "                                        " << std::endl;
        unsigned int dim = (unsigned int)std::atoi(argv[1]);
        std::string parameter_file = argv[2];
        if (dim == 2)
        {
            NavierStokes<2> ns(parameter_file);
            ns.solve();
        }
        else if (dim == 3)
        {
            NavierStokes<3> ns(parameter_file);
            ns.solve();
        }
        else
        {
            std::runtime_error("Invalid dimension. Only 2D and 3D simulations are supported.");
        }
    }
    catch (std::exception &exc)
    {
        std::cerr << std::endl << std::string(30,'-') << std::endl
                  << "Caught exception: " << std::endl
                  << exc.what() << std::endl
                  << "Aborting." << std::endl
                  << std::endl << std::string(30,'-') << std::endl;
        return 1;
    }
    catch (...)
    {
        std::cerr << std::endl << std::string(30,'-') << std::endl
                  << "Unknown exception. " << "Aborting." << std::endl
                  << std::endl << std::string(30,'-') << std::endl;
        return 1;
    }

    return 0;
}
