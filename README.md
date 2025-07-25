```
                  _   _    _____        
                 | \ | |  / ____|       
    ___    __ _  |  \| | | (___    ___  
   / __|  / _` | | . ` |  \___ \  / __| 
   \__ \ | (_| | | |\  |  ____) | \__ \ 
   |___/  \__,_| |_| \_| |_____/  |___/ 
                                        
   still another Navier - Stokes solver 
                                        
```


This is a Finite Element solver for the Navier-Stokes equations.

saNSs implements Galerkin Least Squares stabilisation for the time-dependent equations, and employs an implicit Euler scheme for time integration. The nonlinear system is solved with a monolithic appoach, via Newton iterations.

More details as well as an example are provided in [details.md](details.md).

## Getting started
> [!IMPORTANT]
> saNSs is built on top of [deal.ii](https://github.com/dealii/dealii), hence to use this code you will need an installation of deal.ii as well.

To build saNSs executable you can follow these steps
``` bash
# git clone git@github.com:Gabriele-Codega/sanss.git
cd sanss

mkdir build
cd build

cmake ..

make
```

This will produce an executable in the root directory `sanss`.
> [!NOTE]
> CMake might default to Debug build type. You can change this behaviour by doing `cmake -DCMAKE_BUILD_TYPE=Release ..`. Compiling in release mode makes the code noticeably faster.

Once you built the executable, you can run the solver with
``` bash
mpirun -np <number_of_processes> ./sanss <spacedim> <parameter_file>
```
where `spacedim` is the spatial dimension of your problem (2 or 3), and `parameter_file` is a `.prm` file that contains parameters for the simulation. Currently the code can only be used on three 2-dimensional test cases, whose parameter files are provided in `examples/`.
