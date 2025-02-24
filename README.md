# Solar wind MHD simulator

Based on [PLUTO](https://plutocode.ph.unito.it/) 4.4-patch1.

## Authors

* Samvel Arutyunyan (boundary conditions)
* Alexander Kodukov (parallelization)
* Dmitry Pavlov (initial and boundary conditions)
* Maksim Subbotin (parallelization, scripting)

## Build

```
make
````

## Run

### Stationary background mode

1. Place a `bnd.nc` file in the `bnd` directory.

2. Run a simulation from 10 days before the starting date to the starting date:

```
mpirun -n 32 ./pluto -i pluto_b.ini
```

3. Run a simulation from the starting date for 5 days:

```
mpirun -n 32 ./pluto -i pluto.ini -restart 10
```

### Evolving background mode

1. Place a `bnd.nc` file, and also `bnd-1.nc`..`bnd-10.nc` in the `bnd` directory.

2. Run a simulation from 10 days before the starting date to the starting date:

```
mpirun -n 32 ./pluto -i pluto_b_daily.ini
```

3. Run a simulation from the starting date for 5 days:

```
mpirun -n 32 ./pluto -i pluto.ini -restart 10
```

## Results

Enjoy `dbl` files (internal format of PLUTO) and `vtk` files (openable in Paraview).

`data.0010.dbl` and `data.0010.vtk` are initial conditions. After
number 10, a file is saved every 1 hour of simulation, so
`data.0130.dbl` and `data.0130.vtk` are 5 days after the starting date.

The coordinate system is HEEQ.