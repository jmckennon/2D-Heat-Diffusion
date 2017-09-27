# Simulation of Heat Diffusion using MPI Java

This program simulates heat diffusion over a 2D area in a parralelized manner, using MPI Java.

## To run:

Install MPI Java (FAQ: https://www.open-mpi.org/faq/?category=java)

Set up runtime environment with multiple CPUs/machines (works easily with VMs)

To run the program:

```
javac Heat2D_mpi.java
mpdboot -n <#cpus> -v
prunjava <#cpus> Heat2D_mpi <size> <max_time> <heat_time> <interval> > output.txt
mpdallexit
```
