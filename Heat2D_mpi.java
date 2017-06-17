//Author: Jake McKennon
//Heat Diffusion using MPI Java
//Last Edited: 16 June 2017

// -----------------------------------------------------------------------------
// This program is designed to simulate 2D heat diffusion by producing
// a text-based output. It implements the Forward Euler method to accomplish the
// mathematical aspects of the simulation. It is optimized for a parallel system
// and uses MPI Java to allow the CPUs to communicate with each other.
// "CPU 0" is considered the master CPU and will typically do slightly more
// work than other CPUs.
//
// usage: "prunjava #cpus Heat2D_mpi size max_time heat_time interval"
//
// -----------------------------------------------------------------------------

import java.util.*;
import mpi.*;

public class Heat2D_mpi {
    private static double a = 1.0;          // heat speed
    private static double dt = 1.0;         // time quantum
    private static double dd = 2.0;         // change in system
    private static int size = 0;            // size of 2d space
    private static int max_time = 0;        // total time simulation runs
    private static int heat_time = 0;       // total time heating during sim
    private static int interval = 0;        // interval at which to print

    static int nprocs = 0;                  // total processes in the system
    final static int mtype = 0;             // type tag for Send/Recv

    public static void main(String[] args) throws MPIException {
        size = Integer.parseInt(args[1]);
        max_time = Integer.parseInt(args[2]);
        heat_time = Integer.parseInt(args[3]);
        interval  = Integer.parseInt(args[4]);

        MPI.Init(args);

        // get rank and system size
        int myrank = MPI.COMM_WORLD.Rank();
        nprocs = MPI.COMM_WORLD.Size();

        double r = a * dt / ( dd * dd );

        // create a space
        double[][][] z = new double[2][size][size];
        for(int p = 0; p < 2; p++)
            for(int x = 0; x < size; x++)
                for(int y = 0; y < size; y++)
                    z[p][x][y] = 0.0; // no heat or cold

        // create rows and offsets, finds any extra rows and gives them to rank 0
        int averows = size / nprocs;
        int rows = size / nprocs;
        int extra = size % nprocs;
        int offset = 0;
        int dSize = 2 * size * size;
        if(myrank == 0) {
            rows += extra;
        }
        else {
            offset = (myrank * averows) + extra;
        }

        Date startTime = new Date();

        // simulate heat diffusion
        for(int t = 0; t < max_time; t++) {
            int p = t % 2; // p = 0 or 1: indicates the phase

            // two left-most and two right-most columns are identical
            for(int y = 0; y < size; y++) {
                z[p][0][y] = z[p][1][y];
                z[p][size - 1][y] = z[p][size - 2][y];
            }

            // two upper and lower rows are identical
            for(int x = 0; x < size; x++) {
                z[p][x][0] = z[p][x][1];
                z[p][x][size - 1] = z[p][x][size - 2];
            }

            // keep heating the bottom until t < heat_time
            if(t < heat_time) {
                for(int x = size /3; x < size / 3 * 2; x++)
                z[p][x][0] = 19.0; // heat
            }

            if(nprocs > 1) { //only do send/receives if more than one node
                //arrays to send and receive boundary data
                double[] leftBoundarySend = new double[size];
                double[] rightBoundarySend = new double[size];
                double[] leftBoundaryRecv= new double[size];
                double[] rightBoundaryRecv = new double[size];

                for(int i = 0; i < size; i++) {
                    leftBoundarySend[i] = z[p][offset][i];
                    rightBoundarySend[i] = z[p][offset + rows - 1][i];
                }

                //odds and evens send and receive, and vice versa
                for(int i = 0; i < 2; i++) {
                    if(myrank % 2 == i) {
                        if(myrank != 0) {
                            MPI.COMM_WORLD.Send(leftBoundarySend,
                            0,
                            size,
                            MPI.DOUBLE,
                            myrank - 1,
                            mtype);
                        }
                        if(myrank != (nprocs - 1)) {
                            MPI.COMM_WORLD.Send(rightBoundarySend,
                            0,
                            size,
                            MPI.DOUBLE,
                            myrank + 1,
                            mtype);
                        }
                    }
                    else {
                        if(myrank != (nprocs - 1)) {
                            MPI.COMM_WORLD.Recv(rightBoundaryRecv,
                            0,
                            size,
                            MPI.DOUBLE,
                            myrank + 1,
                            mtype);
                        }
                        if(myrank != 0) {
                            MPI.COMM_WORLD.Recv(leftBoundaryRecv,
                            0,
                            size,
                            MPI.DOUBLE,
                            myrank - 1,
                            mtype);
                        }
                    }
                }

                // display intermediate results
                if(interval != 0 && (t % interval == 0 || t == max_time - 1)) {
                    //ranks > 0 send their data to rank 0 for printing
                    double[] stripe = new double[dSize];
                    if(myrank != 0) {
                        for(int x = offset; x < offset + averows; x++) {
                            for(int y = 0; y < size; y++) {
                                stripe[p * size * size + x * size + y] =
                                    z[p][x][y];
                            }
                        }
                        MPI.COMM_WORLD.Send(stripe,
                        0,
                        dSize,
                        MPI.DOUBLE,
                        0,
                        mtype);
                    }
                    else { // rank 0 receives
                        for(int i = 1; i < nprocs; i++) {
                            MPI.COMM_WORLD.Recv(stripe,
                            0,
                            dSize,
                            MPI.DOUBLE,
                            i,
                            mtype);

                            for(int x = (i * averows) + extra;
                            x < (i * averows) + averows + extra;
                            x++) {
                                for(int y = 0; y < size; y++) {
                                    z[p][x][y] =
                                    stripe[p * size * size + x * size + y];
                                }
                            }
                        }
                    }
                    if(myrank == 0) { //only rank 0 prints
                        System.out.println("time = " + t);
                        for(int y = 0; y < size; y++) {
                            for(int x = 0; x < size; x++)
                            System.out.print((int)(Math.floor(z[p][x][y] / 2)) +
                                " ");
                            System.out.println();
                        }
                        System.out.println();
                    }
                }

                int p2 = (p + 1) % 2;
                //calculate left boundaries
                if(myrank != 0) {
                    for(int y = 1; y < size - 1; y++) {
                        z[p2][offset][y] = z[p][offset][y] +
                        r * (z[p][offset + 1][y] -
                        2 * z[p][offset][y] + leftBoundaryRecv[y]) +
                        r * (z[p][offset][y + 1] -
                        2 * z[p][offset][y] + z[p][offset][y - 1]);
                    }
                }
                //calculate right boundaries
                if(myrank != (nprocs - 1)) {
                    for(int y = 1; y < size - 1; y++) {
                        z[p2][offset + rows - 1][y] =
                        z[p][offset + rows - 1][y] +
                        r * (rightBoundaryRecv[y] -
                        2 * z[p][offset + rows - 1][y] +
                        z[p][offset + rows - 1 - 1][y]) +
                        r * (z[p][offset + rows - 1][y + 1] -
                        2 * z[p][offset + rows - 1][y] +
                        z[p][offset + rows - 1][y - 1]);
                    }
                }
                //calculate the rest of the data
                for(int x = offset + 1; x < offset + rows - 1; x++) {
                    for(int y = 1; y < size - 1; y++) {
                        z[p2][x][y] = z[p][x][y] +
                        r * (z[p][x + 1][y] - 2 * z[p][x][y] + z[p][x - 1][y]) +
                        r * (z[p][x][y + 1] - 2 * z[p][x][y] + z[p][x][y - 1]);
                    }
                }
            }
            else {//only one node, so dont do send/recv
                // display intermediate results
                if(interval != 0 && (t % interval == 0 || t == max_time - 1)) {
                    System.out.println( "time = " + t );
                    for(int y = 0; y < size; y++) {
                        for(int x = 0; x < size; x++)
                        System.out.print(
                            (int)(Math.floor(z[p][x][y] / 2)) + " ");
                        System.out.println();
                    }
                    System.out.println();
                }

                // perform forward Euler method
                int p2 = (p + 1) % 2;
                for(int x = 1; x < size - 1; x++)
                for(int y = 1; y < size - 1; y++)
                z[p2][x][y] = z[p][x][y] +
                r * ( z[p][x+1][y] - 2 * z[p][x][y] + z[p][x-1][y] ) +
                r * ( z[p][x][y+1] - 2 * z[p][x][y] + z[p][x][y-1] );
            }
        } // end of simulation

        // finish the timer
        if(myrank == 0) {
            Date endTime = new Date();
            System.out.println("Elapsed time = " +
                (endTime.getTime() - startTime.getTime()));
        }

        MPI.Finalize();
    }
}
