/*=========================================================/
/  Comms routines for PX457 assignment 4. Contains all     /
/  routines which interact with MPI libraries. Many of     /
/  these are currently incomplete and will work only in    /
/  serial. You will need to correct this.                  /
/                                                          /
/  Original code created by B. Morgan  - November 2025     /
/  (based on previous code by N. Hine, D. Quigley)         /
/=========================================================*/
#include "comms.h"

#include "grid.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

int p;       /* Number of processor       */
int my_rank; /* Rank of current processor */

MPI_Comm cart_comm; /* Cartesian communicator    */

/* Coordinates of current rank in the processor grid */
int my_rank_coords[2];

/* Ranks of neighbours to the current processor (left, right, down, up) */
int my_rank_neighbours[4];

/* Time and initialisation and shutdown */
double t1, t2;

void comms_initialise(int argc, char** argv)
{
    /*==================================================================/
    / Function to initialise MPI, get the communicator size p and the   /
    / rank my_rank of the current process within that communicator.     /
    /-------------------------------------------------------------------/
    / N. Hine (based on code by D. Quigley) - Univ. Warwick             /
    /==================================================================*/
    int proot; /* square root of p */

    /* Remove these lines once you have added calls to initialise MPI */
    /* and to retreive my_rank and p, the number of processors        */
    MPI_Init(&argc, &argv); // Initialise the MPI environment
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); // Get processor rank
    MPI_Comm_size(MPI_COMM_WORLD, &p); // Get total processor num.

    /* Start the timer. Set t1 using MPI_Wtime() which returns a double */
    t1 = MPI_Wtime();
    /* Check that we have a square number of processors */
    proot = (int)sqrt((double)p + 0.5);
    if (proot * proot != p)
    {
        if (my_rank == 0)
        {
            printf("Number of processors must be an exact square!\n");
            exit(EXIT_FAILURE);
        }
    }

    return;
}

void comms_processor_map()
{
    /*====================================================================/
    / Function to map our p processors into a 2D Cartesian grid of      /
    / dimension proot by proot where proot = sqrt(p).                   /
    /                                                                   /
    / Should populate the arrays my_rank_cooords, which contains the    /
    / location of the current MPI task within the processor grid, and   /
    / my_rank_neighbours, which contains (in the order left, right,     /
    / down and up) ranks of neighbouring MPI tasks on the grid with     /
    / which the current task will need to communicate.                  /
    /-------------------------------------------------------------------/
    / N. Hine (based on code by D. Quigley) - University of Warwick     /
    /==================================================================*/

    /* Information for setting up a Cartesian communicator */
    int ndims = 2;
    int reorder = 0;
    int pbc[2] = {1, 1};
    int dims[2];

    /* Local variables */
    int proot; /* square root of p */

    /* Square root of number of processors */
    proot = (int)sqrt((double)p + 0.5);

    /* Dimensions of Cartesian communicator */
    dims[x] = proot;
    dims[y] = proot;

    /* Remove these lines when you have added appropriate MPI calls to
       create a Cartesian communicator from MPI_COMM_WORLD, stored the
       coordinates of the current MPI task in my_rank_coords, and set
       the array my_rank_neighbours to hold the tank of neighbouring
       tasks in the directions left, right, down and up.               */
        // create Cartesian communicator
        MPI_Cart_create(MPI_COMM_WORLD, 2, dims, pbc, reorder, &cart_comm);

        // get processor coords.
        MPI_Cart_coords(cart_comm, my_rank, 2, my_rank_coords);

        // find neighbours: left/right/down/up 
        MPI_Cart_shift(cart_comm, 0, 1, &my_rank_neighbours[left],
                        &my_rank_neighbours[right]);

        MPI_Cart_shift(cart_comm, 1, 1, &my_rank_neighbours[down],
                    &my_rank_neighbours[up]);

        printf("Rank %d coords: (%d, %d) neighbours: Left=%d Right=%d Down=%d Up=%d\n",
            my_rank,my_rank_coords[x], my_rank_coords[y],
            my_rank_neighbours[left], my_rank_neighbours[right],
            my_rank_neighbours[down], my_rank_neighbours[up]);
        

    return;
}

void comms_get_global_mag(double local_mag, double* global_mag)
{
    /*==================================================================/
    / Function to compute the glocal magnetisation of the grid by       /
    / averaging over all values of local_mag, and storing the result    /
    / in global_mag.                                                    /
    /-------------------------------------------------------------------/
    / N. Hine (based on code by D. Quigley) - University of Warwick     /
    /==================================================================*/

    /* This is only correct on one processor. You will need        */
    /* to use a collective communication routine to correct this.  */
    
    MPI_Reduce(&local_mag, global_mag, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    /* Insert collective communication operation here */
    if (my_rank == 0) {
        *global_mag = *global_mag / (double)p;
    }
}

void comms_halo_swaps()
{
    /*==================================================================/
    / Function to send boundary spins on each side of the local grid    /
    / to neighbour processors, and to receive from those processors the /
    / halo information needed to perform computations involving spins   /
    / on the boundary processors grid.                                  /
    /-------------------------------------------------------------------/
    / N. Hine (based on code by D. Quigley) - University of Warwick     /
    /==================================================================*/

    /* Send and receive buffers */
    int *sendbuf, *recvbuf;

    /* MPI Status */
    MPI_Status status;

    int ix, iy; /* Loop counters */

    /* If running on 1 processor copy boundary elements into opposite halo */
    if (p == 1)
    {
        for (iy = 0; iy < grid_domain_size; iy++)
        {
            grid_halo[right][iy] = grid_spin[iy][0];
            grid_halo[left][iy] = grid_spin[iy][grid_domain_size - 1];
        }

        for (ix = 0; ix < grid_domain_size; ix++)
        {
            grid_halo[up][ix] = grid_spin[0][ix];
            grid_halo[down][ix] = grid_spin[grid_domain_size - 1][ix];
        }

        return; /* Do not do any comms */
    }

    /* Allocate buffers */
    sendbuf = (int*)malloc(grid_domain_size * sizeof(int));
    if (sendbuf == NULL)
    {
        printf("Error allocating sendbuf in comms_halo_swaps\n");
        exit(EXIT_FAILURE);
    }
    recvbuf = (int*)malloc(grid_domain_size * sizeof(int));
    if (recvbuf == NULL)
    {
        printf("Error allocating recvbuf in comms_halo_swaps\n");
        exit(EXIT_FAILURE);
    }

    /* Send left hand boundary elements of grid_spin to my_rank_neighbours[left]
       and receive from my_rank_neighbours[right] into the appropriate part
       of grid_halo. Remember to use the appropriate communicator. */
    for (iy = 0; iy < grid_domain_size; iy++) {
        sendbuf[iy] = grid_spin[iy][0];
    }
    /* Insert MPI calls here to implement this swap. Use sendbuf and recvbuf */
    MPI_Sendrecv(sendbuf, grid_domain_size, MPI_INT, my_rank_neighbours[left], 0,
                 recvbuf, grid_domain_size, MPI_INT, my_rank_neighbours[right], 0,
                 cart_comm, &status);

    for (iy = 0; iy < grid_domain_size; iy++) {
        grid_halo[right][iy] = recvbuf[iy];
    }
    /* Send right hand boundary elements of grid_spin to my_rank_neighbours[right]
       and receive from my_rank_neighbours[left] into the appropriate part
       of grid_halo. Remember to use the appropriate communicator. */
    for (iy = 0; iy < grid_domain_size; iy++) {
        sendbuf[iy] = grid_spin[iy][grid_domain_size - 1];
    }
    /* Insert MPI calls here to implement this swap. Use sendbuf and recvbuf */
    MPI_Sendrecv(sendbuf, grid_domain_size, MPI_INT, my_rank_neighbours[right], 1,
                 recvbuf, grid_domain_size, MPI_INT, my_rank_neighbours[left], 1,
                 cart_comm, &status);

    for (iy = 0; iy < grid_domain_size; iy++) {
        grid_halo[left][iy] = recvbuf[iy];
    }
    /* Send bottom boundary elements of grid_spin to my_rank_neighbours[down]
       and receive from my_rank_neighbours[up] into the appropriate part
       of grid halo. Remember to use the appropriate communicator.  */

    /* Insert MPI calls here to implement this swap. Use sendbuf and recvbuf */
    for (ix = 0; ix < grid_domain_size; ix++) {
        sendbuf[ix] = grid_spin[grid_domain_size - 1][ix];
    }

    MPI_Sendrecv(sendbuf, grid_domain_size, MPI_INT, my_rank_neighbours[down], 2,
                 recvbuf, grid_domain_size, MPI_INT, my_rank_neighbours[up], 2,
                 cart_comm, &status);

    for (ix = 0; ix < grid_domain_size; ix++) {
        grid_halo[up][ix] = recvbuf[ix];
    }
    /* Send top boundary elements of grid_spin to my_rank_neighbours[up]
       and receive from my_rank_neighbours[down] into the appropriate part
       of grid halo. Remember to use the appropriate communicator. */
    for (ix = 0; ix < grid_domain_size; ix++) {
        sendbuf[ix] = grid_spin[0][ix];
    }
    /* Insert MPI call or calls here to implement this swap. Use sendbuf and recvbuf */
    MPI_Sendrecv(sendbuf, grid_domain_size, MPI_INT, my_rank_neighbours[up], 3,
                 recvbuf, grid_domain_size, MPI_INT, my_rank_neighbours[down], 3,
                 cart_comm, &status);

    for (ix = 0; ix < grid_domain_size; ix++) {
        grid_halo[down][ix] = recvbuf[ix];
    }
    /* Release memory */
    free(sendbuf);
    free(recvbuf);

    return;
}

void comms_get_global_grid()
{
    /*==================================================================/
    / Function to collect all contributions to the global grid onto     /
    / rank zero for visualisation.                                      /
    /-------------------------------------------------------------------/
    / N. Hine (based on code by D. Quigley) - University of Warwick     /
    /==================================================================*/

    /* comms buffer */
    int* combuff;

    /* MPI Status */
    MPI_Status status;

    /* Information on the remote domain */
    int remote_domain_start[2] = {0, 0};

    /* Loop counters and error flag */
    int ix, iy, ixg, iyg, ip;

    /* Just point at local grid if running on one processor */
    if (p == 1)
    {
        global_grid_spin = grid_spin;
        return;
    }

    if (my_rank == 0)
    {
        /* Rank 0 first fills out its part of the global grid */
        for (iy = 0; iy < grid_domain_size; iy++)
        {
            for (ix = 0; ix < grid_domain_size; ix++)
            {
                /* Global indices */
                ixg = ix + grid_domain_start[x];
                iyg = iy + grid_domain_start[y];

                global_grid_spin[iyg][ixg] = grid_spin[iy][ix];
            }
        }
    }

    /* Remove the following line when you have inserted appropriate MPI calls below */

    /* Allocate buffer */
    combuff = (int*)malloc(grid_domain_size * sizeof(int));
    if (combuff == NULL)
    {
        printf("Error allocating combuff in comms_get_global_grid\n");
        exit(EXIT_FAILURE);
    }
    const int TAG_START = 0;
    const int TAG_ROW = 1;

    if (my_rank == 0)
    {
        /* Now loops over all other ranks receiving their data */
        for (ip = 1; ip < p; ip++)
        {
            /* First receive remote_domain_start from rank ip */
            /* Insert an appropriate MPI call here */
            MPI_Recv(remote_domain_start, 2, MPI_INT, ip, TAG_START, MPI_COMM_WORLD, &status);
            /* Loop over rows within a domain */
            for (iy = 0; iy < grid_domain_size; iy++)
            {
                /* Receive this row from rank ip */
                /* Insert appropriate MPI call here */
                MPI_Recv(combuff, grid_domain_size, MPI_INT, ip,
                         TAG_ROW, MPI_COMM_WORLD, &status);
                for (ix = 0; ix < grid_domain_size; ix++)
                {
                    /* Global indices */
                    ixg = ix + remote_domain_start[x];
                    iyg = iy + remote_domain_start[y];

                    /* Store in global_grid_spin */
                    global_grid_spin[iyg][ixg] = combuff[ix];

                } /* elements in row */

            } /* rows */

        } /* processors */
    }
    else
    {
        /* All other processors must send the data rank 0 needs */

        /* Send grid_domain_start to rank 0 */
        /* Insert appropriate MPI call here */
        MPI_Send(grid_domain_start, 2, MPI_INT, 0, TAG_START, MPI_COMM_WORLD);
        /* Loop over rows in the domain, sending them to rank 0*/
        for (iy = 0; iy < grid_domain_size; iy++)
        {
            /* Insert appropriate MPI call here */
            MPI_Send(grid_spin[iy], grid_domain_size, MPI_INT, 0,
                     TAG_ROW, MPI_COMM_WORLD);
        }
    }

    /* Free memory */
    free(combuff);

    return;
}

void comms_finalise()
{
    /*==================================================================/
    / Function to finalise MPI functionality and exit cleanly           /
    /-------------------------------------------------------------------/
    / N. Hine (based on code by D. Quigley) - University of Warwick     /
    /==================================================================*/

    /* Remove the following line when you have inserted appropriate MPI calls below */
    
    /* Measure the time t2 using MPI_Wtime() which returns a double */
    t2 = MPI_Wtime();
    if (my_rank == 0 && p > 1)
    {
        printf("Total time elapsed since MPI initialised :  %12.6f s\n", t2 - t1);
    }

    /* Shutdown MPI - insert appropriate call here */
    MPI_Finalize();
    return;
}
