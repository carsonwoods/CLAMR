/* L7_Update_Perf - A performance test of the L7_Update Operation
 * from the L7 irregular communication library
 *
 * Background: L7 provides the abstraction of a global linear index space,
 * where each node owns a subset of the indices and can wants to receive a
 * subset of the indices owned by other nodes. The L7_Update operation
 * implements the irregular graph communication plan that moves remote
 * indices to each destination node. This operation is similar to one used
 * in a wide range of applications, including cell-based AMR codes (CLAMR/xRage),
 * particle transport codes (LAMMPS, Cabana proxy apps), linear solvers
 * (Trilinos/TPetra, hypre AMG), and many others. This benchmark seeks to
 * understand the performance of different implementations of that primitive
 * on different systems and system software packages.
 *
 * Benchmark structure:
 * 1. `numpes` processes, each of which owns `nowned` indices
 * 2. Each process requests a set of indices from `nneighbors` neighbors,
 *    equally divide between the nearest nodes above and below it in the
 *    rank space (the rank space is treated as circular - rank 0 with 4 neighbors
 *    will talk with ranks 1, 2, numpes-2, and numpes-1.
 * 3. Each process requests `nremote` indices in total from its neighbors, and
 *    the indices are requests in blocks of size `blocksz` with `stride` indexes
 *    between adjacent blocks. Note that this means that senders can either do
 *    contiguous or scattered sends.
 * 4. `niterations` L7_Update calls are performed, and for each update iteration the
 *    maximum time for any node to complete the update call and the effective
 *    communication bandwidth are computed. The benchmark reports the
 *    minimum, median, mean, and maximum update time and bandwidth for each
 *    iteration.
 * 5. Data being communicated can live either on the host, in CUDA memory,
 *    on OpenCL accelerator memory, or in OpenMP 4.5 accelerator memory.
 *
 * Default argument values:
 * 1. `typesize = 8` - We're moving around doubles
 * 2. `nowned = 2^28/typesize` - Each process has 256MB of data. Note that if this
 *    gets too large, we run the risk of running out of indexes.
 * 3. `nneighbors = sqrt(numpes)` - Each process has a modest number of neighbors
 *    that slowly increases as problem increases
 * 4. `nremote = nowned/64` - Each process needs to receive about 1% of the amount
 *    of data it owns
 * 4. `blocksz = nowned/16384` - Data communicated is *mostly* contiguous, but
 *    there are gaps
 * 5. `stride = 16` - Gaps are small
 * 6. `memspace = MEMSPACE_HOST` - data being communicated is on the host
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

// facilitates CLI-arg parsing
#include <getopt.h>

// Loads the MPI library of choice
#include <mpi.h>

// Loads the bundled L7 library
#include "l7/l7.h"

// If OpenCL is present, load the EZCL library
// The EZCL library provides OpenCL API endpoints
#ifdef HAVE_OPENCL
#include "ezcl/ezcl.h"
#endif

extern void initialize_data_host(void **odata, int nowned, int nremote, int type_size, int start);

// If CUDA is available, initialize data on CUDA device
#ifdef HAVE_CUDA
extern void initialize_data_cuda(void **odata, int nowned, int nremote, int type_size, int start);
#endif

// If OpenCL is available, initialize data on OpenCL device
#ifdef HAVE_OPENCL
extern void initialize_data_opencl(void **odata, int nowned, int nremote, int type_size, int start);
#endif

// If OpenMP is available and of sufficient version,
// initialize data on OpenCL device
#if defined(_OPENMP) && _OPENMP >= 201511
extern void initialize_data_openmp(void **odata, int nowned, int nremote, int type_size, int start);
#endif

// specifies the possible memory spaces (depends on support device)
// also maps members of the memspace enum to the memspace_t type
enum memspace {
  MEMSPACE_HOST,
  MEMSPACE_CUDA,
  MEMSPACE_OPENCL,
  MEMSPACE_OPENMP
};
typedef enum memspace memspace_t;

/*
 * We declare these settings as global
 * because they're used by basically every
 * function and passing them would be a
 * chore. I guess in theory we could make
 * a struct to contain them, but this is
 * good enough for this simple benchmark
*/

 /*
   CW: This is *probably* where the size of the data and the values can be modified.
   More specifically, this isn't where it will happen, but these values are what
   could be used to create irregular communication needs.

   For specific explanations for these values, see lines: 15-45 (benchmark structure section)
 */
static int typesize = 8;
static int numpes = 0;
static int niterations = 100;
static int nneighbors = -1;
static int nowned = -1;
static int nremote = -1;
static int blocksz = -1;
static int stride = -1;
static memspace_t memspace = MEMSPACE_HOST;

/*
 * these are options that are specified using the
 * type defined in the getopt.h file.
 *
 * TODO: I need to ask what is going on
 *       with this block of code
 * This might be helpful:
 * https://en.wikipedia.org/wiki/Getopt#Example_2_.28using_GNU_extension_getopt_long.29
*/
static struct option long_options[] = {
    /* These options set a flag. */
    {"typesize",   required_argument, 0, 't'},
    {"iterations", required_argument, 0, 'i'},
    {"owned",      required_argument, 0, 'o'},
    {"neighbors",  required_argument, 0, 'n'},
    {"remote",     required_argument, 0, 'r'},
    {"blocksize",  required_argument, 0, 'b'},
    {"stride",     required_argument, 0, 's'},
    {"memspace",   required_argument, 0, 'm'},
    {0, 0, 0, 0}
};

// offers a help string to define parameters in the CLI
void
usage(char *exename)
{
   fprintf(stderr, "usage: %s [-t typesize] [-i iterations] [-n neighbors] [-o owned] [-r remote] [-b blocksize] [-s stride] [-m memspace]\n", exename);
   exit(-1);
}

void parse_arguments(int argc, char **argv)
{
    int c;
    while (1) {
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "t:i:n:o:r:b:s:m:",
                       long_options, &option_index);
        if (c == -1) {
            break;
        }

        switch (c) {
        case 'i':
            // used to set custom iteration's value
            niterations = atoi(optarg);
            if (niterations < 0) usage(argv[0]);
            break;
        case 't':
            // used to set typesize value
            typesize = atoi(optarg);
            if (typesize < 1 || typesize > 8) usage(argv[0]);
            break;
        case 'n':
            // used to set nneighbors value
            nneighbors = atoi(optarg);
            if (nneighbors < 0) usage(argv[0]);
            break;
        case 'o':
            // used to set nowned value
            nowned = atoi(optarg);
            if (nowned < 0) usage(argv[0]);
            break;
        case 'r':
            // used to set nremote value
            nremote = atoi(optarg);
            if (nremote < 0) usage(argv[0]);
            break;
        case 'b':
            // used to set blocksz value
            blocksz = atoi(optarg);
            if (blocksz < 0) usage(argv[0]);
            break;
        case 's':
            // used to set stride size value
            stride = atoi(optarg);
            if (stride < 0) usage(argv[0]);
            break;
        case 'm':
            // used to set memory space
            // valid options include host, cuda, opencl, or openmp
            if (strcmp(optarg, "host") == 0) {
                memspace = MEMSPACE_HOST;
            } else if (strcmp(optarg, "cuda") == 0) {
                #ifdef HAVE_CUDA
                    memspace = MEMSPACE_CUDA;
                #else
                    fprintf(stderr, "No compiler support for CUDA accelerator memory\n");
                    usage(argv[0]);
                #endif
            } else if (strcmp(optarg, "opencl") == 0) {
                #ifdef HAVE_OPENCL
                    memspace = MEMSPACE_OPENCL;
                #else
                    fprintf(stderr, "No compiler support for OpenCL accelerator memory\n");
                    usage(argv[0]);
                #endif
            } else if (strcmp(optarg, "openmp") == 0) {
                #if defined(_OPENMP) && _OPENMP >= 201511
                    memspace = MEMSPACE_OPENMP;
                #else
                    fprintf(stderr, "No compiler support for OpenMP accelerator memory\n");
                    usage(argv[0]);
                #endif
            } else {
                fprintf(stderr, "Invalid memory space %s\n", optarg);
                usage(argv[0]);
            }
            break;
        case '?':
            usage(argv[0]);
            break;
        }
    }
    /* Now compute any unset values from the defaults */
    if (niterations < 0)
        niterations = 100;
    if (nneighbors < 0)
        nneighbors = sqrt(numpes);
    if (nowned < 0)
        nowned = (1<<28) / typesize;
    if (nremote < 0)
        nremote = nowned/(1 << 6);
    if (blocksz < 0) {
        int x = gauss_dist_int(1<<15, 10);
        printf("%d\n", x);
        blocksz = nowned/(1 << 15);
    }
    if (stride < 0)
        stride = 16;

    return;
}

// implementation of the following
// https://en.wikipedia.org/wiki/Box–Muller_transform
int gauss_dist_int(double mean, double stdev) {
    // generates two random numbers that form the seeds
    // of the transform
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;

    // generates the R and Theta values from the above
    // documentation
    double r = sqrt(-2.*log(u1));
    double theta = (2*M_PI*u2);

    // an additional number can be generated in the
    // same distribution using the alternate form
    // ((r*sin(theta)) * stdev) + mean

    printf("%f", (r*cos(theta)));

    return ( (r*cos(theta)) * stdev ) + mean;
}

// casts two void pointers (va, vb) into doubles (a, b)
// if a < b, return -1
// if a > b, return 1
// if a == b, return 0
int double_compare(const void *va, const void *vb)
{
    const double a = *(const double *)va, b = *(const double *)vb;
    if (a < b) return -1;
    if (a > b) return 1;
    return 0;
}

// casts two void pointers (va, vb) into integers (a, b)
// returns a - b
int int_compare(const void *va, const void *vb)
{
    const int a = *(const int *)va, b = *(const int *)vb;
    return a - b;
}

// returns L7_Datatype based on type_size integer
// based on L7_Datatype enum
enum  L7_Datatype typesize_to_l7type(int type_size)
{
    switch(type_size) {
        case 1:
            return L7_CHAR;
        case 2:
            return L7_SHORT;
        case 4:
            return L7_INT;
        case 8:
            return L7_DOUBLE;
        default:
            return -1;
   }
}

void
report_results_update(int penum, double *time_total_pe, int count_updated_pe, int num_timings, int type_size)
{
    int i, count_updated_global, bytes_updated, remainder;

    // pointer which will store all timing values
    double *time_total_global;
    // allocates enough space for num_timings # of doubles
    time_total_global = (double *)malloc(num_timings*sizeof(double));

    /* Gather iteration data from each node */
    /* The maximum iteration time */
    L7_Array_Max(time_total_pe, num_timings, L7_DOUBLE, time_total_global);
    /* The total number of items received */
    L7_Sum(&count_updated_pe, 1, L7_INT, &count_updated_global);
    bytes_updated = count_updated_global*type_size;

    // if the process is the first process
    // (sets logic to run on first process)
    if (penum == 0)
    {
        double *bandwidth_global = (double *)malloc(num_timings*sizeof(double));
        double latency_med, bandwidth_med;
        double latency_mean = 0.0, bandwidth_mean = 0.0;
        /* Compute the effective bandwidth for each iteration and the mean
         * time and bandwidth */
        // iterates through all timings taken
        for (i = 0; i < num_timings; i++) {
            if (time_total_global[i] != 0.0)
                // global bandwidth is calculated by bytes_updated/time_elapsed
                bandwidth_global[i] = bytes_updated / time_total_global[i];
            else
                // sets global bandwidth to 0 since cannot divide by zero
                bandwidth_global[i] = 0.0;

            // adds previously calculated latency (total time elapsed per iteration) to sum
            latency_mean += time_total_global[i];
            // adds previously calculated bandwidth to sum
            bandwidth_mean += bandwidth_global[i];
        }

        // calculates mean based on number of timings and previously calculated sums
        latency_mean /= num_timings;
        bandwidth_mean /= num_timings;

        /* Sort the arrays to compute the min, median, and max latency/bandwidth */
        qsort(time_total_global, num_timings, sizeof(double), double_compare);
        qsort(bandwidth_global, num_timings, sizeof(double), double_compare);

        /* Properly compute the median */
        if (num_timings % 2) {
            // runs if odd number of timings
            latency_med = time_total_global[num_timings/2];
            bandwidth_med = bandwidth_global[num_timings/2];
        } else {
            // runs if even number of timings (runs average on two middle values)
            latency_med = (time_total_global[num_timings/2]
                          + time_total_global[num_timings/2 - 1]) / 2;
            bandwidth_med = (bandwidth_global[num_timings/2]
                            + bandwidth_global[num_timings/2 - 1]) / 2;
        }

        /* Print results */
        printf("nPEs\tMem\tType\tnOwned\t\tnRemote\tBlockSz\tStride\tnIter");
        printf("\tLat(avg/min/med/max)\t\t\tBW(avg/min/med/max)\n");
        printf("%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,",
               numpes, memspace, typesize, nowned,
               nremote, blocksz, stride, num_timings);
        printf("\t%f/%f/%f/%f,\t%f/%f/%f/%f\n",
               latency_mean, time_total_global[0], latency_med,
               time_total_global[num_timings-1], bandwidth_mean,
               bandwidth_global[0], bandwidth_med, bandwidth_global[num_timings-1]);


        free(bandwidth_global);
    }

   /*
    * Testing complete
    */
   free(time_total_global);

   return;
}

int main(int argc, char *argv[])
{
    int penum = 0, ierr;

    int i, j, my_start_index,  remainder,
        num_partners_lo, num_partners_hi,
        offset, num_indices_per_partner, num_indices_offpe,
        inum, l7_id, gtime, count_updated_pe, iout;

    enum L7_Datatype l7type;

    // stores what time benchmark starts and stops (per process)
    double time_start, time_stop;

    // will store total time elapsed for benchmark (per process)
    double *time_total_pe;

    // will be array of partner process numbers
    int *partner_pe;
    int *needed_indices;

    void *data;

    // initialize L7
    // L7_Init(process_number, number_of_processes, argc, argv, do_quo_setup, lttrace_on)
    // I think L7_Init places process number in penum here
    ierr = L7_Init(&penum, &numpes, &argc, argv, 0, 0);

    // parse CLI arguments
    parse_arguments(argc, argv);

    // if compiled with OpenCL support
    #ifdef HAVE_OPENCL
    // attempts to initialize OpenCL on GPU
    ierr = ezcl_devtype_init(CL_DEVICE_TYPE_GPU);
    if (ierr == EZCL_NODEVICE) {
        ierr = ezcl_devtype_init(CL_DEVICE_TYPE_CPU);
    }

    // if EZCL fails to initialize, exit
    if (ierr != EZCL_SUCCESS) {
        printf("No opencl device available -- aborting\n");
        exit(-1);
    }

    // if OpenCL was initialized successfully,
    // initialize L7 device
    L7_Dev_Init();
    #endif

    // initialize pointer (type double)
    // to hold niteration # of doubles
    time_total_pe = (double *)malloc(niterations * sizeof(double));

    // sets the start index of each process
    // offset is process number * nowned
    my_start_index = penum * nowned;

    /* Compute which neighbors we talk to */
    // initialize pointer (type int)
    // to hold nneighbors # of ints
    partner_pe = (int *)malloc(nneighbors * sizeof(int));
    remainder = nneighbors % 2;

    // assigns each process a high and low process to
    // serve as a communication partner
    if (penum < (nneighbors /2) ) {
        num_partners_lo = nneighbors / 2;
        num_partners_hi = nneighbors / 2 + remainder;
    }
    else {
        num_partners_lo = nneighbors / 2 + remainder;
        num_partners_hi = nneighbors / 2;
    }

    /*
     * Indices below this PE - make sure to go in ascending PE order so the
     * indices are in ascending order!
    */
    offset = 0;
    for (i=1; i <= num_partners_lo; i++) {
        int partner = (i > penum) ? (numpes + penum - i) : penum - 1;
        //printf("[pe %d] offset %d penum %d i %d \n",penum, offset, penum, i);
        partner_pe[offset] = partner;
        offset++;
   }

    /* Indices above this PE */
    for (i=1; i<=num_partners_hi; i++){
        int partner = (i + penum >= numpes) ? i + penum - numpes : i + penum;
        //printf("[pe %d] offset %d penum %d i %d \n",penum, offset, penum, i);
        partner_pe[offset] = partner;
        offset++;
   }

    /*
     * Note that neighbors need to be an ascending order. We do this by sorting
     * our list of neighbors after computing it
    */
   qsort(partner_pe, nneighbors, sizeof(int), int_compare);

    if (nneighbors != 0) {
        num_indices_per_partner = nremote / nneighbors;
    }
    else {
        nremote = 0;
        num_indices_per_partner = 0;
    }

    /*
     * Generate needed indices
    */
    needed_indices = (int *)malloc(nremote * sizeof(int));
    num_indices_offpe = 0;
    for (i=0; i<nneighbors; i++) {
       int k;
       inum = 0;

       for (j=0, k = 0; j<num_indices_per_partner; j++, k++) {
            /* Detect end of block */
            if (k >= blocksz) {
                inum += (1 + stride);
                k = 0;
            } else {
                inum++;
            }

            /* Detect if we would walk off the end of the remote */
            if (inum >= nowned)
                break;

            needed_indices[num_indices_offpe] = partner_pe[i] * nowned + inum;
            num_indices_offpe++;
      }
   }

   /* Detect if nremote is actually less than we thought due to block/striding
    * and adjust appropriately so the eventual bandwidth calculation is right */
   nremote = num_indices_offpe;

   /*
    * Allocate data arrays on device and wait for initialization to complete
    */
    l7type = typesize_to_l7type(typesize);
    unsigned long data_size = nowned + nremote;
    switch (memspace) {
        case MEMSPACE_HOST:
            initialize_data_host(&data, nowned, nremote, typesize, my_start_index);
            break;
        #if defined(HAVE_CUDA) && defined(L7_CUDA_OFFLOAD)
        case MEMSPACE_CUDA:
            initialize_data_cuda(&data, nowned, nremote, typesize, my_start_index);
            break;
        #endif
        #ifdef HAVE_OPENCL
        case MEMSPACE_OPENCL:
            initialize_data_opencl(&data, nowned, nremote, typesize, my_start_index);
            break;
        #endif
        #if defined(_OPENMP) && defined(L7_OPENMP_OFFLOAD) && _OPENMP >= 201511
        case MEMSPACE_OPENMP:
            initialize_data_openmp(&data, nowned, nremote, typesize, my_start_index);
            break;
        #endif
        default:
            fprintf(stderr, "Unsupported memory space to initialize.\n");
            exit(-1);
            break;
    }

    //printf("[pe %d] Finished array initialization.\n", penum);

    /*
     * Register decomposition with L7
    */
    // This is the ID linked to the database
    l7_id = 0;

    /*
     * Register decomposition with L7
    */
    #ifdef HAVE_OPENCL
    /*
     * L7 uses a different operation for the opencl setup.
     * if OpenCL is used, the if statement is instantiated
    */
    if (memspace == MEMSPACE_OPENCL) {
        L7_Dev_Setup(0, my_start_index, nowned, needed_indices, nremote, &l7_id);
    } else
    #endif
    {
        L7_Setup(0, my_start_index, nowned, needed_indices, nremote, &l7_id);
    }

    /*
     * Begin updating data
     * This is what is being "benchmarked", how quickly the data amount
     * of data specified is being updated across all members of the database
    */
    gtime = 1;

    for (i=0; i<niterations; i++) {
        time_start = L7_Wtime();

        #ifdef HAVE_OPENCL
        /*
         * L7 uses a different operation for the opencl update which
         * explicitly copies the data to the host
         */
        if (memspace == MEMSPACE_OPENCL) {
            L7_Dev_Update(data, l7type, l7_id);
        } else
        #endif
        {
            L7_Update(data, l7type, l7_id);
        }

        time_stop = L7_Wtime();
        time_total_pe[i] = time_stop - time_start;
    }

    /*
     * Report results
    */
    count_updated_pe = nremote;

    // This involves L7 collectives, so we can't shut down L7 or MPI yet.
    report_results_update(penum, time_total_pe, count_updated_pe, niterations, typesize);

    // gracefully clean memory
    free(time_total_pe);
    free(partner_pe);
    free(needed_indices);

    // free's memory allocated through L7
    L7_Free(&l7_id);

    out:
       L7_Terminate();

   return 0;
}

