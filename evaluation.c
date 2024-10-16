#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <unistd.h> 

double f(double x) {
    return x * x; 
}

int main(int argc, char** argv) {
    int rank, size, n = 0;
    double a = 0.0, b = 1.0, h, local_a, local_b;
    double local_integral = 0.0, total_integral = 0.0;
    int local_n;
    double start_time, end_time, execution_time;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        printf("Enter the number of intervals: ");
        fflush(stdout);  
        scanf("%d", &n);
        if (n <= 0) {
            printf("Number of intervals must be positive.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    start_time = MPI_Wtime();

    h = (b - a) / n; 
    local_n = n / size; 
    local_a = a + rank * local_n * h; 
    local_b = local_a + local_n * h;  

    if (rank == 0) {
        usleep(300000); 
    } else if (rank == 1) {
        usleep(200000); 
    } else if (rank == 2) {
        usleep(100000); 
    }

    local_integral = (f(local_a) + f(local_b)) / 2.0; 
    for (int i = 1; i < local_n; i++) {
        local_integral += f(local_a + i * h);
    }
    local_integral = local_integral * h;

    MPI_Reduce(&local_integral, &total_integral, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    end_time = MPI_Wtime();
    execution_time = end_time - start_time;

    if (rank == 0) {
        printf("The total area under the curve is: %f\n", total_integral);
        printf("Execution time: %f seconds\n", execution_time);
    }

    MPI_Finalize();
    return 0;
}
