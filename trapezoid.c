#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

// Function to evaluate the curve (y = f(x))
float f(float x) {
    return x * x ; // Example: y = x^2
}

// Function to compute the area of a trapezoid
float trapezoid_area(float a, float b, float d) {
    float area = 0;
    for (float x = a; x < b; x+=d) {
        area += f(x) + f(x+d);
    }
    return area * d / 2.0f;
}

int main(int argc, char** argv) {
    int rank, size;
    float a = 0.0f, b = 1.0f;
    int n;
    float start, end, local_area, total_area;
    double T1, Tp, parallel_start, parallel_end, speedup, efficiency;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        // Get the number of intervals from the user
        printf("Enter the number of intervals: ");
        scanf("%d", &n);

        // Measure the serial execution time (T1)
        double serial_start = MPI_Wtime();
        float serial_area = trapezoid_area(a, b, (b - a) / n);
        double serial_end = MPI_Wtime();
        T1 = serial_end - serial_start;
    }

    // Broadcast n to all processes
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Calculate parallel execution time (Tp)
    parallel_start = MPI_Wtime();
    float d = (b - a) / n;
    float region = (b - a) / size;
    start = a + rank * region;
    end = start + region;
    local_area = trapezoid_area(start, end, d);
    MPI_Reduce(&local_area, &total_area, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    parallel_end = MPI_Wtime();
    Tp = parallel_end - parallel_start;

    if (rank == 0) {
        speedup = T1 / Tp;
        efficiency = speedup / size;
        printf("Total area: %f\n", total_area);
        printf("T1 (serial time): %f\n", T1);
        printf("Tp (parallel time): %f\n", Tp);
        printf("Speedup: %f\n", speedup);
        printf("Efficiency: %f\n", efficiency);
    }

    MPI_Finalize();
    return 0;
}

