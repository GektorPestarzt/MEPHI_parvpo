#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#define MAX_THREADS 20
#define NUMBER_OF_ATTEMPTS 100
#define ARRAY_SIZE 10000000

#define RANDOM_SEED 920215

#ifdef _OPENMP
double find_max_parallel(int *array, long array_size, int number_of_threads);
#endif  // _OPENMP

double find_max_linear(int *array, long array_size, int);

double *clocking_threads(int max_threads, int number_of_attempts, long array_size);

void output(FILE *fd, double *avarage_time_result, int thread_max, int number_of_attempts, long array_size);

int main() {
    int max_threads = 1;
#ifdef _OPENMP
    max_threads = MAX_THREADS;
#endif // _OPENMP

    double start = omp_get_wtime();

    // double *avarage_time_result;
    // avarage_time_result =
    clocking_threads(max_threads, NUMBER_OF_ATTEMPTS, ARRAY_SIZE);

    double end = omp_get_wtime();
    printf("Completed: %f\n", end - start);

    return 0;
}

double *check_different_threads(int max_threads, long array_size, double (*find_max)(int *, long, int));

double *clocking_threads(int max_threads, int number_of_attempts, long array_size) {
    double *avarage_time_result = (double *) calloc(max_threads, sizeof(double));

    for (int i = 0; i < number_of_attempts; ++i) {
        double *attempt;

#ifdef _OPENMP
        attempt = check_different_threads(max_threads, array_size, find_max_parallel);
#else 
        attempt = check_different_threads(max_threads, array_size, find_max_linear);
#endif  // _OPENMP

        for (int i = 0; i < max_threads; ++i) {
            avarage_time_result[i] += attempt[i];
        }

        free(attempt);
    }

    for (int i = 0; i < max_threads; ++i) {
        avarage_time_result[i] /= number_of_attempts;
    }

    return avarage_time_result;
}

double *check_different_threads(int max_threads, long array_size, double (*find_max)(int *, long, int)) {
    /* Initialize the RNG */
    srand(RANDOM_SEED);

    int *array = (int *) malloc(sizeof(int) * array_size);
    for (long i = 0; i < array_size; ++i) {
        array[i] = rand();
    }

    double *time_result = (double *) calloc(max_threads, sizeof(double));
    for (int threads = 1; threads <= max_threads; ++threads) {
        double time = find_max(array, array_size, threads);
        time_result[threads - 1] += time;
    }

    free(array);
    return time_result;
}

#ifdef _OPENMP
double find_max_parallel(int *array, long array_size, int number_of_threads) {
    int  max   = -1;

    double start = omp_get_wtime();
    // #pragma omp parallel num_threads(number_of_threads) shared(array, array_size) reduction(max: max) default(none)
    #pragma omp parallel num_threads(number_of_threads) shared(array, array_size) reduction(max: max) default(none)
    {
        #pragma omp for
        for (long i = 0; i < array_size; ++i)
        {
            if (array[i] > max) { max = array[i]; };
        }
    }

    #pragma omp barrier
    double end = omp_get_wtime();

    return end - start;
}
#endif  // _OPENMP

double find_max_linear(int *array, long array_size, int number_of_threads) {
    number_of_threads += number_of_threads; // unused parametr

    int  max   = -1;
    clock_t start, end;

    start = clock();
    for (long i = 0; i < array_size; ++i)
    {
        if (array[i] > max) { max = array[i]; };
    }
    end = clock();

    double execution_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    return execution_time;
}