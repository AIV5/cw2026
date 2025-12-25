#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <omp.h>

typedef uint64_t codes_t;
#include "check.c"

enum state_t {START, RUN, STOP};

void print_time(double s) {
    double m = s / 60;
    double h = m / 60;
    double d = h / 24;
    if (s < 100) {
        printf("%.1f s", s);
    } else if (m < 100) {
        printf("%.1f m", m);
    } else if (h < 100) {
        printf("%.1f h", h);
    } else {
        printf("%.1f d", d);
    }
}

int main(int argc, char** argv) {
    if (argc != 2) {
        return 1;
    }
    FILE* fp = fopen(argv[1], "rb");
    if (!fp) {
        return 2;
    }
    uint32_t signature, size;
    fread(&signature, sizeof(signature), 1, fp);
    if (signature != SIGNATURE) {
        printf("Wrong signature: got %u, expected %u\n", signature, SIGNATURE);
        fclose(fp);
        return 1;
    }
    fread(&size, sizeof(size), 1, fp);
    codes_t * codes_array = calloc(size, sizeof(codes_t));
    fread(codes_array, sizeof(codes_t), size, fp);
    fclose(fp);
    printf("read %"PRIu32" codes\n", size);
    double start = omp_get_wtime();
    int solutions = 0;
    int i = 0, processed = 0;
    #pragma omp parallel shared(i, processed) reduction(+: solutions)
    {
        #pragma omp master
        printf("got %d threads\n", omp_get_num_threads());
        enum state_t state = START;
        int j0, j1;
        while (state != STOP) {
            #pragma omp critical
            {
                if (state != START) {
                    processed += j1 - j0;
                    double elapsed = omp_get_wtime() - start;
                    double ratio = (double) processed / size;
                    double time_estimate = elapsed / ratio;
                    printf("%.2f%%\t", ratio * 100.0);
                    print_time(elapsed);
                    printf(" + ");
                    print_time(time_estimate - elapsed);
                    printf(" = ");
                    print_time(time_estimate);
                    printf("          \r");
                    fflush(stdout);
                } else {
                    state = RUN;
                }
                if (i == size) {
                    state = STOP;
                    j0 = size;
                    j1 = size;
                } else {
                    j0 = i;
                    j1 = j0 + 16;
                    if (j1 > size) {
                        j1 = size;
                    }
                    i = j1;
                }
            }
            for (int j = j0; j < j1; ++j) {
	        int t = check(codes_array[j]);
                solutions += t;
                if (t) {
                    printf("found "PRIROWT"          \n", codes_array[j]);
                }
            }
        }
    }
    double end = omp_get_wtime();
    printf("processed in %.2f s          \nsolutions: %d\n", end - start, solutions);
    free(codes_array);
}
