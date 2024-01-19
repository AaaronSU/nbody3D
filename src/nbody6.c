//
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "utils.c"
#include <math.h>

//
typedef struct
{
    f32 *restrict x;
    f32 *restrict y;
    f32 *restrict z;
    f32 *restrict vx;
    f32 *restrict vy;
    f32 *restrict vz;

} particle_t;

//
void init(particle_t *p, u64 n)
{
    p->x = malloc(sizeof(f32) * n);
    p->y = malloc(sizeof(f32) * n);
    p->z = malloc(sizeof(f32) * n);
    p->vx = malloc(sizeof(f32) * n);
    p->vy = malloc(sizeof(f32) * n);
    p->vz = malloc(sizeof(f32) * n);
    for (u64 i = 0; i < n; i++)
    {
        //
        u64 r1 = (u64)rand();
        u64 r2 = (u64)rand();
        f32 sign = (r1 > r2) ? 1 : -1;

        //
        p->x[i] = sign * (f32)rand() / (f32)RAND_MAX;
        p->y[i] = (f32)rand() / (f32)RAND_MAX;
        p->z[i] = sign * (f32)rand() / (f32)RAND_MAX;

        //
        p->vx[i] = (f32)rand() / (f32)RAND_MAX;
        p->vy[i] = sign * (f32)rand() / (f32)RAND_MAX;
        p->vz[i] = (f32)rand() / (f32)RAND_MAX;
    }
}

void delete_particule_table(particle_t *p)
{
    free(p->x);
    free(p->y);
    free(p->z);
    free(p->vx);
    free(p->vy);
    free(p->vz);
}

//
void move_particles(particle_t *p, const f32 dt, u64 n)
{
    // Used to avoid division by 0 when comparing a particle to itself
    const f32 softening = 1e-20;

// For all particles
#pragma omp parallel for simd
    for (u64 i = 0; i < n; i++)
    {
        //
        f32 fx = 0.0;
        f32 fy = 0.0;
        f32 fz = 0.0;

        f32 x_i = p->x[i];
        f32 y_i = p->y[i];
        f32 z_i = p->z[i];

        // Newton's law: 17 FLOPs (Floating-Point Operations) per iteration
        for (u64 j = 0; j < n; j++)
        {
            // 3 FLOPs (Floating-Point Operations)
            const f32 dx = p->x[j] - x_i; // 1 (sub)
            const f32 dy = p->y[j] - y_i; // 2 (sub)
            const f32 dz = p->z[j] - z_i; // 3 (sub)

            // Compute the distance between particle i and j: 6 FLOPs
            const f32 d_2 = (dx * dx) + (dy * dy) + (dz * dz) + softening; // 9 (mul, add)

            // 2 FLOPs (here, we consider pow to be 1 operation)
            const f32 d_3_over_2 = 1 / (d_2 * (f32)(*(int *)&d_2 >> 1)); // 11 (pow, div)

            // Calculate net force: 6 FLOPs
            fx += dx * d_3_over_2; // 13 (add, div)
            fy += dy * d_3_over_2; // 15 (add, div)
            fz += dz * d_3_over_2; // 17 (add, div)
        }

        // Update particle velocities using the previously computed net force: 6 FLOPs
        p->vx[i] += dt * fx; // 19 (mul, add)
        p->vy[i] += dt * fy; // 21 (mul, add)
        p->vz[i] += dt * fz; // 23 (mul, add)
    }
    // Update positions: 6 FLOPs
#pragma omp parallel for simd
    for (u64 i = 0; i < n; i++)
    {
        p->x[i] += dt * p->vx[i];
        p->y[i] += dt * p->vy[i];
        p->z[i] += dt * p->vz[i];
    }
}

//
int main(int argc, char **argv)
{
    // Number of particles to simulate
    const u64 n = (argc > 1) ? atoll(argv[1]) : 16384;

    // Number of experiments
    const u64 steps = 13;

    // Time step
    const f32 dt = 0.01;

    //
    f64 rate = 0.0, drate = 0.0;
    f64 rate_ref = 0.0;
    f64 speed_up;

    // Steps to skip for warm up
    const u64 warmup = 3;

    //
    particle_t p;

    //
    init(&p, n);

    const u64 s = sizeof(f32) * n;

    printf("\n\033[1mTotal memory size:\033[0m %llu B, %llu KiB, %llu MiB\n\n", s, s >> 10, s >> 20);

    //
    printf("\033[1m%5s %10s %10s %8s\033[0m\n", "Step", "Time, s", "Interact/s", "GFLOP/s");
    fflush(stdout);

    //
    for (u64 i = 0; i < steps; i++)
    {
        // Measure
        const f64 start = omp_get_wtime();
        move_particles(&p, dt, n);
        const f64 end = omp_get_wtime();

        // Number of interactions/iteration
        const f32 h1 = (f32)(n) * (f32)(n);

        // Number of GFLOPs
        // Innermost loop (Newton's law)   : 17 FLOPs x n (innermost trip count) x n (outermost trip count)
        // Velocity update (outermost body):  6 FLOPs x n (outermost trip count)
        // Positions update                :  6 FLOPs x n
        const f32 h2 = (17.0 * h1 + 6.0 * (f32)n + 6.0 * (f32)n) * 1e-9;

        // Do not take warm up runs into account
        if (i >= warmup)
        {
            rate += h2 / (f32)(end - start);
            drate += (h2 * h2) / (f32)((end - start) * (end - start));
        }

        //
        printf("%5llu %10.3e %10.3e %8.1f %s\n",
               i,
               (end - start),
               h1 / (end - start),
               h2 / (end - start),
               (i < warmup) ? "(warm up)" : "");

        fflush(stdout);
    }
    // Average GFLOPs/s
    rate /= (f64)(steps - warmup);

    // Deviation in GFLOPs/s
    drate = sqrt(drate / (f64)(steps - warmup) - (rate * rate));

    printf("-----------------------------------------------------\n");
    printf("\033[1m%s %4s \033[42m%10.1lf +- %.1lf GFLOP/s\033[0m\n",
           "Average performance:", "", rate, drate);
    printf("-----------------------------------------------------\n");

    //
    delete_particule_table(&p);

    //
    return 0;
}