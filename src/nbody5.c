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
    f32 *x;
    f32 *y;
    f32 *z;
    f32 *vx;
    f32 *vy;
    f32 *vz;

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

//
void pcopy(particle_t *p_ref, particle_t *p, u64 n)
{
    p->x = malloc(sizeof(f32) * n);
    p->y = malloc(sizeof(f32) * n);
    p->z = malloc(sizeof(f32) * n);
    p->vx = malloc(sizeof(f32) * n);
    p->vy = malloc(sizeof(f32) * n);
    p->vz = malloc(sizeof(f32) * n);
    for (u64 i = 0; i < n; i++)
    {
        p->x[i] = p_ref->x[i];
        p->y[i] = p_ref->y[i];
        p->z[i] = p_ref->z[i];

        //
        p->vx[i] = p_ref->vx[i];
        p->vy[i] = p_ref->vy[i];
        p->vz[i] = p_ref->vz[i];
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

delta compute_delta(particle_t *p_ref, particle_t *p, u64 n)
{
    delta d = {0.0, 0.0, 0.0};

    for (u64 i = 0; i < n; i++)
    {
        d.x += (p_ref->x[i] - p->x[i]);
        d.y += (p_ref->y[i] - p->y[i]);
        d.z += (p_ref->z[i] - p->z[i]);
    }

    d.x /= (f64)n;
    d.y /= (f64)n;
    d.z /= (f64)n;

    return d;
}

//
void move_particles_ref(particle_t *p, const f32 dt, u64 n)
{
    // Used to avoid division by 0 when comparing a particle to itself
    const f32 softening = 1e-20;

    // For all particles
    for (u64 i = 0; i < n; i++)
    {
        //
        f32 fx = 0.0;
        f32 fy = 0.0;
        f32 fz = 0.0;

        // Newton's law: 17 FLOPs (Floating-Point Operations) per iteration
        for (u64 j = 0; j < n; j++)
        {
            // 3 FLOPs (Floating-Point Operations)
            const f32 dx = p->x[j] - p->x[i]; // 1 (sub)
            const f32 dy = p->y[j] - p->y[i]; // 2 (sub)
            const f32 dz = p->z[j] - p->z[i]; // 3 (sub)

            // Compute the distance between particle i and j: 6 FLOPs
            const f32 d_2 = (dx * dx) + (dy * dy) + (dz * dz) + softening; // 9 (mul, add)

            // 2 FLOPs (here, we consider pow to be 1 operation)
            const f32 d_3_over_2 = pow(d_2, 3.0 / 2.0); // 11 (pow, div)

            // Calculate net force: 6 FLOPs
            fx += dx / d_3_over_2; // 13 (add, div)
            fy += dy / d_3_over_2; // 15 (add, div)
            fz += dz / d_3_over_2; // 17 (add, div)
        }

        // Update particle velocities using the previously computed net force: 6 FLOPs
        p->vx[i] += dt * fx; // 19 (mul, add)
        p->vy[i] += dt * fy; // 21 (mul, add)
        p->vz[i] += dt * fz; // 23 (mul, add)
    }

    // Update positions: 6 FLOPs
    for (u64 i = 0; i < n; i++)
    {
        p->x[i] += dt * p->vx[i];
        p->y[i] += dt * p->vy[i];
        p->z[i] += dt * p->vz[i];
    }
}

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
    f64 speed_up = 0.0;

    // Steps to skip for warm up
    const u64 warmup = 3;

    //
    particle_t p;
    particle_t p_ref;

    //
    init(&p_ref, n);
    pcopy(&p_ref, &p, n);

    const u64 s = sizeof(f32) * n;

    printf("\n\033[1mTotal memory size:\033[0m %llu B, %llu KiB, %llu MiB\n\n", s, s >> 10, s >> 20);

    //
    printf("\033[1m%5s %10s %10s %8s\033[0m\n", "Step", "Time, s", "Interact/s", "GFLOP/s");
    fflush(stdout);

    const f32 h1 = (f32)(n) * (f32)(n);
    const f32 h2 = (17.0 * h1 + 6.0 * (f32)n + 6.0 * (f32)n) * 1e-9;
    //
    for (u64 i = 0; i < steps; i++)
    {
        // Measure
        const f64 start = omp_get_wtime();
        move_particles(&p, dt, n);
        const f64 end = omp_get_wtime();

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
    for (u64 i = 0; i < steps; i++)
    {
        const f64 start_ref = omp_get_wtime();
        move_particles_ref(&p_ref, dt, n);
        const f64 end_ref = omp_get_wtime();

        // Do not take warm up runs into account
        if (i >= warmup)
        {
            rate_ref += h2 / (f32)(end_ref - start_ref);
        }
        printf("%5llu %10.3e %10.3e %8.1f %s\n",
               i,
               (end_ref - start_ref),
               h1 / (end_ref - start_ref),
               h2 / (end_ref - start_ref),
               (i < warmup) ? "(warm up)" : "");

        fflush(stdout);
    }

    delta d = compute_delta(&p_ref, &p, n);
    speed_up = rate / rate_ref;

    printf("Delta x : %f\n", d.x);
    printf("Delta y : %f\n", d.y);
    printf("Delta z : %f\n", d.z);

    printf("Speed Up de %f\n", speed_up);

    //
    delete_particule_table(&p);
    delete_particule_table(&p_ref);

    //
    return 0;
}