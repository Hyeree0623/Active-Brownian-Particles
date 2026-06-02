/*
==========================================================
Active Brownian Particle (ABP) Simulation
==========================================================

Model:
    Overdamped Active Brownian Particles

Interactions:
    Weeks-Chandler-Andersen (WCA) potential

Boundary Conditions:
    Periodic Boundary Conditions (PBC)

Features:
    - Rotational diffusion
    - Self propulsion
    - Brownian noise
    - Cell-list acceleration
    - GRO trajectory output

==========================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

// -------------------- PARAMETERS --------------------
#define N_PARTICLES   5000
#define NDIMS         2

#define BOX_X 100.0
#define BOX_Y 100.0

#define DT 0.00001
#define N_STEPS 10000000

#define DIFF 0.05
#define DIFF_R 0.15

#define EPSILON 1.0
#define SIGMA   1.0
#define CUTOFF  1.1224620483

#define V0 5.0

#define CELL_SIZE 2.5
#define NX 40
#define NY 40

#define OUTFILE "coord.gro"
#define BUF_SIZE 8192

// -------------------- STRUCTS --------------------
typedef struct {
    double x, y;
} Particle;

typedef struct Node {
    int id;
    Particle p;
    struct Node *next;
} Node;

// -------------------- GLOBAL STATE --------------------
static Particle pos[N_PARTICLES];
static Particle pos_old[N_PARTICLES];
static Particle force[N_PARTICLES];
static Particle active[N_PARTICLES];
static double theta[N_PARTICLES];

static Node node_pool[200000];
static int node_ptr = 0;
static Node *cell[NX][NY];

// -------------------- UTIL --------------------
static double randn() {
    double u1, u2, s;
    do {
        u1 = 2.0 * rand() / RAND_MAX - 1.0;
        u2 = 2.0 * rand() / RAND_MAX - 1.0;
        s = u1*u1 + u2*u2;
    } while (s >= 1.0 || s == 0.0);
    return u1 * sqrt(-2.0 * log(s) / s);
}

static Node* new_node() {
    return (node_ptr < 200000) ? &node_pool[node_ptr++] : NULL;
}

// -------------------- INIT --------------------
static void init_positions() {
    for (int i = 0; i < N_PARTICLES; i++) {
        int ok;
        do {
            ok = 1;
            pos[i].x = BOX_X * rand() / RAND_MAX;
            pos[i].y = BOX_Y * rand() / RAND_MAX;

            for (int j = 0; j < i; j++) {
                double dx = pos[i].x - pos[j].x;
                double dy = pos[i].y - pos[j].y;

                dx -= BOX_X * round(dx / BOX_X);
                dy -= BOX_Y * round(dy / BOX_Y);

                if (dx*dx + dy*dy < 1.2) {
                    ok = 0;
                    break;
                }
            }
        } while (!ok);

        pos_old[i] = pos[i];
        theta[i] = 2.0 * M_PI * rand() / RAND_MAX;
    }
}

// -------------------- PBC --------------------
static void apply_pbc(Particle *p) {
    if (p->x < 0) p->x += BOX_X;
    if (p->x >= BOX_X) p->x -= BOX_X;
    if (p->y < 0) p->y += BOX_Y;
    if (p->y >= BOX_Y) p->y -= BOX_Y;
}

// -------------------- CELL LIST --------------------
static void clear_cells() {
    memset(cell, 0, sizeof(cell));
    node_ptr = 0;
}

static void insert(int i) {
    int cx = (int)(pos[i].x / CELL_SIZE);
    int cy = (int)(pos[i].y / CELL_SIZE);

    if (cx < 0 || cx >= NX || cy < 0 || cy >= NY) return;

    Node *n = new_node();
    if (!n) return;

    n->id = i;
    n->p = pos[i];
    n->next = cell[cx][cy];
    cell[cx][cy] = n;
}

// -------------------- FORCE --------------------
static inline void compute_LJ() {
    memset(force, 0, sizeof(force));

    int dx8[5] = {0,-1,-1,0,-1};
    int dy8[5] = {0, 0,-1,-1, 1};

    double rc2 = CUTOFF * CUTOFF;

    for (int cx = 0; cx < NX; cx++) {
        for (int cy = 0; cy < NY; cy++) {

            for (Node *a = cell[cx][cy]; a; a = a->next) {
                int i = a->id;

                for (int k = 0; k < 5; k++) {
                    int nx = (cx + dx8[k] + NX) % NX;
                    int ny = (cy + dy8[k] + NY) % NY;

                    for (Node *b = cell[nx][ny]; b; b = b->next) {
                        int j = b->id;
                        if (i >= j) continue;

                        double dx = pos[j].x - pos[i].x;
                        double dy = pos[j].y - pos[i].y;

                        dx -= BOX_X * round(dx / BOX_X);
                        dy -= BOX_Y * round(dy / BOX_Y);

                        double r2 = dx*dx + dy*dy;
                        if (r2 > rc2) continue;

                        double r = sqrt(r2);
                        double sr = SIGMA / r;
                        double sr6 = pow(sr, 6);
                        double sr12 = sr6 * sr6;

                        double f = 24.0 * EPSILON * (2*sr12 - sr6) / r;

                        double fx = f * dx / r;
                        double fy = f * dy / r;

                        force[i].x += fx;
                        force[i].y += fy;
                        force[j].x -= fx;
                        force[j].y -= fy;
                    }
                }
            }
        }
    }
}

// -------------------- ACTIVE BROWNIAN STEP --------------------
static void step() {
    double dr = sqrt(2.0 * DIFF * DT);
    double rot = sqrt(2.0 * DIFF_R * DT);

    for (int i = 0; i < N_PARTICLES; i++) {

        theta[i] += rot * randn();

        active[i].x = V0 * cos(theta[i]) * DT;
        active[i].y = V0 * sin(theta[i]) * DT;

        pos_old[i] = pos[i];

        double dx = DT * force[i].x + dr * randn() + active[i].x;
        double dy = DT * force[i].y + dr * randn() + active[i].y;

        pos[i].x += dx;
        pos[i].y += dy;

        apply_pbc(&pos[i]);
    }
}

// -------------------- OUTPUT --------------------
static void write_frame(FILE *fp, int step) {

    fprintf(fp, "t=%d\n", step);
    fprintf(fp, "%d\n", N_PARTICLES);

    for (int i = 0; i < N_PARTICLES; i++) {
        fprintf(fp, "%d MOL S %d %lf %lf %lf\n",
                i, i+1, pos[i].x, pos[i].y, 0.0);
    }

    fprintf(fp, "%lf %lf 0.0\n", BOX_X, BOX_Y);
}

// -------------------- MAIN --------------------
int main() {

    srand(time(NULL));

    FILE *fp = fopen(OUTFILE, "w");

    init_positions();

    for (int step_i = 0; step_i < N_STEPS; step_i++) {

        clear_cells();

        for (int i = 0; i < N_PARTICLES; i++)
            insert(i);

        compute_LJ();
        step();

        if (step_i % 1000 == 0) {
            write_frame(fp, step_i);
        }
    }

    fclose(fp);
    return 0;
}
