#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <sstream>
using namespace std;

const double PI = acos(-1.0);
int N_ens = 3;               // Number of ensemble runs

double interaction_radius = 1.0;
double v0 = 0.03;
int L = 20;
double noise_level = 0.6;
double dt = 1.0;
int equilibration_steps = 2000;
int measure_steps = 2000;
int total_steps = equilibration_steps + measure_steps;

struct Boid {
    double x, y;
    double theta;
};

double periodic(double x) {
    while (x < 0) x += L;
    while (x >= L) x -= L;
    return x;
}

double minimum_image(double dx) {
    if (dx >  L/2) return dx - L;
    if (dx < -L/2) return dx + L;
    return dx;
}

/* ================= CELL LIST STRUCTURES ================= */

double cell_size = interaction_radius;
int ncell;
vector<vector<int>> cells;

inline int cell_index(int ix, int iy) {
    ix = (ix + ncell) % ncell;
    iy = (iy + ncell) % ncell;
    return ix + iy * ncell;
}

void build_cell_list(const vector<Boid>& boids, int N) {
    for (auto &c : cells) c.clear();

    for (int i = 0; i < N; i++) {
        int cx = int(boids[i].x / cell_size);
        int cy = int(boids[i].y / cell_size);
        cells[cell_index(cx, cy)].push_back(i);
    }
}

/* ======================================================== */

vector<double> linspace(double start, double end, int num) {
    vector<double> v;
    if (num == 1) return {start};
    double dx = (end - start) / (num - 1);
    for (int i = 0; i < num; i++) v.push_back(start + i * dx);
    return v;
}

vector<double> density_values = linspace(1.0, 12.0, 10);

int main() {
    stringstream fname;
    fname << "C:/Users/Maitreya/mips_simulation/week3/task3/phi_vs_rho.dat";
    ofstream fout(fname.str());

    for (double rho : density_values) {
        int N = int(rho * L * L);

        ncell = int(L / cell_size);
        cells.assign(ncell * ncell, vector<int>());

        double phi_ens_accum = 0.0;

        for (int ens = 0; ens < N_ens; ens++) {
            mt19937 rng_local(42 + ens + int(100 * rho));
            uniform_real_distribution<double> uniform_dist(0.0, 1.0);

            vector<Boid> boids(N);
            for (int i = 0; i < N; i++) {
                boids[i].x = uniform_dist(rng_local) * L;
                boids[i].y = uniform_dist(rng_local) * L;
                boids[i].theta = uniform_dist(rng_local) * 2 * PI;
            }

            build_cell_list(boids, N);

            double phi_accum = 0.0;
            int phi_count = 0;

            for (int t = 0; t < total_steps; t++) {
                vector<double> new_thetas(N);

                for (int i = 0; i < N; i++) {
                    double sx = cos(boids[i].theta);
                    double sy = sin(boids[i].theta);

                    int cx = int(boids[i].x / cell_size);
                    int cy = int(boids[i].y / cell_size);

                    for (int dx = -1; dx <= 1; dx++) {
                        for (int dy = -1; dy <= 1; dy++) {
                            int c = cell_index(cx + dx, cy + dy);
                            for (int j : cells[c]) {
                                double dxp = minimum_image(boids[j].x - boids[i].x);
                                double dyp = minimum_image(boids[j].y - boids[i].y);
                                if (dxp*dxp + dyp*dyp < interaction_radius*interaction_radius) {
                                    sx += cos(boids[j].theta);
                                    sy += sin(boids[j].theta);
                                }
                            }
                        }
                    }

                    new_thetas[i] =
                        atan2(sy, sx)
                      + (uniform_dist(rng_local) - 0.5) * noise_level * 2.0 * PI;
                }

                for (int i = 0; i < N; i++) {
                    boids[i].theta = new_thetas[i];
                    boids[i].x = periodic(boids[i].x + v0 * cos(boids[i].theta) * dt);
                    boids[i].y = periodic(boids[i].y + v0 * sin(boids[i].theta) * dt);
                }

                build_cell_list(boids, N);

                if (t >= equilibration_steps) {
                    double vx_sum = 0.0, vy_sum = 0.0;
                    for (int i = 0; i < N; i++) {
                        vx_sum += cos(boids[i].theta);
                        vy_sum += sin(boids[i].theta);
                    }
                    phi_accum += sqrt(vx_sum*vx_sum + vy_sum*vy_sum) / N;
                    phi_count++;
                }
            }

            phi_ens_accum += phi_accum / phi_count;
        }

        fout << rho << " " << (phi_ens_accum / N_ens) << endl;
    }

    fout.close();
}
