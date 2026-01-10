#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <sstream>
using namespace std;

const double PI = acos(-1.0);

int N = 1000;     // Number of boids
double interaction_radius = 1.0;        // Interaction radius
double v0 = 0.03;        // Preferred speed
double phi = 1.0;       // Density
double L = sqrt(N / phi);             // Size of cell
double dt = 1.0;          // Time step
int total_steps = 10000; // Total simulation steps

double r_cutoff = interaction_radius;
double r_skin = 1.0 * r_cutoff;
double r_list = r_cutoff + r_skin;

struct Boid {
    double x, y;      // Position
    double x_ref, y_ref; // Reference position for neighbor list
    double theta;    // Direction angle
};

double periodic(double x) {
    while (x < 0) x = x + L;
    while (x >= L) x = x - L;
    return x;
}

double minimum_image(double dx) {
    if (dx >  L/2) return dx - L;
    if (dx < -L/2) return dx + L;
    return dx;
}

vector<vector<int>> neighbor_list;

void build_neighbor_list(const vector<Boid>& boids) {
    neighbor_list.clear();
    neighbor_list.resize(N);

    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double dx = minimum_image(boids[i].x - boids[j].x);
            double dy = minimum_image(boids[i].y - boids[j].y);
            double r2 = dx*dx + dy*dy;

            if (r2 < r_list * r_list) {
                neighbor_list[i].push_back(j);
                neighbor_list[j].push_back(i);
            }
        }
    }
}

mt19937 rng(21);                  // Random number generator with fixed seed
uniform_real_distribution<double> uniform_dist(0.0, 1.0);

int main() {
    vector<double> noise_levels = {0.2, 0.4, 0.6, 0.8}; // Different noise levels

    for (double eta : noise_levels) {
        cout << "Simulating for noise level: " << eta << endl;

        vector<Boid> boids(N);       // Initialize boids
        for (int i = 0; i < N; i++) {
            boids[i].x = uniform_dist(rng) * L;
            boids[i].y = uniform_dist(rng) * L;
            boids[i].theta = uniform_dist(rng) * 2 * PI;

            boids[i].x_ref = boids[i].x;
            boids[i].y_ref = boids[i].y;
        }

        build_neighbor_list(boids);

        for (int t = 0; t < total_steps; t++) {
            vector<double> new_thetas(N);

            for (int i = 0; i < N; i++) {
                double sx = cos(boids[i].theta), sy = sin(boids[i].theta);

                for (int j : neighbor_list[i]) {
                    double dx = minimum_image(boids[j].x - boids[i].x);
                    double dy = minimum_image(boids[j].y - boids[i].y);

                    if (dx*dx + dy*dy < interaction_radius*interaction_radius) {
                        sx += cos(boids[j].theta);
                        sy += sin(boids[j].theta);
                    }
                }
                new_thetas[i] = atan2(sy, sx) + (uniform_dist(rng) - 0.5) * eta * 2.0 * PI;
            }
            double max_displacement2 = 0.0;

            for (int i = 0; i < N; i++) {
                boids[i].theta = new_thetas[i];
                boids[i].x = periodic(boids[i].x + v0 * cos(boids[i].theta) * dt);
                boids[i].y = periodic(boids[i].y + v0 * sin(boids[i].theta) * dt);

                double dx_ref = minimum_image(boids[i].x - boids[i].x_ref);
                double dy_ref = minimum_image(boids[i].y - boids[i].y_ref);

                max_displacement2 = max(max_displacement2, dx_ref*dx_ref + dy_ref*dy_ref);
            }
            if (max_displacement2 > r_skin * r_skin) {
                for (int i = 0; i < N; i++) {
                    boids[i].x_ref = boids[i].x;
                    boids[i].y_ref = boids[i].y;
                }
                build_neighbor_list(boids);
            }
        }
        stringstream fname;
        fname << "C:/Users/Maitreya/mips_simulation/week3/task1/snapshot_eta_" << int(eta * 10) << ".dat";
        ofstream file(fname.str());
        if (!file.is_open()) {
            cerr << "ERROR: Could not open file:\n" << fname.str() << endl;
            return 1;
        }
        for (int i = 0; i < N; i++) {
            file << boids[i].x << " " << boids[i].y << " " << boids[i].theta << "\n";
        }

        file.close();
    }
}