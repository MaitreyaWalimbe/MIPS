#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>             // For file output if needed
using namespace std;

int N = 500;                   // Number of particles
double phi = 0.4;              // Packing fraction
double v0 = 1.0;               // Self-propulsion speed
double L = sqrt(N / phi);      // Box size
double Dr = 0.1;               // Rotational diffusion coefficient
double dt = 0.01;              // Time step
int steps = 10000;             // Number of simulation steps
int output_interval = 100;     // Output interval

const double PI = acos(-1.0);

mt19937 rng(42);               // Random number generator with fixed seed
normal_distribution<double> gauss(0.0, 1.0);        // Defines the distribution only
uniform_real_distribution<double> uni(0.0, 1.0);    // Numbers are not generated yet

struct Particle {
    double x, y;       // Position
    double x_actual, y_actual; // Actual position without periodic boundaries
    double theta;      // Orientation
};

double periodic(double coord) {                // Periodic boundary condition
    if (coord < 0) return coord + L;
    if (coord >= L) return coord - L;
    return coord;
}

int main() {
    vector<Particle> particles(N);
    vector<double> x_actual(N);                // To store actual positions
    vector<double> y_actual(N);

    for (int i = 0; i < N; ++i) {              // Initialization of particles
        particles[i].x = uni(rng) * L;
        particles[i].y = uni(rng) * L;
        particles[i].x_actual = particles[i].x;
        particles[i].y_actual = particles[i].y;
        particles[i].theta = uni(rng) * 2 * PI;

        x_actual[i] = particles[i].x_actual;
        y_actual[i] = particles[i].y_actual;
    }

    ofstream msd_file("masd.dat");             // File to store MSD data

    for (int t = 0; t < steps; t++) {
        for (int i = 0; i < N; ++i) {
            particles[i].theta += sqrt(2 * Dr * dt) * gauss(rng);  // Update orientation

            double delta_x = v0 * cos(particles[i].theta) * dt;    // Update position x
            double delta_y = v0 * sin(particles[i].theta) * dt;    // Update position y

            particles[i].x_actual += delta_x;  // Update actual positions
            particles[i].y_actual += delta_y;

            particles[i].x += delta_x;
            particles[i].y += delta_y;

            particles[i].x = periodic(particles[i].x);    // Apply periodic boundary condition
            particles[i].y = periodic(particles[i].y);
        }

        if (t % output_interval == 0) {
            double msd = 0.0;
            for (int i = 0; i < N; ++i) {
                double dx = particles[i].x_actual - x_actual[i];
                double dy = particles[i].y_actual - y_actual[i];
                msd += dx * dx + dy * dy;
            }
            msd /= N;
            msd_file << t * dt << " " << msd << endl;    // Output time and MSD
        }
    }
    msd_file.close();
}