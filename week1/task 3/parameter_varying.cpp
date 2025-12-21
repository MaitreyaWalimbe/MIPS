#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>             // For file output if needed
#include <sstream>             // For string stream if needed
using namespace std;

int N = 500;                   // Number of particles
double a = 1.0;                // Particle radius
double Dr = 0.1;               // Rotational diffusion coefficient
double dt = 0.01;              // Time step
int steps = 10000;             // Number of simulation steps
int output_interval = 100;     // Output interval

double r_cutoff = 2.0 * a;        // Cutoff distance for repulsion
double r_skin = 0.2 * a;          // Skin distance for neighbor list
double r_list = r_cutoff + r_skin; // Neighbor list cutoff

const double PI = acos(-1.0);

mt19937 rng(42);               // Random number generator with fixed seed
normal_distribution<double> gauss(0.0, 1.0);        // Defines the distribution only
uniform_real_distribution<double> uni(0.0, 1.0);    // Numbers are not generated yet

struct Particle {
    double x, y;       // Position
    double x_reference, y_reference; // Neighbor list reference position for rebuild
    double x_actual, y_actual; // Actual position without periodic boundaries
    double theta;      // Orientation
};

double L;                      // Box size

double periodic(double coord) {                // Periodic boundary condition
    if (coord < 0) return coord + L;
    if (coord >= L) return coord - L;
    return coord;
}

double minimum_image(double dx) {              // Minimum image convention
    if (dx > L / 2) return dx - L;
    if (dx < -L / 2) return dx + L;
    return dx;
}

vector<vector<int>> neighbor_list;            // Neighbor list

void build_neighbor_list(const vector<Particle>& particles) {
    neighbor_list.clear();
    neighbor_list.resize(N);

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            if (i==j) continue;
            double dx = minimum_image(particles[i].x-particles[j].x);
            double dy = minimum_image(particles[i].y-particles[j].y);
            double r = sqrt(dx*dx + dy*dy);
            if (r < r_list) {
                neighbor_list[i].push_back(j);
                neighbor_list[j].push_back(i);
            }
        }
    }
}

int main() {
    vector<double> phi_list = {0.2, 0.4, 0.6}; // Different packing fractions to simulate
    vector<double> pe_list = {5.0, 20.0, 50.0, 100.0}; // Different Péclet numbers to simulate

    for (double phi : phi_list) {
        L = sqrt(N / phi);

        for (double Pe : pe_list) {
            double v0 = Pe * Dr * a;   // Self-propulsion speed

            cout << "Running simulation for phi=" << phi << ", Pe=" << Pe << endl;

            vector<Particle> particles(N);
            vector<double> x_actual(N);                // To store actual positions
            vector<double> y_actual(N);

            for (int i = 0; i < N; ++i) {              // Initialization of particles
                particles[i].x = uni(rng) * L;
                particles[i].y = uni(rng) * L;
                particles[i].x_reference = particles[i].x;
                particles[i].y_reference = particles[i].y;
                particles[i].x_actual = particles[i].x;
                particles[i].y_actual = particles[i].y;
                particles[i].theta = uni(rng) * 2 * PI;

                x_actual[i] = particles[i].x_actual;
                y_actual[i] = particles[i].y_actual;
            }

            build_neighbor_list(particles);          // Initial neighbor list

            stringstream fname;
            fname << "C:/Users/Maitreya/mips_simulation/week1/task3/"<< "msd_phi" << int(phi*10) << "_Pe" << int(Pe) << ".dat";
            cout << "Attempting to write to:\n" << fname.str() << endl;
            ofstream msd_file(fname.str());

            for (int t = 0; t < steps; t++) {
                vector<double> fx(N, 0.0);    // Forces in x
                vector<double> fy(N, 0.0);    // Forces in y

                for (int i = 0; i < N; ++i) {
                    for (int j : neighbor_list[i]) {
                        if (j <= i) continue; // Avoid double counting
                        double dx = minimum_image(particles[i].x - particles[j].x);
                        double dy = minimum_image(particles[i].y - particles[j].y);
                        double r = sqrt(dx * dx + dy * dy);

                        if (r < r_cutoff && r > 1e-12) { // Avoid division by zero
                            double overlap = r_cutoff - r;
                            double force_mag = overlap; // Harmonic repulsion
                            fx[i] += force_mag * (dx / r);
                            fy[i] += force_mag * (dy / r);
                            fx[j] -= force_mag * (dx / r);        // To uphold Newton's third law
                            fy[j] -= force_mag * (dy / r);
                        }
                    }
                
                }

                double max_displacement = 0.0;

                for (int i = 0; i < N; ++i) {
                    particles[i].theta += sqrt(2 * Dr * dt) * gauss(rng);  // Update orientation

                    double delta_x = (v0 * cos(particles[i].theta) + fx[i]) * dt;    // Update position x
                    double delta_y = (v0 * sin(particles[i].theta) + fy[i]) * dt;    // Update position y
                    
                    particles[i].x_actual += delta_x;  // Update actual positions
                    particles[i].y_actual += delta_y;

                    particles[i].x += delta_x;
                    particles[i].y += delta_y;

                    particles[i].x = periodic(particles[i].x);    // Apply periodic boundary condition
                    particles[i].y = periodic(particles[i].y);

                    double drx = minimum_image(particles[i].x - particles[i].x_reference);
                    double dry = minimum_image(particles[i].y - particles[i].y_reference);
                    max_displacement = max(max_displacement, sqrt(drx*drx + dry*dry));
                }

                if (max_displacement > r_skin / 2.0) {
                    for (int i = 0; i < N; ++i) {
                        particles[i].x_reference = particles[i].x;
                        particles[i].y_reference = particles[i].y;
                    }
                    build_neighbor_list(particles);      // Rebuild neighbor list
                }
                if (t % output_interval == 0) {
                    double msd_val = 0.0;
                    for (int i = 0; i < N; ++i) {
                        double dx = particles[i].x_actual - x_actual[i];
                        double dy = particles[i].y_actual - y_actual[i];
                        msd_val += dx * dx + dy * dy;
                    }
                    msd_val /= N;
                    msd_file << t * dt << " " << msd_val << endl;    // Output time and MSD
                }
            }
            msd_file.close();
        }
    }
    cout << "All simulations completed." << endl;
}