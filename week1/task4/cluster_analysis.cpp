#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>             // For file output if needed
#include <sstream>            // For string stream operations
using namespace std;

const double PI = acos(-1.0);
int N_ens = 5;               // Number of ensembles


int N = 500;                   // Number of particles
double a = 1.0;                // Particle diameter
double phi = 0.4;              // Packing fraction
double Dr = 0.001;               // Rotational diffusion coefficient
double Dt = 0.0;               // Translational diffusion coefficient
double dt = 0.005;              // Time step
int steps = 20000;             // Number of simulation steps
int equilibration_steps = 10000; // Equilibration steps
int cluster_analysis_interval = 2000; // Cluster analysis interval

double L = sqrt(N * PI * a * a / (phi * 4));      // Box size

double r_cutoff = a;        // Cutoff distance for repulsion
double r_skin = 0.1 * a;          // Skin distance for neighbor list   // 0.1a
double r_list = r_cutoff + r_skin; // Neighbor list cutoff

double r_cluster = 1.2 * a;    // 1.2 * diameter for cluster analysis

double k_rep = 100.0;          // Repulsion strength (changed from 100 to better observe mips curve)
double mu = 1.0;               // Mobility

mt19937 rng(42);               // Random number generator with fixed seed
normal_distribution<double> gauss(0.0, 1.0);        // Defines the distribution only
uniform_real_distribution<double> uni(0.0, 1.0);    // Numbers are not generated yet

struct Particle {
    double x, y;       // Position
    double x_reference, y_reference; // Neighbor list reference position for rebuild
    double theta;      // Orientation
};

double periodic(double coord) {                // Periodic boundary condition
    while (coord < 0) coord += L;
    while (coord >= L) coord -= L;
    return coord;
}

double minimum_image(double dx) {              // Minimum image convention
    if (dx > L / 2) return dx - L;
    if (dx < -L / 2) return dx + L;
    return dx;
}

vector<vector<int>> neighbor_list;            // Neighbor list
vector<vector<int>> cluster_neighbor_list;    // Cluster neighbor list

void build_neighbor_list(const vector<Particle>& particles) {
    neighbor_list.clear();
    neighbor_list.resize(N);

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
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

void build_cluster_neighbor_list(const vector<Particle>& particles) {    // Separate neighbor list because of different cutoff
    cluster_neighbor_list.clear();
    cluster_neighbor_list.resize(N);


    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            if (i==j) continue;
            double dx = minimum_image(particles[i].x-particles[j].x);
            double dy = minimum_image(particles[i].y-particles[j].y);
            double r2 = dx*dx + dy*dy;
            if (r2 < r_cluster * r_cluster) {
                cluster_neighbor_list[i].push_back(j);
                cluster_neighbor_list[j].push_back(i);
            }
        }
    }
}

int largest_cluster_size(const vector<Particle>& particles) {
    vector<bool> visited(N, false);
    int max_cluster_size = 0;

    for (int i = 0; i < N; ++i) {
        if (visited[i]) continue;

        int cluster_size = 0;
        vector<int> stack = {i};
        visited[i] = true;

        while (!stack.empty()) {
            int curr = stack.back();
            stack.pop_back();
            cluster_size++;

            for (int j : cluster_neighbor_list[curr]) {
                if (visited[j]) continue;
                visited[j] = true;
                stack.push_back(j);
            }
        }
        max_cluster_size = max(max_cluster_size, cluster_size);
    }
    return max_cluster_size;
}

int main() {
    vector<double> pe_list = {2, 5, 8,
                              10, 12, 15, 18, 22, 26, 30, 35, 40, 45, 50 , 55, 60,
                              65, 75, 90, 100};

    ofstream fout("C:/Users/Maitreya/mips_simulation/week1/task4/fmax_vs_Pe.dat");
    fout << "# Pe <f_max>\n";
    
    for (double Pe : pe_list) {
        double fmax_ensemble_sum = 0.0;

        for (int ens = 0; ens < N_ens; ++ens) {
            mt19937 rng_local(42 + ens * 1000 + int(Pe * 10)); // Different seed for each ensemble
            
            double v0 = Pe * Dr * a;   // Self-propulsion speed

            cout << "Running simulation for Pe=" << Pe << endl;

            vector<Particle> particles(N);

            for (int i = 0; i < N; ) {    // Initialization with no overlap
                double x = uni(rng_local) * L;
                double y = uni(rng_local) * L;

                bool overlap = false;
                for (int j = 0; j < i; ++j) {
                    double dx = minimum_image(x - particles[j].x);
                    double dy = minimum_image(y - particles[j].y);
                    if (dx*dx + dy*dy < a*a) {
                        overlap = true;
                        break;
                    }
                }

                if (!overlap) {
                    particles[i].x = x;
                    particles[i].y = y;
                    particles[i].x_reference = x;
                    particles[i].y_reference = y;
                    particles[i].theta = uni(rng_local) * 2 * PI;
                    i++;
                }

}

            build_neighbor_list(particles);          // Initial neighbor list
            build_cluster_neighbor_list(particles);  // Initial cluster neighbor list

            double f_max_accum = 0.0;
            int f_max_count = 0;

            vector<double> fx(N, 0.0);    // Forces in x    // Defined outside loop to reduce time complexity
            vector<double> fy(N, 0.0);    // Forces in y

            for (int t = 0; t < steps; t++) {
                fill(fx.begin(), fx.end(), 0.0);   // Reset forces
                fill(fy.begin(), fy.end(), 0.0);

                for (int i = 0; i < N; ++i) {
                    for (int j : neighbor_list[i]) {
                        if (j <= i) continue; // Avoid double counting
                        double dx = minimum_image(particles[i].x - particles[j].x);
                        double dy = minimum_image(particles[i].y - particles[j].y);
                        double r2 = dx * dx + dy * dy;

                        if (r2 < r_cutoff * r_cutoff && r2 > 1e-12) { // Avoid division by zero
                            double r = sqrt(r2);
                            double overlap = r_cutoff - r;
                            double force_mag = k_rep * overlap; // Harmonic repulsion
                            fx[i] += force_mag * (dx / r);
                            fy[i] += force_mag * (dy / r);
                            fx[j] -= force_mag * (dx / r);        // To uphold Newton's third law
                            fy[j] -= force_mag * (dy / r);
                        }
                    }
                }

                double max_displacement = 0.0;

                for (int i = 0; i < N; ++i) {
                    particles[i].theta += sqrt(2 * Dr * dt) * gauss(rng_local);  // Update orientation

                    double delta_x = (v0 * cos(particles[i].theta) + mu * fx[i]) * dt;    // Update position x
                    double delta_y = (v0 * sin(particles[i].theta) + mu * fy[i]) * dt;    // Update position y
                        
                    particles[i].x += delta_x + sqrt(2 * Dt * dt) * gauss(rng_local); // Add translational diffusion
                    particles[i].y += delta_y + sqrt(2 * Dt * dt) * gauss(rng_local);

                    particles[i].x = periodic(particles[i].x);    // Apply periodic boundary condition
                    particles[i].y = periodic(particles[i].y);

                    double drx = minimum_image(particles[i].x - particles[i].x_reference);
                    double dry = minimum_image(particles[i].y - particles[i].y_reference);
                    max_displacement = max(max_displacement, sqrt(drx*drx + dry*dry));
                }

                if (max_displacement > r_skin) {
                    for (int i = 0; i < N; ++i) {
                        particles[i].x_reference = particles[i].x;
                        particles[i].y_reference = particles[i].y;
                    }
                    build_neighbor_list(particles);      // Rebuild neighbor list
                }

                if (t == equilibration_steps) {
                    for (int i = 0; i < N; ++i) {
                        particles[i].x_reference = particles[i].x;
                        particles[i].y_reference = particles[i].y;
                    }
                }

                if (t > equilibration_steps && t % cluster_analysis_interval == 0) {
                    build_cluster_neighbor_list(particles);  // Update cluster neighbor list
                    f_max_accum += double(largest_cluster_size(particles)) / N;
                    f_max_count++;
                }

            }
            double f_max_avg = f_max_accum / f_max_count;
            fmax_ensemble_sum += f_max_avg;
        }
        double fmax_final = fmax_ensemble_sum / N_ens;
        fout << Pe << " " << fmax_final << "\n";
        cout << "Pe=" << Pe << " Ensemble averaged <f_max>=" << fmax_final << endl;
    }
    fout.close();
}