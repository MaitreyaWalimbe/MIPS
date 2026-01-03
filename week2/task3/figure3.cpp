#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
using namespace std;

const int L = 300;                   // No. of lattice sites
double hop_rate = 1.0;               // Hop rate
double tumble_rate = 0.01;                  // Tumble rate
double simulation_time = 100.0;     // Total simulation time
double N = 60;

mt19937 rng(21);                  // Random number generator with fixed seed
uniform_real_distribution<double> uniform_dist(0.0, 1.0);

struct Particle {
    int position;
    int orientation;
};

int periodic(int x) {
    while (x < 0) x += L;
    while (x >= L) x -= L;
    return x;
}

bool is_occupied(int pos, const vector<Particle>& particles) {
    for (const auto& p : particles) {
        if (p.position == pos) return true;
    }
    return false;
}

int main() {
    {
    vector<Particle> Particles(N);
    for (int i = 0; i < N; ++i) {
        Particles[i].position = i * (L / N);
        Particles[i].orientation = (uniform_dist(rng) < 0.5) ? 1 : -1;
    }

    double current_time = 0.0;

    ofstream traj_file("C:/Users/Maitreya/mips_simulation/week2/task3/trajectory_rat.dat");

    while (current_time < simulation_time) {
        vector<double> rates;
        vector<int> events;

        for (int i = 0; i < N; ++i) {
            int x_next = periodic(Particles[i].position + Particles[i].orientation);
            if (!is_occupied(x_next, Particles)) {
                rates.push_back(hop_rate);
                events.push_back(i); // Particle i hops
            }

            rates.push_back(tumble_rate);
            events.push_back(N + i); // Particle i tumbles
        }

        double total_rate = 0.0;
        for (double r : rates) total_rate += r;

        double dt = -log(uniform_dist(rng)) / total_rate;
        current_time += dt;

        double r = uniform_dist(rng) * total_rate;
        double cumulative_rate = 0.0;

        for (int i = 0; i < rates.size(); ++i) {
            cumulative_rate += rates[i];
            if (r < cumulative_rate) {
                int event = events[i];
                if (event < N) {
                    // Hop event
                    int idx = event;
                    Particles[idx].position = periodic(Particles[idx].position + Particles[idx].orientation);
                } else {
                    // Tumble event
                    int idx = event - N;
                    Particles[idx].orientation *= -1;
                }
                break;
            }
        }
        traj_file << current_time;
        for (const auto& p : Particles) {
            traj_file << " " << p.position;
        }
        traj_file << "\n";
    }
    traj_file.close();
    }

    {
    vector<Particle> Particles(N);
    for (int i = 0; i < N; ++i) {
        Particles[i].position = i * (L / N);
    }
    double current_time = 0.0;
    ofstream traj_file("C:/Users/Maitreya/mips_simulation/week2/task3/trajectory_srw.dat");

    while (current_time < simulation_time) {
        vector<double> rates;
        vector<int> events;

        for (int i = 0; i < N; ++i) {
            int x_next_right = periodic(Particles[i].position + 1);
            if (!is_occupied(x_next_right, Particles)) {
                rates.push_back(hop_rate/2.0);
                events.push_back(2 * i); // Particle i hops right
            }

            int x_next_left = periodic(Particles[i].position - 1);
            if (!is_occupied(x_next_left, Particles)) {
                rates.push_back(hop_rate/2.0);
                events.push_back(2 * i + 1); // Particle i hops left
            }
        }

        double total_rate = 0.0;
        for (double r : rates) total_rate += r;

        double dt = -log(uniform_dist(rng)) / total_rate;
        current_time += dt;

        double r = uniform_dist(rng) * total_rate;
        double cumulative_rate = 0.0;

        for (int i = 0; i < rates.size(); ++i) {
            cumulative_rate += rates[i];
            if (r < cumulative_rate) {
                int event = events[i];
                int idx = event / 2;
                if (event % 2 == 0) {
                    // Hop right
                    Particles[idx].position = periodic(Particles[idx].position + 1);
                } else {
                    // Hop left
                    Particles[idx].position = periodic(Particles[idx].position - 1);
                }
                break;
            }
        }
        traj_file << current_time;
        for (const auto& p : Particles) {
            traj_file << " " << p.position;
        }
        traj_file << "\n";
    }
    traj_file.close();
    }
}

