#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
using namespace std;

const int L = 100;                   // No. of lattice sites
double hop_rate = 1.0;               // Hop rate
double tumble_rate = 0.01;                  // Tumble rate
double simulation_time = 100.0;     // Total simulation time

mt19937 rng(42);                  // Random number generator with fixed seed
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

int main() {
    Particle p1, p2;

    p1.position = 40;
    p1.orientation = 1;  // Right

    p2.position = 90;
    p2.orientation = -1; // Left

    double current_time = 0.0;

    ofstream traj_file("C:/Users/Maitreya/mips_simulation/week2/task1/trajectory.dat");

    while (current_time < simulation_time) {
        vector<double> rates;
        vector<int> events;

        int x1_next = periodic(p1.position + p1.orientation);
        if (x1_next != p2.position) {
            rates.push_back(hop_rate);
            events.push_back(1); // Particle 1 hops
        }

        int x2_next = periodic(p2.position + p2.orientation);
        if (x2_next != p1.position) {
            rates.push_back(hop_rate);
            events.push_back(2); // Particle 2 hops
        }

        rates.push_back(tumble_rate);
        events.push_back(3); // Particle 1 tumbles

        rates.push_back(tumble_rate);
        events.push_back(4); // Particle 2 tumbles

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
                if (event == 1) {
                    p1.position = periodic(p1.position + p1.orientation);
                } else if (event == 2) {
                    p2.position = periodic(p2.position + p2.orientation);
                } else if (event == 3) {
                    p1.orientation *= -1;
                } else if (event == 4) {
                    p2.orientation *= -1;
                }
                break;
            }
        }
        traj_file << current_time << " " << p1.position << " " << p2.position << endl;
    }
    traj_file.close();
}

