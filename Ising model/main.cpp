#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

// Function prototypes
void display(const vector<vector<int>> &);
void initialize_lattice(vector<vector<int>> &, int);
double energy_diff_single_site(const vector<vector<int>> &, int, int, int, double, double);
void site_flip(vector<vector<int>> &, int, int, double, double);
void monte_carlo_move(vector<vector<int>> &, int, double, double, double);
double system_hamiltonian(const vector<vector<int>> &, int, double, double);
double system_magnetization(const vector<vector<int>> &, int);
void display_vector(const vector<double>);


// Prints 2D spin lattice, showing each site (+1 or -1) in matrix form
void display(const vector<vector<int>> &v){
    cout << showpos;
    for (auto i : v){
        for (auto j : i){
            cout <<  j << " ";
        }
        cout << endl;
    }
    cout << noshowpos << endl;
}

//Creates a N x N lattice of random spin values of +1 or -1
void initialize_lattice(vector<vector<int>> &spins, int N){
    vector<int> temp {};
    int rand_spin_state {};
    for (int i = 0; i < N; i++){
        temp.clear();
        for (int j = 0; j < N; j++){
            rand_spin_state = 2*(rand() % 2) -1;
            temp.push_back(rand_spin_state);
        }
        spins.push_back(temp);
    }
}
// Computes The energy difference of initial energy and energy after a flip at site (a,b) happens
double energy_diff_single_site(const vector<vector<int>> &spins, int a, int b, int N, double J, double h){
    int s {};
    int nn_sum {0};
    s = spins.at(a).at(b);
    nn_sum = spins.at((a+1) % N).at(b) + spins.at(a).at((b+1) % N)
                + spins.at((a-1 + N) % N).at(b) + spins.at(a).at((b-1 +N) % N);
    
    double interaction_energy = 2*J*s*nn_sum;
    double field_energy = 2*h*s;

    return interaction_energy + field_energy;
}

//Implements the metropolis acceptance rule
void site_flip(vector<vector<int>> &spins, int a, int b, double beta, double energy_diff){
    if (energy_diff < 0){
        spins.at(a).at(b) *= -1;
    }else{
        double r = ((double) rand() / (RAND_MAX));
        double prob_of_flip = exp((-1)*beta*(energy_diff));
        if (r < prob_of_flip){
            spins.at(a).at(b) *= -1;
        }
    }
}

// Performs one monte carlo sweep. random site is selected N x N times and
// a flip is attempted based on metropolis algorithm.
void monte_carlo_move(vector<vector<int>> &spins, int N, double J, double beta, double h){
        int a {};
        int b {};
    for (int i = 0; i < N*N; i++){
        a = rand() % N;
        b = rand() % N;
        double energy_diff = energy_diff_single_site(spins, a, b, N, J, h);
        site_flip(spins, a, b, beta, energy_diff);
    }
}

//Calculates the total energy of the lattice
double system_hamiltonian(const vector<vector<int>> &spins, int N, double J, double h){
    double total_energy {0};
    int s {};
    int nn_sum {0};
    for (int i = 0; i <N; i++){
        for (int j = 0; j < N; j++){
            s = spins.at(i).at(j);
            nn_sum = spins.at((i+1) % N).at(j) + spins.at(i).at((j+1) % N)
                + spins.at((i-1 + N) % N).at(j) + spins.at(i).at((j-1 +N) % N);
            //Internal field energy
            total_energy += -0.5*J*s*nn_sum;
            //External field energy
            total_energy += -h*s;
        }
    }
    return total_energy;

}

//Calculates the average magnetization of the lattice
double system_magnetization(const vector<vector<int>> &spins, int N){
    double mag {0};
    for (auto row: spins){
        for (auto spin: row){
            mag += spin;
        }
    }
    return mag /(N*N);
}

// prints all elements of vector
void display_vector(const vector<double> v){
    for (auto i : v){
        cout << i << " ";
    }
    cout << endl;
}


//This is a numerical simulation of Ising model using metropolis algorithm.
int main(){
    srand(time(nullptr));

    // lattice size is N x N
    int N {20};

    // Interaction term
    double J {10};
    // External field
    double h {10};

    vector<vector<int>> spins {};
    vector<double> temperature_values {};
    vector<double> magnetization_values {};
    vector<double> energy_values {};

    initialize_lattice(spins, N);

    double hamiltonian {};
    double magnetization {};
    double beta {};

    for (int i = 0; i < 20; i++){
        // Inverse temperature beta = 1/T
        beta = 0.01 + i*0.002;

        // for each beta performs 1000 monte carlo sweeps to thermalize the system
        for (int j = 0; j < 1000; j++){
            monte_carlo_move(spins, N, J, beta, h);
            hamiltonian = system_hamiltonian(spins, N, J, h);
            magnetization = system_magnetization(spins, N);
        }
        // Stores the magnetization and energy of the last sweep
        temperature_values.push_back(beta);
        magnetization_values.push_back(magnetization);
        energy_values.push_back(hamiltonian);

    }


    cout << "Temperature: " << endl;
    display_vector(temperature_values);
    cout << "\n Magnetization: " << endl;
    display_vector(magnetization_values);
    cout << "\n Energy: " << endl;
    display_vector(energy_values);
    cout << endl;
    //display_vector(energy_values);

    display(spins);

    return 0;
}