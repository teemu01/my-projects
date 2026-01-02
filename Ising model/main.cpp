#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>

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
void compute_thermodynamical_quantities(int, double, double, int, int,
                                const vector<vector<int>> &,
                                const vector<double> &,
                                vector<double> &, vector<double> &,
                                vector<double> &, vector<double> &);




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
void display_vector(const vector<double> &v){
    for (auto i : v){
        cout << i << " ";
    }
    cout << endl;
}


void compute_thermodynamical_quantities(int N, double J, double h, int thermalisation_number, int mean_number,
                                vector<vector<int>> &spins,
                                const vector<double> &temperature_values,
                                vector<double> &energy_values,
                                vector<double> &magnetization_values,
                                vector<double> &susceptibility_values,
                                vector<double> &specific_heat_values){



    for (auto temperature: temperature_values){
        double beta = 1/temperature;

        //Thermalisation of the Ising model
        for (int j = 0; j < thermalisation_number; j++){
            monte_carlo_move(spins, N, J, beta, h);
        }

        double sum_magnetization {0.0};
        double sum_magnetization_squared {0.0};
        double sum_energy {0.0};
        double sum_energy_squared {0.0};

        //Loop for thermalised system to calculate the average energy and magnetization
        for (int j = 0; j < mean_number; j++){
            monte_carlo_move(spins, N, J, beta, h);

            double E = system_hamiltonian(spins, N, J, h);
            sum_energy += E;
            sum_energy_squared += E*E;
            double mag = system_magnetization(spins, N);
            sum_magnetization += abs(mag);
            sum_magnetization_squared += mag * mag;
        }

        double energy = sum_energy / mean_number;
        double energy_squared = sum_energy_squared / mean_number;
        double magnetization = sum_magnetization / mean_number;
        double magnetization_squared = sum_magnetization_squared / mean_number;
        double susceptibility = beta*N*N*(magnetization_squared - magnetization*magnetization);
        double specific_heat = beta*beta*N*N*(energy_squared - energy*energy);

        magnetization_values.push_back(magnetization);
        energy_values.push_back(energy);
        susceptibility_values.push_back(susceptibility);
        specific_heat_values.push_back(specific_heat);

    }
}




//This is a numerical simulation of Ising model using metropolis algorithm.
int main(){

    // lattice size is N x N
    int N {30};
    // Interaction term
    double J {1};
    // External field
    double h {0};
    
    double T_initial {1.5};
    int thermalisation_number {5000};
    int mean_number {3000};


    //Initialize lattice with random spin values (+1 or -1) for each site
    vector<vector<int>> spins {};
    srand(time(nullptr));
    initialize_lattice(spins, N);


    //List of different thermodynamic quantities to derive
    vector<double> magnetization_values {};
    vector<double> susceptibility_values {};
    vector<double> specific_heat_values {};
    vector<double> energy_values {};


    vector<double> temperature_values {};
    for (int i = 0; i < 60; i++){
        double temperature = T_initial + 0.02*i;
        temperature_values.push_back(temperature);
    }

    compute_thermodynamical_quantities(N, J, h, thermalisation_number, mean_number,
                                spins, temperature_values, energy_values,
                                magnetization_values, susceptibility_values,
                                specific_heat_values);



    //#####################

    TGraph *gr1 = new TGraph(temperature_values.size(),
                 temperature_values.data(),
                 magnetization_values.data());

    TCanvas *c1 = new TCanvas("c1", "Graph Draw Options", 200, 10, 800, 600);


    gr1->SetTitle("Magnetization vs Temperature;T;M");
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(1.2);
    gr1->SetLineWidth(2);

    //TF1 *fitFcn = new TF1("fitFcn",fitCurieFunction,1.5,2.3,2);

    //gr1->Fit("fitFcn");

    gr1-> Draw("AC*");

    c1->SaveAs("magnetization_vs_temperature.png");

    //#####################

    TGraph *gr2 = new TGraph(temperature_values.size(),
                 temperature_values.data(),
                 energy_values.data());

    TCanvas *c2 = new TCanvas("c2", "Graph Draw Options", 200, 10, 800, 600);

    gr2->SetTitle("Energy vs Temperature;T;E");
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerSize(1.2);
    gr2->SetLineWidth(2);

    gr2-> Draw("AC*");

    c2->SaveAs("energy_vs_temperature.png");

    //###############################

    TGraph *gr3 = new TGraph(temperature_values.size(),
                 temperature_values.data(),
                 susceptibility_values.data());

    TCanvas *c3 = new TCanvas("c3", "Graph Draw Options", 200, 10, 800, 600);

    gr3->SetTitle("Susceptibility vs Temperature;T;X");
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerSize(1.2);
    gr3->SetLineWidth(2);

    gr3-> Draw("AC*");

    c3->SaveAs("susceptibility_vs_temperature.png");

    //#######################################
    TGraph *gr4 = new TGraph(temperature_values.size(),
                 temperature_values.data(),
                 specific_heat_values.data());

    TCanvas *c4 = new TCanvas("c4", "Graph Draw Options", 200, 10, 800, 600);

    gr4->SetTitle("Specific Heat vs Temperature;T;C");
    gr4->SetMarkerStyle(20);
    gr4->SetMarkerSize(1.2);
    gr4->SetLineWidth(2);


    gr4-> Draw("AC*");

    c4->SaveAs("specific_heat_vs_temperature.png");



    return 0;
}
