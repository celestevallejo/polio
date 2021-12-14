//
//  PDE_Simulator.cpp
//  
//
//  Created by Celeste on 1/3/18.
//
//

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <assert.h>
#include <vector>
#include <array>
#include <deque>
#include <fstream>
//#include <gsl/gsl_odeiv2.h>
//#include <gsl/gsl_math.h>
//#include "SIR_H.h"

using namespace std;

void usage() {
    
    cerr << "\n\tUsage: ./pde <max_simulationTimeSteps> <waningTimeSteps> <equilibrium_epsilon> <equilibrium_window>\n\n";
    exit(-1);
    
}

int main(int argc, char** argv){
    
    if (argc!=5) usage();
    /*ofstream myfile;
     ofstream myfile1;
     ofstream myfile2;
     ofstream myfile3;
     myfile.open("/Users/Celeste/Desktop/pde_test_Ir.csv");
     myfile1.open("/Users/Celeste/Desktop/pde_test_I1.csv");
     myfile2.open("/Users/Celeste/Desktop/pde_test_S.csv");
     myfile3.open("/Users/Celeste/Desktop/pde_test_R.csv");*/
    //within-host parameters
    const double mu0 = 0.1841; //pathogen growth rate
    const double mu = 0.2; //antibody growth rate
    const double b0 = 1; //initial pathogen concentration
    const double y0 = 1; //initial antibody concentration
    const double y1 = pow(10,2);//peak antibody concentration
    const double t1 = (1.0/mu)*log(y1/y0); //duration of infection (days)
    const double c = b0*((mu-mu0)/(y0*(exp((mu-mu0)*t1)-1))); //pathogen clearance rate
    const double r =1.7; //immunity shape parameter
    const double nu = 2; //immunity waning rate
    
    //between-host parameters
    const double K1 = (100/(double)365); //I1 daily contact rate
    const double K2 = (50/(double)365); //Ir daily contact rate
    const double alpha1 = (13/(double)365); //I1 daily recovery parameter
    const double alpha2 = (31/(double)365); //Ir daily recovery parameter
    const double omega = (0.2/(double)365); //daily waning parameter
    const double C = 0.0001; //saturation constant
    const double C1 = 0.0001; //constant for recovery functions
    const double C2 = 0.0001; //constant for waning function
    const double delta = (0.02/(double)365); //daily birth/death rate
    //const double totalPop = 101; //total population size
    
    //time parameters
    const double dt = 0.01; //time step in y direction (simulation time direction)
    const double dtau = 0.01; //time step in x direction (time since infection direction)
    const double dwane = 0.01; //time step in x direction (time since recovery direction)
    const unsigned int simulationTimeSteps = atoi(argv[1])/dt; //delta t * simulationTimeSteps = final time of interest T (days)
    const unsigned int recoveryTimeSteps = t1/dtau; //delta tau * recoveryTimeSteps = final infectious time of interest
    const unsigned int waningTimeSteps = atoi(argv[2])/dt; //delta t * waningTimeSteps = final waning time of interest
    const unsigned int withinHostTimeSteps = recoveryTimeSteps + waningTimeSteps; //delta t * withinHostTimeSteps = final time since infection tau (days)
    
    //convergence parameters
    const double epsilon = atof(argv[3]);
    const unsigned int eq_interval = atoi(argv[4])/dt;
    pair<unsigned int, double> obs_min;
    pair<unsigned int, double> obs_max;
    deque<double> obs_deque;
    //const double T = dt*N;
    //const double tau = dt*M;
    
    //vectors for linking functions
    vector<double>  beta1(recoveryTimeSteps,0.0);
    vector<double>  beta2(recoveryTimeSteps,0.0);
    vector<double> gamma1(recoveryTimeSteps,0.0);
    vector<double> gamma2(recoveryTimeSteps,0.0);
    vector<double>    rho(waningTimeSteps,0.0);
    
    int rhoCounter = 0;
    //vectors for within-host solution
    //Construct linking functions
    for (unsigned int j = 0; j < withinHostTimeSteps; ++j) {
        const double t = dt*j;
        double pathogen;
        double antibody;
        if(t<=t1){
            pathogen = b0*exp(mu0*t) - (c*y0*(exp(mu*t)-exp(mu0*t))/(mu-mu0));
            antibody = y0*exp(mu*t);
        }
        else{
            pathogen = 0;
            antibody = y1*pow((1+(r-1)*pow(y1,r-1)*nu*(t-t1)),-(1/(r-1)));
        }
        
        if(j<=t1/dt){//integrals of these quantities only defined up until recovery time
            const double beta1Link = (K1*pathogen)/(pathogen+C);
            const double beta2Link = (K2*pathogen)/(pathogen+C);
            const double gamma1Link = (alpha1*antibody)/(pathogen+C1);
            const double gamma2Link = (alpha2*antibody)/(pathogen+C1);
            
            gamma1[j] = gamma1Link;
            gamma2[j] = gamma2Link;
            beta1[j]  = beta1Link;
            beta2[j]  = beta2Link;
        }
        
        else{//only keeps track of waning rate during waning period
            rhoCounter = j - recoveryTimeSteps-1;
            const double rhoLink = omega/(antibody+C2);
            rho[rhoCounter] = rhoLink;
        }
        
    }
    
    //vectors for between-host solution
    cerr << "sT: " << simulationTimeSteps << "\n";
    vector<double> S(simulationTimeSteps,0.0);
    cerr << "b\n";
    //initialize between-host compartments
    S[0] = 1000;
    vector<vector<double> > I1(2,vector<double>(recoveryTimeSteps,0.0));
    vector<vector<double> > R(2,vector<double>(waningTimeSteps,0.0));
    vector<vector<double> > Ir(2,vector<double>(recoveryTimeSteps,0.0));
    vector<double> symptomaticIncidence (simulationTimeSteps,0.0);
    vector<double> asymptomaticIncidence (simulationTimeSteps,0.0);
    I1[0][0] = 100;
    
    double I1pop = 0;
    double Rpop  = 0;
    double Irpop = 0;
    for(unsigned int j = 0; j < recoveryTimeSteps; j++){
        I1pop += I1[0][j];
        Irpop += Ir[0][j];
    }
    for(unsigned int j = 0; j< waningTimeSteps; j++){
        Rpop +=   R[0][j];
    }
    double totalPop = S[0] + I1pop + Rpop + Irpop;
    symptomaticIncidence[0] = S[0]*(beta1[0]*I1[0][0]+beta2[0]*Ir[0][0])/totalPop;
    obs_min = {0,symptomaticIncidence[0]};
    obs_max = {0,symptomaticIncidence[0]};
    
    //Finite Difference Method
    // use backward Euler difference quotient to approximate time derivatives
    // approximate integrals using right end point rule
    
    for(unsigned int k=1; k<simulationTimeSteps; k++){//looping through time (rows)
        if (k % 1000 == 0) cerr << k << " " << obs_max.second - obs_min.second << " | " <<  obs_max.first << "," << obs_max.second << endl;
        //calculate the total population
        double I1pop = 0;
        double Rpop  = 0;
        double Irpop = 0;
        for(unsigned int j = 0; j < recoveryTimeSteps; j++){
            I1pop += I1[0][j];
            Irpop += Ir[0][j];
        }
        for(unsigned int j = 0; j< waningTimeSteps; j++){
            Rpop +=   R[0][j];
        }
        double totalPop = S[k-1] + I1pop + Rpop + Irpop;
        assert(totalPop>0);
        double intSum = 0;
        double intSum1 = 0;
        for(unsigned int j=0; j<recoveryTimeSteps; j++){
            //linearize to get 0 index
            intSum  += dt * (beta1[j] * I1[0][j] + beta2[j] * Ir[0][j]);
            intSum1 += dt * (gamma1[j] * I1[0][j] + gamma2[j] * Ir[0][j]);
        }
        S[k] = (S[k-1] + dt*delta*totalPop)/(1+dt*(intSum/totalPop+delta));
        double intSum2 =0;
        
        for(unsigned int j=0; j<waningTimeSteps; j++){//columns
            if(j==0) {
                R[1][j] = intSum1;//boundary condition
            } else {
                R[1][j] = (R[0][j]+(dt/dwane)*R[1][j-1])/(1+dt*((1.0/dwane)+(rho[j]/totalPop)*intSum+delta));
            }
            intSum2+= dt*rho[j]*R[1][j];
        }
        double intSum3=0;
        for(unsigned int j=0; j<recoveryTimeSteps; j++){//columns
            if(j==0) {
                //boundary conditions
                I1[1][j] = intSum*S[k]/totalPop;
                Ir[1][j] = (1.0/totalPop)*intSum2*intSum;
            } else {
                I1[1][j] = (I1[0][j]+(dt/dtau)*I1[1][j-1])/(1+dt*((1.0/dtau)+gamma1[j]+delta));
                Ir[1][j] = (Ir[0][j]+(dt/dtau)*Ir[1][j-1])/(1+dt*((1.0/dtau)+gamma2[j]+delta));
            }
            intSum3 += dt*(beta1[j]*I1[0][j]+beta2[j]*Ir[0][j]);
        }
        const double sI = S[k]*intSum3/totalPop;
        symptomaticIncidence[k] = sI;
        asymptomaticIncidence[k] = (1.0/totalPop)*intSum2*intSum3;
        
        
        // convergence bookkeeping
        if (sI >= obs_max.second) obs_max = {k, sI};       // update min and
        if (sI <= obs_min.second) obs_min = {k, sI};       // max obs if necessary
        if (k<eq_interval) {
            // haven't generated enough data yet to know whether we've converged
            obs_deque.push_back(sI);
        } else {
            // have enough data; update min, max as needed and break if max-min < epsilon
            obs_deque.pop_front();
            if (obs_max.first == k - eq_interval) {        // max val is old--need to find a new one
                int offset = eq_interval - 1;
                obs_max = {k - offset, obs_deque.front()}; // forced update
                for (auto e: obs_deque) {                  // scan for a better max
                    if (obs_max.second < e) {
                        obs_max = {k-offset,e};
                    }
                    --offset;
                }
            } else if (obs_min.first == k - eq_interval) { // min val is old--need to find a new one
                int offset = eq_interval - 1;
                obs_min = {k - offset, obs_deque.front()}; // foced update
                for (auto e: obs_deque) {                  // scan for a better min
                    if (obs_min.second > e) {
                        obs_min = {k-offset,e};
                    }
                    --offset;
                }
            }
            if (obs_max.second - obs_min.second < epsilon) {
                // woohoo!
                symptomaticIncidence.resize(k+1);
                asymptomaticIncidence.resize(k+1);
                break;
            } else {
                // aww ...
                obs_deque.push_back(sI);
            }
        }
        // shift just calculated row, making space for next iteration
        for(unsigned int j = 0; j < recoveryTimeSteps; j++){
            I1[0][j] = I1[1][j];
            Ir[0][j] = Ir[1][j];
        }
        for(unsigned int j = 0; j < waningTimeSteps; j++){
            R[0][j]  = R[1][j];
        }
        
    }
    
    /*    cout<<"\nR vec size "<<R.size()<<"\n";
     for (auto v: R) { for (auto e: v) cout << e << " "; }
     cout<<"\nI1 vec size "<<I1.size()<<"\n";
     for (auto v: I1) { for (auto e: v) cout << e << " "; }
     cout<<"\nIr vec size "<<Ir.size()<<"\n";
     for (auto v: Ir) { for (auto e: v) cout << e << " "; }
     cout<<"\nS vec size "<<S.size()<<"\n";
     for (auto e: S) { cout << e << " "; }*/
    //cout<<"\nsymptomatic incidence size "<<symptomaticIncidence.size()<<"\n";
    //for (auto e: symptomaticIncidence) { cout << e << endl; }
    //for (auto e: asymptomaticIncidence) { cout << e << endl; }
    //cout<<endl;
    cerr << "Converged at [" << obs_min.second << ", " << obs_max.second << "] after " << symptomaticIncidence.size() << " observations\n";
    
    return 0;
}
