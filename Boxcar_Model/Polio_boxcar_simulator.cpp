//
//  Polio_boxcar_simulator.cpp
//  
//
//  Created by Celeste on 12/10/17.
//
//

#include <stdio.h>
#include "adjustable_boxcar_model.h"
//#include "Polio_boxcar_model_extended_waning.h" //uncomment for explicit ode model
#include "DIFFEQ_SIM.h"
#include <math.h>
#include "Polio_boxcar_parameters.h"
#include <vector>

void usage(){
    cout<<"\n\t Usage: Input model initialization <S> <I1> <I2> <I3> <sum of all R compartments>\n";
    exit(-1);
}

int main(int argc, char** argv){
    if(argc!= 6) usage();
    
    vector<double> initialValues(argc-1);
    for(int i = 1; i<argc; i++){
        initialValues[i-1] = atof(argv[i]);
    }
    
    //Number of compartments set in parameter file
    
    //set initial antibody level
    setIntAnt();
    
    //set peak antibody level
    setPeakAnt();
    
    //set infectious period
    setInfectPeriod();
    
    //vectors for parameters
    vector<double> beta(sbins + rbins-1);//1 recovered compartment has complete immunity
    vector<double> gamma(ibins);
    vector<double> fastWane(ibins);
    vector<double> slowWane(ibins);
    vector<double> recovery(rbins);
    
    for(unsigned int i = 0; i < beta.size(); i++){
        //beta[i] = numContacts*probInfection(pow(10,i));
        //beta[i] = numContacts*((i+1)/2.0);
        beta[i]=10;//toy number for simulation validation
    }
    for(unsigned int i = 0; i < gamma.size(); i++){
        gamma[i] = 1.0/getInfectPeriod(i);
        //const double baseline_recovery_rate = 365.0/28;
        //gamma[i] = baseline_recovery_rate*(1+i); // gamma increases linearly with sequential infections (i)
        //gamma[i] = 1;//toy number for simulation validation
    }
    for(unsigned int i = 0; i < fastWane.size(); i++){
        fastWane[i] = 0.01;//toy number for simulation validation
        //fastWane[i] = fastWaningRate(i);
    }
    for(unsigned int i = 0; i < slowWane.size(); i++){
        slowWane[i] = 0.01;//toy number for simulation validation
        //slowWane[i] = slowWaningRate(i);
    }
    for(unsigned int i =0; i < recovery.size(); i++){
        recovery[i] = 0.01;//toy number for simulation validation
    }
    
    
    
    //Constuct model and parameters
    adjustable_boxcar_model model(birthRate,beta,gamma,recovery,fastWane,slowWane,nbins,ibins,rbins,subrbins);
    //Polio_boxcar_model model(birthRate,beta,gamma,recovery,nbins); //uncomment for explicit ode model
    
    
    //Initialize model using user inputs
    model.initialize(initialValues);
    //model.initialize(initialValues,rbins);//uncomment for explicit ode model
    
    //Run model
    model.run_simulation();
    
    //Output simulation results
    model.printX();
    
    
    return 0;
}
