//
//  SIR_H.h
//  
//
//  Created by Celeste on 10/1/17.
//
//

#ifndef Polio_boxcar_model_h
#define Polio_boxcar_model_h

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include "DIFFEQ_SIM.h"
#include <math.h>
#include <vector>
#include <array>

using namespace std;

class Polio_boxcar_model : public DiffEq_Sim {
    
private:
    //birth/death rate
    const double mu;
    //effective contact rates
    vector<double> b;
    //recovery rates
    vector<double> g;
    //waning rates
    vector<double> r;
    
    
public:
    Polio_boxcar_model(): mu(0.0){
        nbins = 0;
        b = {0};
        g = {0};
        r = {0};
    }
    Polio_boxcar_model(double bd, vector<double> infect, vector<double> recover, vector<double> wane, int totalBins): mu(bd){
        nbins = totalBins;
        b.resize(infect.size());
        g.resize(recover.size());
        r.resize(wane.size());
        for(unsigned int i = 0; i < infect.size(); i++){
            b[i] = infect[i];
        }
        for(unsigned int i = 0; i < recover.size(); i++){
            g[i] = recover[i];
        }
        for(unsigned int i = 0; i < wane.size(); i++){
            r[i] = wane[i];
        }
        
    }
    ~Polio_boxcar_model() {};
    
    void initialize(vector<double>initialValues, int rbins) {
        x = new double[nbins];
        for(unsigned int i = 0; i < nbins; i++){
            if(i >= (initialValues.size()-1)){
                x[i] = initialValues.back()/(double)rbins;//last element in initial values vector is sum of all R compartments
            }
            else{
                x[i] = initialValues[i];
            }
        }
    }
    
    void derivative(double const x[], double dxdt[]) {
        double totalPop = 0;
        double infectPop = 0;
        for(unsigned int i = 0; i < nbins; i++){
            totalPop += x[i];
            if(i >0 and i < 5){
                infectPop += x[i];
            }
        }
        dxdt[0] = mu*totalPop - (b[0]*infectPop/totalPop + mu)*x[0];
        dxdt[1] = b[0]*x[0]*infectPop/totalPop - (g[0] + mu)*x[1];
        dxdt[2] = b[1]*x[5]*infectPop/totalPop - (g[1] + mu)*x[2];
        dxdt[3] = b[2]*x[6]*infectPop/totalPop - (g[2] + mu)*x[3];
        dxdt[4] = b[3]*x[7]*infectPop/totalPop - (g[3] + mu)*x[4];
        dxdt[5] = g[0]*x[1] - (b[1]*infectPop/totalPop + mu)*x[5] + r[0]*x[6];
        dxdt[6] = g[1]*x[2] - (b[2]*infectPop/totalPop + r[0] + mu)*x[6] + r[1]*x[7];
        dxdt[7] = g[2]*x[3] - (b[3]*infectPop/totalPop + r[1] + mu)*x[7] + r[2]*x[8];
        dxdt[8] = g[3]*x[4] - (r[2] + mu)*x[8];
        
    }
    void printX(){for(unsigned int i=0; i < nbins; i++) { cout << x[i] << " ";} cout << endl; }
    
};


#endif /* Polio_boxcar_model_h */
