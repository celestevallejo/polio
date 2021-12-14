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
#include <cassert>

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
    Polio_boxcar_model(double bd, vector<double> infect, vector<double> recover, vector<double> wane, int totalBins): DiffEq_Sim(totalBins), mu(bd), b(infect), g(recover), r(wane){}
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
    int counter = 0;
    void derivative(double const x[], double dxdt[]) {
        //key:
        //S -> x[0]     R_1' -> x[7]
        //I_1 -> x[1]   R_2'-> x[8]
        //I_2 -> x[2]   R_3'-> x[9]
        //I_3 -> x[3]   R_1'' -> x[10]
        //R_1 -> x[4]   R_2''-> x[11]
        //R_2-> x[5]   R_3''-> x[12]
        //R_3-> x[6]
        double totalPop = 0;
        double infectPop = 0;
        for(unsigned int i = 0; i < nbins; i++){
            assert(x[i] >=0);
            totalPop += x[i];
            if(i >0 and i < 4){
                infectPop += x[i];
            }
        }
        dxdt[0] = mu*totalPop - (b[0]*infectPop/totalPop + mu)*x[0];
        dxdt[1] = (b[0]*x[0] + b[3]*x[7] + b[6]*x[10])*infectPop/totalPop - (g[0] + mu)*x[1];
        dxdt[2] = (b[1]*x[4] + b[4]*x[8] + b[7]*x[11])*infectPop/totalPop - (g[1] + mu)*x[2];
        dxdt[3] = (b[2]*x[5] + b[5]*x[9] + b[8]*x[12])*infectPop/totalPop - (g[2] + mu)*x[3];
        dxdt[4] = g[0]*x[1] - (b[1]*infectPop/totalPop + r[0] + mu)*x[4];
        dxdt[5] = g[1]*x[2] - (b[2]*infectPop/totalPop + r[1] + mu)*x[5];
        dxdt[6] = g[2]*x[3] - (r[2] + mu)*x[6];
        dxdt[7] = r[0]*x[4] - (b[3]*infectPop/totalPop + r[3] + mu)*x[7];
        dxdt[8] = r[1]*x[5] - (b[4]*infectPop/totalPop + r[4] + mu)*x[8];
        dxdt[9] = r[2]*x[6] - (b[5]*infectPop/totalPop + r[5] + mu)*x[9];
        dxdt[10] = r[3]*x[7] - (b[6]*infectPop/totalPop + mu)*x[10];
        dxdt[11] = r[4]*x[8] - (b[7]*infectPop/totalPop + mu)*x[11];
        dxdt[12] = r[5]*x[9] - (b[8]*infectPop/totalPop + mu)*x[12];
        
    }
    void printX(){
        for(unsigned int i=0; i < nbins; i++){
            cout << x[i] << " ";
        }
        cout << endl;
        delete [] x;
        x=NULL;
    }
    
};


#endif /* Polio_boxcar_model_h */
