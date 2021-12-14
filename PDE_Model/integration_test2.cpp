#include <cmath>
#include <iostream>

#include <gsl/gsl_integration.h>

// compiled on roc using:
// g++ -O2 -std=c++11 -Wall -Wextra -Wno-deprecated-declarations --pedantic -I/home/tjhladish/work/AbcSmc/gsl_local/include/ integration_test2.cpp -o integ_test -lm -L/home/tjhladish/work/AbcSmc/gsl_local/lib/ -lgsl -lgslcblas -lpthread -ldl

using namespace std;

struct f_params {
    double y0;
    double b0;
    double mu;
    double mu0;
    double c;
    double K;
    double delta;
    double r;
    double nu;
    double t1;
    double contact1;
    double contact2;
    double recover1;
    double recover2;
    double waneRate;
    double y1;
};

double f(double tau, string str, void *params) {
    f_params &p= *reinterpret_cast<f_params *>(params);
    double boostPat = p.b0*exp(p.mu0*tau) - (p.c*p.y0*(exp(p.mu*tau)-exp(p.mu0*tau))/(p.mu-p.mu0));
    double boostAnt = p.y0*exp(p.mu*tau);
    double waneAnt = p.y1*pow((1+(p.r-1)*pow(p.y1,(p.r-1)*p.nu*(tau-p.t1)),(-1/(p.r-1))));
    
    double beta1 = p.contact1*boostPat/(boostPat*p.K);
    double beta2 = p.contact2*boostPat/(boostPat*p.K);
    double gamma1 = p.recover1*boostAnt/(boostPat+p.K);
    double gamma2 = p.recover2*boostAnt/(boostPat+p.K);
    double wane = p.waneRate/(waneAnt+p.K);
                              
    if(str == "R0"){
        return beta1*exp(-integ(tau,"pi1"));
                                  
    }
    else if(str=="pi1"){
        return gamma1 + p.delta;
                                  
    }
                              
    //return p.alpha*p.y0*exp(p.mu*tau)/(p.b0*exp(p.mu0*tau)-(p.c*p.y0*(exp(p.mu*tau)-exp(p.mu0*tau)/(p.mu-p.mu0)))+p.K);
                
}

void integ(const double tau,const string str) {
    f_params params;
    params.y0    = 1.0;
    params.b0    = 1.0;
    params.mu    = 0.318;
    params.mu0   = 0.16;
    params.c     = 0.01;
    params.K     = 0.0001;
    params.delta = 0.02;
    params.recover1 = 13.0/365;
    params.contact1=135.0/365;

    gsl_function F;
    F.function = &f;
    F.params = reinterpret_cast<void *>(&params);

    double result, error;
    size_t neval;

    const double xlow=0;
    const double xhigh=tau;
    const double epsabs=1e-4;
    const double epsrel=1e-4;

    int code=gsl_integration_qng (&F, 
            xlow,
            xhigh,
            epsabs, 
            epsrel, 
            &result,
            &error, 
            &neval);
    if(code) {
        cerr<< "There was a problem with integration: code " << code << endl;
    } else {
        cout << "tau, Result " << tau << ", " << exp(-result) << " +/- " << exp(-error) << ", " << error << " from " << neval << " evaluations" << endl;
    }
}

int main(void) {
    const double tau_max = 28;
    const double step = 1.0;
    for (double tau = 0; tau < tau_max; tau += step) integ(tau,"R0");
}
