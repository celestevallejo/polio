#include <cmath>
#include <iostream>
#include <string>
#include <fstream>

#include <gsl/gsl_integration.h>

// compiled on roc using:
// g++ -O2 -std=c++11 -Wall -Wextra -Wno-deprecated-declarations --pedantic -I/home/tjhladish/work/AbcSmc/gsl_local/include/ integration_test2.cpp -o integ_test -lm -L/home/tjhladish/work/AbcSmc/gsl_local/lib/ -lgsl -lgslcblas -lpthread -ldl

using namespace std;

double integ_pi1(const double tau_max);
double integ_pi2(const double tau_max);
double integ_pi3(const double tau_max);
double integ_R0(const double tau_max);
double integ_gamma1pi1(const double tau_max);
double integ_gamma2pi3(const double tau_max);
double R0(double tau, void *params);
double pi1(double tau, void *params);
double pi2(double tau, void *params);
double pi3(double tau,void *params);

struct f_params {
    double y0=1;
    double b0=1;
    double mu=.021;
    double mu0=.02;
    double K=0.0001;
    double delta=0.02;
    double contactRate1 =197.0/365;
    double contactRate2 = 150.0/365;
    double recover1=10.0/365;
    double recover2=15.0/365;
    double waneRate=0.02/365;
    double y1=100;
    double r=1.7;
    double nu=15;
    double t1=(1/mu)*log(y1/y0);
    double c = b0*(mu-mu0)/(y0*(exp((mu-mu0)*t1)-1));
    double popSize = 10;
    double T = t1;//infectivity time
    double S = 200; //susceptiblility time

};

double beta1(double tau, void *params){
    f_params &p= *reinterpret_cast<f_params *>(params);
    return p.contactRate1*(p.b0*exp(p.mu0*tau)-(p.c*p.y0*(exp(p.mu*tau)-exp(p.mu0*tau))/(p.mu-p.mu0)))/((p.b0*exp(p.mu0*tau)-(p.c*p.y0*(exp(p.mu*tau)-exp(p.mu0*tau))/(p.mu-p.mu0)))+p.K);
}
double gamma1(double tau, void *params){
    f_params &p= *reinterpret_cast<f_params *>(params);
    return p.recover1*p.y0*exp(p.mu*tau)/(p.b0*exp(p.mu0*tau)-(p.c*p.y0*(exp(p.mu*tau)-exp(p.mu0*tau)/(p.mu-p.mu0)))+p.K);
}
double gamma2(double tau, void *params){
    f_params &p= *reinterpret_cast<f_params *>(params);
    return p.recover2*p.y0*exp(p.mu*tau)/(p.b0*exp(p.mu0*tau)-(p.c*p.y0*(exp(p.mu*tau)-exp(p.mu0*tau)/(p.mu-p.mu0)))+p.K);
}
double gamma1pi1(double tau,void *params){
    f_params &p= *reinterpret_cast<f_params *>(params);
    return gamma1(tau,&p)*pi1(tau,&p);
}
double gamma2pi3(double tau,void *params){
    f_params &p= *reinterpret_cast<f_params *>(params);
    return gamma2(tau,&p)*pi3(tau,&p);
}
double pi1(double tau, void *params) {
    f_params &p= *reinterpret_cast<f_params *>(params);
    //solved the integral below symbollically using mathematica
    return exp(-(p.K*p.recover1*tau+tau*p.delta+(p.recover1*(tau*p.mu0-log(-p.c*exp(tau*p.mu)*p.y0+exp(tau*p.mu0)*(p.c*p.y0+p.b0*(p.mu-p.mu0)))+log(p.b0*(p.mu-p.mu0))))/p.c));
}
double eq_I1(double tau, void *params){
    f_params &p = *reinterpret_cast<f_params *>(params);
    return (1.0/10.0)*pi1(tau,&p);
    //return p.mu*p.popSize*(1.0-1.0/integ_R0(p.T))*pi1(tau,&p);
}
double pi2(double tau, void *params){
    f_params &p= *reinterpret_cast<f_params *>(params);
    //approximate the function using a Taylor series and then integrate that
    //integrated taylor series is summed below
    //double nTerm=0;
    //double taylorSum=0;
    //double tol = 1e-2;
    //int n =0;//series counter

    /*int factorial;
    for(int n=0;n<5;n++){
        if(n==0){
            factorial=1;
        }
        else{
            factorial*=n;
        }

        taylorSum += pow(p.waneRate*p.mu,n)*pow(p.y1,n*p.r-1)*pow(p.nu,n)*pow(tau,n+1)/factorial;//this assumes concentration of antibodies is never zero
    }

    return exp(-taylorSum*(R0(tau,&p))+p.delta);*/
    
    //below is the integration of the correct taylor series polynomial for the waning function with constant in the denominator (constant prevents the function from being undefined)...since I am only taking the taylor poly up to the second derivative, the best the accuracy can be is O(step size)^2 (big O notation). I'm breaking it up into terms of the eqn so its easier to read
    
    double firstTerm = p.mu*p.waneRate*(integ_R0(p.T)-1)*tau/(p.y1+p.K);
    double secondTerm = p.mu*p.waneRate*(integ_R0(p.T)-1)*pow(p.y1,p.r)*p.nu*pow(tau,2)/(2*pow(p.y1+p.K,2));
    double thirdTerm = p.mu*p.waneRate*(integ_R0(p.T)-1)*pow(p.y1,2*p.r-1)*pow(p.nu,2)*(2*p.y1-(p.r*(p.y1+p.K)))*pow(tau,3)/(3*2*pow(p.y1+p.K,3));
    //fourthTerm not derived using taylor series since we only need to integrate the constant mu
    double fourthTerm = p.mu*tau;

    return exp(-(firstTerm+secondTerm+thirdTerm+fourthTerm));
    
}
double eq_R(double tau,void *params){
    f_params &p = *reinterpret_cast<f_params *>(params);
    return (1.0/10.0)*pi2(tau,&p);
    //return (p.mu*p.popSize*(1.0-1.0/integ_R0(p.T))*integ_gamma1pi1(p.T)+p.popSize*integ_gamma2pi3(p.T)*((1.0-1.0/integ_R0(p.T))*(1.0-p.mu*(integ_pi1(p.T)+integ_gamma1pi1(p.T)*integ_pi2(p.S))))/integ_gamma2pi3(p.T)*integ_pi2(p.S)+integ_pi3(p.T))*pi2(tau,&p);
}
double pi3(double tau, void *params){
    f_params &p= *reinterpret_cast<f_params *>(params);
    //solved the integral below symbollically using mathematica
    return exp(-(p.K*p.recover2*tau+tau*p.delta+(p.recover2*(tau*p.mu0-log(-p.c*exp(tau*p.mu)*p.y0+exp(tau*p.mu0)*(p.c*p.y0+p.b0*(p.mu-p.mu0)))+log(p.b0*(p.mu-p.mu0))))/p.c));
}
double eq_Ir(double tau, void *params){
    f_params &p= *reinterpret_cast<f_params *>(params);
    return (1.0/10.0)*pi3(tau,&p);
    //return p.popSize*(1.0-1.0/integ_R0(p.T))*(1.0-p.mu*(integ_pi1(p.T)+integ_gamma1pi1(p.T)*integ_pi2(p.S)))/(integ_gamma2pi3(p.T)*integ_pi2(p.S)+integ_pi3(p.T))*pi3(tau,&p);
}

double R0(double tau, void *params){
    f_params &p= *reinterpret_cast<f_params *>(params);
    return beta1(tau,&p)*pi1(tau,&p);//need to integrate this term to get R0
}



double integ_eqI1(const double tau) {
    f_params params;


    gsl_function F;
    F.function = &eq_I1;
    F.params = reinterpret_cast<void *>(&params);

    double result, error;
    size_t neval;

    const double xlow=0;
    const double xhigh=tau;
    const double epsabs=1e-4;
    const double epsrel=1e-4;

    gsl_integration_qng (&F,
            xlow,
            xhigh,
            epsabs,
            epsrel,
            &result,
            &error,
            &neval);
    return result;
}


double integ_eqR(const double tau) {
    f_params params;


    gsl_function F;
    F.function = &eq_R;
    F.params = reinterpret_cast<void *>(&params);

    double result, error;
    size_t neval;

    const double xlow=0;
    const double xhigh=tau;
    const double epsabs=1e-4;
    const double epsrel=1e-4;

    gsl_integration_qng (&F,
            xlow,
            xhigh,
            epsabs,
            epsrel,
            &result,
            &error,
            &neval);
    return result;
}


double integ_eqIr(const double tau) {
    f_params params;

    gsl_function F;
    F.function = &eq_Ir;
    F.params = reinterpret_cast<void *>(&params);

    double result, error;
    size_t neval;

    const double xlow=0;
    const double xhigh=tau;
    const double epsabs=1e-4;
    const double epsrel=0;

    gsl_integration_qng (&F,
            xlow,
            xhigh,
            epsabs,
            epsrel,
            &result,
            &error,
            &neval);
    return result;
}


double integ_gamma1pi1(const double tau) {
    f_params params;

    gsl_function F;
    F.function = &gamma1pi1;
    F.params = reinterpret_cast<void *>(&params);

    double result, error;
    size_t neval;

    const double xlow=0;
    const double xhigh=tau;
    const double epsabs=1e-4;
    const double epsrel=0;

    gsl_integration_qng (&F,
            xlow,
            xhigh,
            epsabs,
            epsrel,
            &result,
            &error,
            &neval);
    return result;
}


double integ_gamma2pi3(const double tau) {
    f_params params;


    gsl_function F;
    F.function = &gamma2pi3;
    F.params = reinterpret_cast<void *>(&params);

    double result, error;
    size_t neval;

    const double xlow=0;
    const double xhigh=tau;
    const double epsabs=1e-4;
    const double epsrel=0;

    gsl_integration_qng (&F,
            xlow,
            xhigh,
            epsabs,
            epsrel,
            &result,
            &error,
            &neval);
    return result;
}


double integ_R0(const double tau) {
    f_params params;

    gsl_function F;
    F.function = &R0;
    F.params = reinterpret_cast<void *>(&params);

    double result, error;
    size_t neval;

    const double xlow=0;
    const double xhigh=tau;
    const double epsabs=1e-1;
    const double epsrel=1e-1;

    gsl_integration_qng (&F,
            xlow,
            xhigh,
            epsabs,
            epsrel,
            &result,
            &error,
            &neval);
    return result;
}


double integ_pi1(const double tau) {
    f_params params;

    gsl_function F;
    F.function = &pi1;
    F.params = reinterpret_cast<void *>(&params);

    double result, error;
    size_t neval;

    const double xlow=0;
    const double xhigh=tau;
    const double epsabs=1e-4;
    const double epsrel=0;

    gsl_integration_qng (&F,
            xlow,
            xhigh,
            epsabs,
            epsrel,
            &result,
            &error,
            &neval);
    return result;
}


double integ_pi2(const double tau) {
    f_params params;

    gsl_function F;
    F.function = &pi2;
    F.params = reinterpret_cast<void *>(&params);

    double result, error;
    size_t neval;

    const double xlow=0;
    const double xhigh=tau;
    const double epsabs=1e-1;
    const double epsrel=1e-1;

    gsl_integration_qng (&F,
            xlow,
            xhigh,
            epsabs,
            epsrel,
            &result,
            &error,
            &neval);
    return result;
}


double integ_pi3(const double tau) {
    f_params params;

    gsl_function F;
    F.function = &pi3;
    F.params = reinterpret_cast<void *>(&params);

    double result, error;
    size_t neval;

    const double xlow=0;
    const double xhigh=tau;
    const double epsabs=1e-1;
    const double epsrel=1e-1;

    gsl_integration_qng (&F,
            xlow,
            xhigh,
            epsabs,
            epsrel,
            &result,
            &error,
            &neval);
    return result;
}


int main(void) {
    ofstream myfile;
    myfile.open("/Users/Celeste/Desktop/polio_pde_test_R0=10.csv",ios_base::app);
    const double step = 1.0;
    f_params params;
    /*for (double tau_max = 0; tau_max < params.T; tau_max += step){
        cout<<"T "<<params.T<<"\n";
        //cout<<"integ_pi1 "<<integ_pi1(params.T)<<"\n";

    }*/
    double i1 = integ_eqI1(params.T);
    double r=integ_eqR(params.T);
    double ir=integ_eqIr(params.T);
    double r0 = integ_R0(params.T);
    //double b1 = beta1(tau_max,&params);
    //myfile<<b1<<"\n";
    //myfile<<"R0"<<" , "<<"I1"<<" , "<<"R"<<" , "<<"Ir "<<" , "<<"mu"<<" , "<<"mu0"<<" , "<<"b0"<<" , "<<"y0"<<" , "<<"c"<<" , "<<"t1"<<" , "<<"y1"<<" , "<<"r"<<" , "<<"nu"<<" , "<<"contact1"<<" , "<<"recover1"<<" , "<<"popsize"<<"\n";
    //myfile<<"\n";
    //myfile<<r0<<" , "<<i1<<" , "<<r<<" , "<<ir<<", "<<params.mu<<" , "<<params.mu0<<" , "<<params.b0<<" , "<<params.y0<<" , "<<params.c<<" , "<<params.t1<<" , "<<params.y1<<" , "<<params.r<<" , "<<params.nu<<" , "<<params.contactRate1<<" , "<<params.recover1<<" , "<<params.popSize<<"\n";
    cout<<"R0 "<<r0<<"\n";
    cout<<"I1 "<<i1<<"\n";
    cout<<"I1(0) "<<params.mu*params.popSize*(1.0-1.0/integ_R0(params.T))<<"\n";
    cout<<"R "<<r<<"\n";
    cout<<"R(0) "<<(params.mu*params.popSize*(1.0-1.0/integ_R0(params.T))*integ_gamma1pi1(params.T)+params.popSize*integ_gamma2pi3(params.T)*((1.0-1.0/integ_R0(params.T))*(1.0-params.mu*(integ_pi1(params.T)+integ_gamma1pi1(params.T)*integ_pi2(params.S))))/integ_gamma2pi3(params.T)*integ_pi2(params.S)+integ_pi3(params.T))<<"\n";
    cout<<"Ir "<<ir<<"\n";
    cout<<"Ir(0) "<<params.popSize*(1.0-1.0/integ_R0(params.T))*(1.0-params.mu*(integ_pi1(params.T)+integ_gamma1pi1(params.T)*integ_pi2(params.S)))/(integ_gamma2pi3(params.T)*integ_pi2(params.S)+integ_pi3(params.T))<<"\n";
    cout<<"S "<<params.popSize/r0<<"\n";
    cout<<(1/r0)+(i1/params.popSize)+(r/params.popSize)+(ir/params.popSize)<<"\n";
    //myfile<<"R0 = "<<integ_R0(28)<<"\n";
    myfile.close();

}
