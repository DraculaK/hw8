#include <iostream>
#include <fstream>
#include <cmath> 

using namespace std;

//----------------------------------------------------------
void sympeuler(double* const p, double* const q, const double& dt, double& H, const double& dim);
//----------------------------------------------------------

int main(){
 
 const int dim=2;
 double q[dim], p[dim], H;
 const double dt = 0.05, pi = 3.141592653589793238, e = 0.6;
 const double t0 = 0, tEnd = 20*pi; 
 const double N = (tEnd-t0)/dt;
 
 //--initial conditions-------------------------------------    
 q[0] = 1-e;
 q[1] = 0;
 p[0] = 0;
 p[1] = sqrt((1+e)/(1-e));
 //---------------------------------------------------------
 double a = sqrt(pow(q[0],2)+pow(q[1],2));
 H = 0.5*(pow(p[0],2)+pow(p[1],2))-(1/a);
 
 ofstream out("data.txt");   
 out << t0 << "\t" << q[0] << "\t" << q[1] << "\t" << p[0] << "\t" << p[1] << "\t" << H << endl;
 
 double t = t0;
 
 for(int i=0; i<N; i++){
  
  sympeuler(p, q, dt, H, dim);
  t += dt;
  out << t << "\t" << q[0] << "\t" << q[1] << "\t" << p[0] << "\t" << p[1] << "\t" << H << endl;
     
 }
 
 out.close();   
 return(0);    
}

//--sympletic Euler method----------------------------------
void sympeuler(double* const p, double* const q, const double& dt, double& H, const double& dim){
 double b = q[0]*q[0]+q[1]*q[1];
 for(int i=0; i<dim; i++){
  p[i] -= dt*q[i]/(pow(b,3.0/2.0));
  q[i] += dt*p[i];
 }
 H = 0.5*(pow(p[0],2)+pow(p[1],2))-(1/sqrt(b));   
}
