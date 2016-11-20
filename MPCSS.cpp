#include<iostream>
#include<cmath>
using namespace std;
#include<stdio.h>

int times = 0;
double m1 = 0.5;
double l1 = 0.15;
int NumInput = 1;
int NumState = 2;
double derutaT = 0.001;
int step =8;
double p = 1000;
double q = 10;
double r = 0.08;
double s = 0.0001;
double xref = 0.0;
double yref = 0.15;
//double ita = 10;
//int opt = 1000;
int period = 150;
  
int main(){
  
  FILE *MPCdata1;
  if((MPCdata1 = fopen("MPCSSdata1.txt", "w+" ))==NULL){
    cout<<"打开文件失败!";  
  }
  FILE *MPCdata2;
  if((MPCdata2 = fopen("MPCSSdata2.txt", "w+" ))==NULL){
    cout<<"打开文件失败!";  
  }
  FILE *MPCdata3;
  if((MPCdata3 = fopen("MPCSSdata3.txt", "w+" ))==NULL){
    cout<<"打开文件失败!";  
  }

  double thetad = 30;
  double omegal = 0;
  double thetal = 1.0*thetad*M_PI/180;
  double u1;
  double U[NumInput][step];
  double HU[NumInput][step];
  for(int i = 0; i< NumInput; i++){
    for(int j = 0; j<step; j++){
      U[i][j] = 0;
      HU[i][j] = 0;
    }
  }
  double fx = l1*cos(thetal);
  double fy = l1*sin(thetal);
  double X[NumState][step+1];
  double Ramuda[NumState][step+1];
 
  for(int i = 0; i< NumState; i++){
    for(int j = 0; j<step+1; j++){
      X[i][j] = 0;
      Ramuda[i][j] = 0;
    }
  }

  double huaix[2];
  double ramuda1;
  double ramuda2;
  double normHU = 1;
  //MPC parameter
  double MPCX[NumState][period+1];
  double MPCHX[NumState][period+1];
  double MPCU[NumInput][period+1]; 
  for(int i = 0; i< NumState; i++){
    for(int j = 0; j<period+1; j++){
      MPCX[i][j] = 0;
      MPCHX[i][j] = 0;
    }
  }
  for(int i = 0; i< NumInput; i++){
    for(int j = 0; j<period+1; j++){
      MPCU[i][j] = 0;
    }
  }
  MPCX[0][0] = thetal;
  MPCX[1][0] = omegal;
  MPCHX[0][0] = fx;
  MPCHX[1][0] = fy;

  //Get Optimal U
  for(int k = 0; k<period; k++){

    X[0][0] = MPCX[0][k];
    X[1][0] = MPCX[1][k];
    thetal = X[0][0];
    omegal = X[1][0];
    for(int j = 0; j<step; j++){
       u1 = U[0][j];
       double f[2];
       f[0] = omegal;
       f[1] = u1/(m1*l1*l1);
       X[0][j+1] = X[0][j]+derutaT*f[0];
       X[1][j+1] = X[1][j]+derutaT*f[1];
       thetal = X[0][j+1];
       omegal = X[1][j+1];
    }
    huaix[0] = -2*(l1*cos(thetal)-xref)*p*l1*sin(thetal)+2*(l1*sin(thetal)-yref)*p*l1*cos(thetal);
    huaix[1] = 2*omegal*s;
    Ramuda[0][step] = huaix[0];
    Ramuda[1][step] = huaix[1];
      
    ramuda1 =  Ramuda[0][step];
    ramuda2 =  Ramuda[1][step];
      
    for(int j = 0; j<step; j++){
	//u1 = U[0][step-j-1];//should care number 1
	thetal = X[0][step-j-1];
	omegal = X[1][step-j-1];
	double Hx[2];
	Hx[0] = -2*(l1*cos(thetal)-xref)*q*l1*sin(thetal)+2*(l1*sin(thetal)-yref)*q*l1*cos(thetal);
	Hx[1] = ramuda1;
	Ramuda[0][step-j-1] = Ramuda[0][step-j]+derutaT*Hx[0];//1212121212
	Ramuda[1][step-j-1] = Ramuda[1][step-j]+derutaT*Hx[1];
	ramuda1 = Ramuda[0][step-j-1];
	ramuda2 = Ramuda[1][step-j-1];
    }

    //SS method
    double SA[8][8];
    double Sb[8];
    for(int i=0; i<step; i++){
      for(int j=0; j<step; j++){
	if(i == j){
	   SA[i][j] = 2.0*r;
	}
	else{
	  SA[i][j] = 0;
	}
      }
    }
    for(int i=0; i<step; i++){
      Sb[i] = 1;
    }
    for(int j = 0; j<step; j++){
       	
       ramuda1 = Ramuda[0][j+1];
       ramuda2 = Ramuda[1][j+1];
       Sb[j] = -Sb[j]*ramuda2/(m1*l1*l1);
    }
    double P_tem[10][10];
    double P[10];
    for(int i=0; i<step; i++){
      for(int j=0; j<step; j++){
	if(i == j){
	  P_tem[i][j] = 1.0/(2.0*r);
	}
	else{
	  P_tem[i][j] = 0;
	}
      }
    }
    for(int i=0; i<step; i++){
      P[i] = P_tem[i][i]*Sb[i];
      U[0][i] = P[i];
    }
    ///////////////////
    u1 = U[0][0];
    omegal = MPCX[1][k];
    double f[2];
    f[0] = omegal;
    f[1] = u1/(m1*l1*l1);
    MPCX[0][k+1] = MPCX[0][k]+derutaT*f[0];
    MPCX[1][k+1] = MPCX[1][k]+derutaT*f[1];
    MPCU[0][k] = u1;
    
    fprintf(MPCdata3, "%f\n", MPCU[0][k]);
  }

  for(int i = 0; i<period+1; i++){
    thetal = MPCX[0][i];
    MPCHX[0][i] = l1*cos(thetal);
    MPCHX[1][i] = l1*sin(thetal);
    fprintf(MPCdata1, "%f\n", MPCHX[0][i]);
    fprintf(MPCdata2, "%f\n", MPCHX[1][i]);
  }
  double NORM[1][period+1];
  for(int i = 0; i<period+1; i++){
    NORM[0][i] = 0;
  }
  
  for(int i = 0; i<period+1; i++){
    thetal = MPCX[0][i];
    NORM[0][i] = sqrt((l1*cos(thetal)-xref)*(l1*cos(thetal)-xref)+(l1*sin(thetal)-yref)*(l1*sin(thetal)-yref));
    //fprintf(MPCdata1, "%f\n", NORM[1][i]);
  }
  cout<<"Finial!"<<times<<endl;
  return 1;
}
