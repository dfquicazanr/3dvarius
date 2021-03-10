#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=128;
const int Ly=128;

const int Q=5;
const double W0=1.0/3;

const double C=0.5; // C < sqrt(1.0/(3(1-W0)))
const double C2=C*C;
const double AUX0=1-3*C2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1.0-Utau;
//--------------------- Clase LatticeBoltzmann ------------

//---------------- class LatticeBoltzmann ----------------
class LatticeBoltzmann{
private:
  double w[Q];
  int v[2][Q]; // v[0][i]=V_ix  , v[1][i]=V_iy
  double f[Lx][Ly][Q],fnew[Lx][Ly][Q]; // f[ix][iy][i] =f_i en la celda ix,iy
public:
  LatticeBoltzmann(void);
  double rho(int ix,int iy,bool UseNew);
  double Jx(int ix,int iy,bool UseNew);
  double Jy(int ix,int iy,bool UseNew);
  double feq(double rho0,double Jx0,double Jy0,int i);
  void Inicie(double rho0,double Jx0,double Jy0);
  void Colisione(void);
  void ImponerCampos(int t);
  void Adveccione(void);
  void Imprimase(const char * NombreArchivo);
};
//--------- Implementación de las funciones de la clase

LatticeBoltzmann::LatticeBoltzmann(void){//OJO: AQUI ESTA EL ERROR !!
  //cargar los pesos
  w[0]=W0; w[1]=w[2]=w[3]=w[4]=(1-w[0])/4;
  //cargar los vectores velocidad
  v[0][0]=0;
  v[1][0]=0;

  v[0][1]=1;  v[0][2]=0;  v[0][3]=-1; v[0][4]=0;
  v[1][1]=0;  v[1][2]=1;  v[1][3]=0;  v[1][4]=-1;
}

double LatticeBoltzmann::rho(int ix,int iy,bool UseNew){
  double suma; int i;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]; else suma+=f[ix][iy][i];
  return suma;
}
double LatticeBoltzmann::Jx(int ix,int iy,bool UseNew){//OJO: AQUI HAY UN ERROR !
  double suma; int i;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=v[0][i]*fnew[ix][iy][i]; else suma+=v[0][i]*f[ix][iy][i];
  return suma;
}
double LatticeBoltzmann::Jy(int ix,int iy,bool UseNew){//OJO: AQUI HAY UN ERROR !
  double suma; int i;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=v[1][i]*fnew[ix][iy][i]; else suma+=v[1][i]*f[ix][iy][i];
  return suma;
}
double LatticeBoltzmann::feq(double rho0,double Jx0,double Jy0,int i){
  double resultado;
  if(i==0)
    resultado=rho0*AUX0;
  else
    resultado=3*w[i]*(C2*rho0+v[0][i]*Jx0+v[1][i]*Jy0);
  return resultado;
}
void LatticeBoltzmann::Inicie(double rho0,double Jx0,double Jy0){
  for(int ix=0;ix<Lx;ix++) //para cada celda
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++) //en cada dirección
        f[ix][iy][i]=feq(rho0,Jx0,Jy0,i); //cargue com funcion el valor de equilibrio
}
void LatticeBoltzmann::Colisione(void){
  int ix,iy,i; double rho0,Jx0,Jy0;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++){
      //calculo los campos
      rho0=rho(ix,iy,false); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
      for(i=0;i<Q;i++) //en cada dirección
        fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(rho0,Jx0,Jy0,i);
    }
}
void LatticeBoltzmann::ImponerCampos(int t){
  int ix,iy,i; double rho0,Jx0,Jy0;
  double lambda=10,omega,A=10;
  omega=2*M_PI*C/lambda;
  rho0=A*sin(omega*t);  Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
  ix=Lx/2; iy=Ly/2;
  for(i=0;i<Q;i++)
    fnew[ix][iy][i]=feq(rho0,Jx0,Jy0,i);
}
void LatticeBoltzmann::Adveccione(void){
  for(int ix=0;ix<Lx;ix++) //para cada celda
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++) //en cada dirección
        f[(ix+v[0][i]+Lx)%Lx][(iy+v[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
}
void LatticeBoltzmann::Imprimase(const char * NombreArchivo){
  ofstream MiArchivo(NombreArchivo); double rho0;
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,true);
      MiArchivo<<ix<<" "<<iy<<" "<<rho0<<endl;
    }
    MiArchivo<<endl;
  }
  MiArchivo.close();
}

//---------------- Funciones Globales ----------------

//---------------- Programa Principal -----------
int main(void){
  LatticeBoltzmann Ondas;
  double rho0=0, Jx0=0, Jy0=0;
  int t,tmax=100;

  Ondas.Inicie(rho0,Jx0,Jy0);
  for(t=0;t<tmax;t++){
    Ondas.Colisione();
    Ondas.ImponerCampos(t);
    Ondas.Adveccione();
  }
  Ondas.Imprimase("Ondas2021-I.dat");

  return 0;
}