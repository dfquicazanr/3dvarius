// Vectores

#include <iostream>
#include <cmath>
using namespace std;
//---------------------- class Vector3D --------------------
struct Vector3D{
  double v[3];
 public:
  void   cargue(double x0, double y0, double z0);
  void   show(void);
  // Funciones de salida de componentes
  double x(void){return v[0];};
  double y(void){return v[1];};
  double z(void){return v[2];};
  //Lectura de Elementos
  double & operator[](int i){return v[i];};

  // Operaciones vectoriales
  Vector3D    operator= (Vector3D v2);
  Vector3D    operator+ (Vector3D v2);
  Vector3D    operator+=(Vector3D v2);
  Vector3D    operator- (Vector3D v2);
  Vector3D    operator-=(Vector3D v2);
  // Producto por escalar
  Vector3D    operator* (double a);
  Vector3D    operator*=(double a);
  friend  Vector3D    operator* (double a,Vector3D v1);	
  // Division por escalar
  Vector3D    operator/ (double a);
  // Producto cruz
  Vector3D    operator^ (Vector3D v2);
  // Producto punto
  double operator* (Vector3D v2);
  // Norma 
  friend  double norma2(Vector3D v1);    
  friend  double norma(Vector3D v1);    
};
// Metodos de la clase Vector3D
void Vector3D::cargue(double x0, double y0, double z0){
  v[0]=x0; v[1]=y0; v[2]=z0;
}
void Vector3D::show(void){
  cout << "(" <<v[0]<< "," <<v[1]<< "," <<v[2]<< ")" << endl;
}
Vector3D Vector3D::operator=(Vector3D v2){
  for(int i=0;i<3;i++)
    v[i] = v2.v[i];
  return *this;
}
Vector3D Vector3D::operator+(Vector3D v2){
  Vector3D total;
  for(int i=0;i<3;i++)
    total.v[i] = v[i] + v2.v[i];
  return total;
}
Vector3D Vector3D::operator+=(Vector3D v2){
  *this = *this + v2;
  return *this;
}
Vector3D Vector3D::operator*(double a){
  Vector3D total;
  for(int i=0;i<3;i++)
    total.v[i] = a*v[i];
  return total;
}
Vector3D Vector3D::operator*=(double a){
  *this = (*this)*a;
  return *this;
}
Vector3D Vector3D::operator/(double a){
  double inver = 1.0/a;
  Vector3D total;
  for(int i=0;i<3;i++)
    total.v[i] = inver*v[i];
  return total;
}
Vector3D Vector3D::operator-(Vector3D v2){
  return *this + v2*(-1); 
}
Vector3D Vector3D::operator-=(Vector3D v2){
  *this = *this - v2;
  return *this;
}
double Vector3D::operator*(Vector3D v2){
  double p=0;
  for(int i=0;i<3;i++)
    p += v[i]*v2.v[i];
  return p;
}
Vector3D Vector3D::operator^(Vector3D v2){
  Vector3D c;
  c.v[0] = v[1]*v2.v[2]-v[2]*v2.v[1];
  c.v[1] = v[2]*v2.v[0]-v[0]*v2.v[2];
  c.v[2] = v[0]*v2.v[1]-v[1]*v2.v[0];
  return c;
}
Vector3D operator*(double a,Vector3D v1){
  Vector3D total;
  total = v1*a;	
  return total;
}
double norma2(Vector3D v1){
  double n=0;
  for(int i=0;i<3;i++)
    n += v1.v[i]*v1.v[i];
  return n;
}
double norma(Vector3D v1){
  return sqrt(norma2(v1));
}
