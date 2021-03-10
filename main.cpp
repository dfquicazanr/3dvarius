#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

const int Lx = 128;
const int Ly = 128;

const int Q = 5;
const double W0 = 1.0 / 3.0;

const double c = 0.5;
const double c2 = c*c;
const double fEquilibriumAuxConst = 1.0 - 3.0 * c2 * (1.0 - W0);

const double tau = 0.5;
const double uTau = 1.0 / tau;
const double uMUTau = 1.0 - uTau;

// Help structs
struct Vector2D {
  double x;
  double y;
};

// Lattice Boltzmann Class
class LatticeBoltzmann {
private:
  double w[Q];
  int v[2][Q];                            // v[0][i]=v_ix, v[1][i]=v_iy
  double f[Lx][Ly][Q], fNew[Lx][Ly][Q];   // f[ix][iy][i] => f_i in cell (ix, iy)
public:
  LatticeBoltzmann(void);
  double rho(int ix, int iy, bool useNew);
  Vector2D j(int ix, int iy, bool useNew);
  double fEquilibrium(double rhoValue, Vector2D jValue, int i);
  void init(double rhoValue, Vector2D jValue);
  void collide(void);
  void advection(void);
  void imposeField(int t);
  void print(const char * fileName);
};
// Lattice Boltzmann Class method definitions
LatticeBoltzmann::LatticeBoltzmann(void) {
  // Load weights
  w[0] = W0;
  w[1] = w[2] = w[3] = w[4] = (1.0 - W0) / 4.0;
  // Load speed vectors
  // v_0          v_1             v_2             v_3             v_4
  v[0][0] = 0;    v[0][1] = 1;    v[0][2] = 0;    v[0][3] = -1;   v[0][4] = 0;    //v_x
  v[1][0] = 0;    v[1][1] = 0;    v[1][2] = 1;    v[1][3] = 0;    v[1][4] = -1;   //v_y
}

double LatticeBoltzmann::rho(int ix, int iy, bool useNew = false) {
  double sum = 0;
  for (int i = 0; i < Q; i++) {
    sum += useNew ? fNew[ix][iy][i] : f[ix][iy][i];
  }
  return sum;
}

//Vector2D LatticeBoltzmann::j(int ix, int iy, bool useNew = false) {
//  Vector2D sum;
//  sum.x = 0.0;
//  sum.y = 0.0;
//  for (int i = 0; i < Q; i++) {
//    sum.x += useNew ? v[0][i] * fNew[ix][iy][i] : v[0][i] * f[ix][iy][i];
//    sum.y += useNew ? v[1][i] * fNew[ix][iy][i] : v[1][i] * f[ix][iy][i];
//  }
//  return sum;
//}

Vector2D LatticeBoltzmann::j(int ix, int iy, bool useNew = false) {
  Vector2D sum;
  sum.x = 0.0;
  sum.y = 0.0;
  for (int i = 0; i < Q; i++) {
    if (useNew) {
      sum.x += v[0][i] * fNew[ix][iy][i];
      sum.y += v[1][i] * fNew[ix][iy][i];
    } else {
      sum.x += v[0][i] * f[ix][iy][i];
      sum.y += v[1][i] * f[ix][iy][i];
    }
  }
  return sum;
}


double LatticeBoltzmann::fEquilibrium(double rhoValue, Vector2D jValue, int i) {
  return i == 0 ?
         rhoValue * fEquilibriumAuxConst :
         3.0 * w[i] * (c2 * rhoValue + v[0][i] * jValue.x + v[1][i] * jValue.y);
}

void LatticeBoltzmann::init(double rhoValue, Vector2D jValue) {
  for (int ix = 0; ix < Lx; ix++) {
    for (int iy = 0; iy < Ly; iy++) {
      for (int i = 0; i < Q; i++) {
        f[ix][iy][i] = fEquilibrium(rhoValue, jValue, i);
      }
    }
  }
}

void LatticeBoltzmann::collide(void) {
  double rhoValue;
  Vector2D jValue;

  for (int ix = 0; ix < Lx; ix++) {
    for (int iy = 0; iy < Ly; iy++) {
      rhoValue = rho(ix, iy);
      jValue = j(ix, iy);
      for (int i = 0; i < Q; i++) {
        fNew[ix][iy][i] = uMUTau * f[ix][iy][i] + uTau * fEquilibrium(rhoValue, jValue, i);
      }
    }
  }
}

void LatticeBoltzmann::advection(void ) {
  for (int ix = 0; ix < Lx; ix++) {
    for (int iy = 0; iy < Ly; iy++) {
      for (int i = 0; i < Q; i++) {
        int nextIx = (ix + v[0][i] + Lx) % Lx;
        int nextIy = (iy + v[1][i] + Ly) % Ly;
        f[nextIx][nextIy][i] = fNew[ix][iy][i];
      }
    }
  }
}

void LatticeBoltzmann::imposeField(int t) {
  int ix = Lx / 2, iy = Ly / 2;
  double waveLength = 10.0;
  double omega = 2.0 * M_PI * c / waveLength;
  double rho0 = 10.0 * sin(omega * t);
  Vector2D j0 = j(ix,iy);
  for(int i = 0; i < Q; i++) {
    fNew[ix][iy][i] = fEquilibrium(rho0, j0, i);
  }
}

void LatticeBoltzmann::print(const char * fileName) {
  ofstream outFile(fileName);
  double rho0;
  for(int ix = 0; ix < Lx; ix++) {
    for(int iy = 0; iy < Ly; iy++) {
      rho0 = rho(ix, iy, true);
      outFile << ix << " " << iy << " " << rho0 << endl;
    }
    outFile << endl;
  }
  outFile.close();
}

int main(void) {
  LatticeBoltzmann waves;
  double rho0 = 0;
  Vector2D j0;
  j0.x = 0.0;
  j0.y = 0.0;
  int tMax = 100;

  waves.init(rho0, j0);
  for (int t = 0; t < tMax; ++t) {
    waves.collide();
    waves.imposeField(t);
    waves.advection();
  }
  waves.print("waves_2021-1.dat");
  return 0;

}