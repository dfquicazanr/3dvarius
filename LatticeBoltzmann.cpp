#include <iostream>

using namespace std;

double wD1Q3[3]       = {4.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0};
int    vD1Q3[1][3]    = {{0, 1, -1}};
double wD2Q5[5]       = {1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0};
int    vD2Q5[2][5]    = {{0, 0, -1, 1, 0},
                         {0, -1, 0, 0, 1}};
double wD2Q9[9]       = {4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
int    vD2Q9[2][9]    = {{0, 0, -1, 1, 0, -1, 1, -1, 1},
                         {0, -1, 0, 0, 1, -1, -1, 1, 1}};
double wD3Q7[7]       = {1.0 / 4.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0};
int    vD3Q7[3][7]    = {{0, 0, 0, -1, 1, 0, 0},
                         {0, 0, -1, 0, 0, 1, 0},
                         {0, -1, 0, 0, 0, 0, 1}};
double wD3Q19[19]     = {1.0 / 3.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
int    vD3Q19[3][19]  = {{0, 0, 0, -1, 1, 0, 0, 0, -1, 0, 1, -1, 0, 0, 1, -1, 1, -1, 1},
                         {0, 0, -1, 0, 0, 1, 0, -1, 0, 1, 0, 0, -1, 1, 0, -1, -1, 1, 1},
                         {0, -1, 0, 0, 0, 0, 1, -1, -1, -1, -1, 1, 1, 1, 1, 0, 0, 0, 0}};

void printArray(const double *array, int size) {
  for (int i = 0; i < size; ++i) {
    cout << array[i] << endl;
  }
}

int[] flattenMatrix() {

}

void printMatrix(const double *array, int size1, int size2) {
  for (int i = 0; i < size1; ++i) {
    for (int j = 0; j < size2; ++j) {
      cout << array[i][i];
    }
    cout << endl;
  }
}

class LatticeBoltzmann {
private:
  double *w;
  int *v;
public:
  LatticeBoltzmann(int dim, int q);
};

LatticeBoltzmann::LatticeBoltzmann(int dim, int q) {
  if (dim == 1 && q == 3) {
    w = wD1Q3;
    v = vD1Q3;
  } else
  if (dim == 2 && q == 5) {
    w = wD2Q5;
    v = vD2Q5;
  } else
  if (dim == 2 && q == 9) {
    w = wD2Q9;
    v = vD2Q9;
  } else
  if (dim == 3 && q == 7) {
    w = wD3Q7;
    v = vD3Q7;
  } else
  if (dim == 4 && q == 19) {
    w = wD3Q19;
    v = vD3Q19;
  }
  printArray(w, q);
  printMatrix(v, dim, q);
}

int main() {
  LatticeBoltzmann * waves = new LatticeBoltzmann(2, 5);
}
//
//
//class LatticeBoltzmann {
//private:
//  double *w;
//  int **v;                          // v[0][i]=v_ix, v[1][i]=v_iy
//  double ***f, ***fNew;   // f[ix][iy][i] => f_i in cell (ix, iy)
//public:
//  LatticeBoltzmann(int dim, int q);
//  double rho(int ix, int iy, bool useNew);
//  Vector2D j(int ix, int iy, bool useNew);
//  double fEquilibrium(double rhoValue, Vector2D jValue, int i);
//  void init(double rhoValue, Vector2D jValue);
//  void collide(void);
//  void advection(void);
//  void imposeField(int t);
//  void print(const char * fileName);
//};
//
//LatticeBoltzmann::LatticeBoltzmann(int dim, int q) {
//  if (dim == 1 && q == 3) {
//    w = wD1Q3;
//
//  } else
//  if (dim == 2 && q == 5) {
//
//  } else
//  if (dim == 2 && q == 9) {
//
//  } else
//  if (dim == 3 && q == 7) {
//
//  } else
//  if (dim == 4 && q == 19) {
//
//  }
//}