typedef unsigned int uint;

typedef const double cdouble;

typedef struct { uint x; uint y; } vect;

typedef vect vect; // typedef to be removed

typedef int[2][2] mat2x2;

typedef int*** mat3d; // typedef to be removed

int main() {
   uint x;
   uint v[3];
   cdouble y1, t2;
   vect v;
   mat2x2 m;
   mat3d M;
   mat3d* T;
}

