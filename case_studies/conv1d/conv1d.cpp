#include <optitrust.h>

void rvm_mld(float *src, float *dest, int N1, int N2) {
  for (int t_i = 0; t_i < 4; t_i++) {
    for (int t_j = 0; t_j < 4; t_j++) {
      dest[MINDEX2(4,4,t_i,t_j)] = src[MINDEX2(N1, N2, t_i, t_j)];
    }
  }
}

void rvm_mst(float* src, float* dest, int N1, int N2) {
  for (int t_i = 0; t_i < 4; t_i++) {
    for (int t_j = 0; t_j < 4; t_j++) {
      dest[MINDEX2(N1, N2, t_i, t_j)] = src[MINDEX2(4, 4, t_i, t_j)];
    }
  }
}

void rvm_mzero(float* dest) {
  for (int t_i = 0; t_i < 4; t_i++) {
    for (int t_j = 0; t_j < 4; t_j++) {
      dest[MINDEX2(4, 4, t_i, t_j)] = 0.0f;
    }
  }
}


void rvm_mma(float* C, float* A, float *B) {
  for (int t_i = 0; t_i < 4; t_i++) {
    for (int t_j = 0; t_j < 4; t_j++) {
      for (int t_k = 0; t_k < 4; t_k++) {
        C[MINDEX2(4, 4, t_i, t_j)] += B[MINDEX2(4, 4, t_j, t_k)] * A[MINDEX2(4, 4, t_i, t_k)];
      }
    }
  }
}

#define KW 4
void conv1d(float* O, float* I, float* W, int IW, int IC, int OC) {
  for (int i = 0; i < OC; i++) {
    for (int j = 0; j < IW; j++) {
      float sum = 0.0f;
      for (int k = 0; k < IC; k++) {
        for (int r = 0; r < KW; r++) {
          float y = 0.0;
          if (j + r < IW) {
          y = I[MINDEX2(IC, IW, k, j + r)];
          sum += y * W[MINDEX3(OC, IC, KW, i, k, r)];
          }
        }
      }
      O[MINDEX2(OC, IC, i, j)] = sum;
    }
  }
}


/*
void conv1d(float* O, float* I, float* W, int IW, int IC, int OC, int KW) {
  float x;
  for (int i = 0; i < OC; i++) {
    for (int j = 0; j < IW; j++) {
      for (int k = 0; k < IC; k++) {
        for (int r = 0; r < KW; r++) {
          float y;
        }
      }
    }
  }
}*/

