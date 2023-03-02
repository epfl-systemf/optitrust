#include "../../include/optitrust.h"
#include "matmul.h"
#include "omp.h"

void mm(float* C, float* A, float* B, int m, int n, int p) {
  float* Bt = (float*)MALLOC4(n / 32, p / 4, 4, 32, sizeof(float));
#pragma omp parallel for
  for (int bj = 0; bj < n; bj += 32) {
    for (int bk = 0; bk < p; bk += 4) {
      for (int k = 0; k < 4; k++) {
#pragma omp simd
        for (int j = 0; j < 32; j++) {
          Bt[j + bj * p + bk * 32 + k * 32] =
              B[j + bj + n * (k + bk)];
        }
      }
    }
  }
#pragma omp parallel for
  for (int bi = 0; bi < m; bi += 32) {
    for (int bj = 0; bj < n; bj += 32) {
      float* sum = (float*)MALLOC2(32, 32, sizeof(float));
      for (int i = 0; i < 32; i++) {
#pragma omp simd
        for (int j = 0; j < 32; j++) {
          sum[i * 32 + j] = 0.;
        }
      }
      for (int bk = 0; bk < p; bk += 4) {
        for (int i = 0; i < 32; i++) {
#pragma omp simd
          for (int j = 0; j < 32; j++) {
            sum[i * 32 + j] +=
                A[0 + bk + p * (i + bi)] * Bt[j + bj * p + bk * 32 + 0 * 32];
          }
#pragma omp simd
          for (int j = 0; j < 32; j++) {
            sum[i * 32 + j] +=
                A[1 + bk + p * (i + bi)] * Bt[j + bj * p + bk * 32 + 1 * 32];
          }
#pragma omp simd
          for (int j = 0; j < 32; j++) {
            sum[i * 32 + j] +=
                A[2 + bk + p * (i + bi)] * Bt[j + bj * p + bk * 32 + 2 * 32];
          }
#pragma omp simd
          for (int j = 0; j < 32; j++) {
            sum[i * 32 + j] +=
                A[3 + bk + p * (i + bi)] * Bt[j + bj * p + bk * 32 + 3 * 32];
          }
        }
      }
      for (int i = 0; i < 32; i++) {
#pragma omp simd
        for (int j = 0; j < 32; j++) {
          C[j + bj + n * (i + bi)] = sum[i * 32 + j];
        }
      }
      free(sum);
    }
  }
  free(Bt);
}
