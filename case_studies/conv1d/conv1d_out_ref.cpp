#include <optitrust.h>

void rvm_mld(float* src, float* dest, int N1, int N2) {
  for (int t_i = 0; t_i < 4; t_i++) {
    for (int t_j = 0; t_j < 4; t_j++) {
      dest[MINDEX2(4, 4, t_i, t_j)] = src[MINDEX2(N1, N2, t_i, t_j)];
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
      dest[MINDEX2(4, 4, t_i, t_j)] = 0.f;
    }
  }
}

void rvm_mma(float* C, float* A, float* B) {
  for (int t_i = 0; t_i < 4; t_i++) {
    for (int t_j = 0; t_j < 4; t_j++) {
      for (int t_k = 0; t_k < 4; t_k++) {
        C[MINDEX2(4, 4, t_i, t_j)] +=
            B[MINDEX2(4, 4, t_j, t_k)] * A[MINDEX2(4, 4, t_i, t_k)];
      }
    }
  }
}

void conv1d(float* O, float* I, float* W, int IW, int IC, int OC) {
  for (int tile_i_hi = 0; tile_i_hi < exact_div(exact_div(OC, 4), 4);
       tile_i_hi++) {
    int data_base = 0;
    for (int tile_j = 0; tile_j < exact_div(IW, 4); tile_j++) {
      {
        float* const sum0 = (float*)MALLOC2(4, 4, sizeof(float));
        float* const sum1 = (float*)MALLOC2(4, 4, sizeof(float));
        float* const sum2 = (float*)MALLOC2(4, 4, sizeof(float));
        float* const sum3 = (float*)MALLOC2(4, 4, sizeof(float));
        rvm_mzero(&sum0[0]);
        rvm_mzero(&sum1[0]);
        rvm_mzero(&sum2[0]);
        rvm_mzero(&sum3[0]);
        int data_row = 0;
        for (int k = 0; k < IC; k++) {
          {
            float* const y = (float*)MALLOC2(4, 4, sizeof(float));
            for (int j = 0; j < 4; j++) {
              int ofs = data_base + j + 0;
              int drow_ofs = data_row + ofs;
              for (int r = 0; r < 4; r++) {
                {
                  {
                    y[0 + j * 4 + r] = 0.f;
                    ;
                    ;
                    if (ofs < IW) {
                      y[0 + j * 4 + r] = I[drow_ofs];
                    }
                  }
                  ofs = ofs + 1;
                }
                drow_ofs = drow_ofs + 1;
              }
            }
            float* const data_tile = (float*)MALLOC2(4, 4, sizeof(float));
            rvm_mld(y, data_tile, 4, 4);
            float* const kernel_tile = (float*)MALLOC2(4, 4, sizeof(float));
            rvm_mld(&W[4 * k + 4 * tile_i_hi * 4 * 4 * IC], kernel_tile, OC,
                    4 * IC);
            rvm_mma(&sum0[0], kernel_tile, data_tile);
            MFREE2(4, 4, kernel_tile);
            float* const kernel_tile1 = (float*)MALLOC2(4, 4, sizeof(float));
            rvm_mld(&W[4 * k + (4 * tile_i_hi + 1) * 4 * 4 * IC], kernel_tile1,
                    OC, 4 * IC);
            rvm_mma(&sum1[0], kernel_tile1, data_tile);
            MFREE2(4, 4, kernel_tile1);
            float* const kernel_tile2 = (float*)MALLOC2(4, 4, sizeof(float));
            rvm_mld(&W[4 * k + (4 * tile_i_hi + 2) * 4 * 4 * IC], kernel_tile2,
                    OC, 4 * IC);
            rvm_mma(&sum2[0], kernel_tile2, data_tile);
            MFREE2(4, 4, kernel_tile2);
            float* const kernel_tile3 = (float*)MALLOC2(4, 4, sizeof(float));
            rvm_mld(&W[4 * k + (4 * tile_i_hi + 3) * 4 * 4 * IC], kernel_tile3,
                    OC, 4 * IC);
            rvm_mma(&sum3[0], kernel_tile3, data_tile);
            MFREE2(4, 4, kernel_tile3);
            MFREE2(4, 4, data_tile);
            MFREE2(4, 4, y);
          }
          data_row = data_row + IW;
        }
        rvm_mst(&sum0[0], &O[4 * tile_j + 4 * tile_i_hi * 4 * IC], OC, IC);
        rvm_mst(&sum1[0], &O[4 * tile_j + (4 * tile_i_hi + 1) * 4 * IC], OC,
                IC);
        rvm_mst(&sum2[0], &O[4 * tile_j + (4 * tile_i_hi + 2) * 4 * IC], OC,
                IC);
        rvm_mst(&sum3[0], &O[4 * tile_j + (4 * tile_i_hi + 3) * 4 * IC], OC,
                IC);
        MFREE2(4, 4, sum0);
        MFREE2(4, 4, sum1);
        MFREE2(4, 4, sum2);
        MFREE2(4, 4, sum3);
      }
      data_base = data_base + 4;
    }
  }
}
