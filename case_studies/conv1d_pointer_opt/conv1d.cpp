#include <optitrust.h>

#define INLINE_ASM_STUB

void conv1d(int* data, int* kernels, int* out ) {
for (int hi = 0; hi < 1; hi++) {
  for (int tile_j = 0; tile_j < 4; tile_j++) {
    #define output_tile_0 "m7"
    #define output_tile_1 "m6"
    #define output_tile_2 "m5"
    #define output_tile_3 "m4"
    INLINE_ASM_STUB;
    INLINE_ASM_STUB;
    INLINE_ASM_STUB;
    INLINE_ASM_STUB;
    for (int c = 0; c < 4; c++) {
      static int __attribute__((section(".xheep_data_interleaved"))) y[4 * 4];
      for (int j = 0; j < 4; j++) {
          for (int r = 0; r < 4; r++) {
            y[j * 4 + r] = 0;
            if (j + r + 4 * tile_j < 16) {
              y[j * 4 + r] = data[c * 16 + (j + r + 4 * tile_j)];
            }
          }
        }
      #define data_tile "m3"
      INLINE_ASM_STUB;
      #define kernel_tile_0 "m2"
      #define kernel_tile_1 "m1"
      #define kernel_tile_2 "m0"
      INLINE_ASM_STUB;
      INLINE_ASM_STUB;
      INLINE_ASM_STUB;
      INLINE_ASM_STUB;
      #undef kernel_tile_1
      INLINE_ASM_STUB;
      INLINE_ASM_STUB;
      #undef kernel_tile_2
      INLINE_ASM_STUB;
      INLINE_ASM_STUB;
      #undef kernel_tile_0
      #undef data_tile
    }
    INLINE_ASM_STUB;
    #undef output_tile_0
    INLINE_ASM_STUB;
    #undef output_tile_1
    INLINE_ASM_STUB;
    #undef output_tile_2
    INLINE_ASM_STUB;
    #undef output_tile_3
  }
}
}


/* relying on the following instruction..."
rvm_mld(dst,src)
asm volatile("mld.w "{dst_int}", (%1), %0" :: "r"(4*({src}.strides[0])), "r"(&{src_data}));
*/

/* relying on the following instruction..."
rvm_mmasa(md,ms1,ms2)
asm volatile("mmasa.w "{md_int}", "{ms1_int}", "{ms2_int});
*/

/* relying on the following instruction..."
rvm_mst(src,dst)
asm volatile("mst.w "{src_int}", (%1), %0" :: "r"(4*({dst}.strides[0])), "r"(&{dst_data}));
*/

/* relying on the following instruction..."
rvm_mzero(dst)
asm volatile("mzero "{dst_int});
*/
