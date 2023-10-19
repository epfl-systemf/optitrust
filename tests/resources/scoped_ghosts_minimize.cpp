#include "../../include/optitrust.h"

int main() {
  __pure();

  int x = 0;

  float* const A = (float* const) MALLOC1(10, sizeof(float));
  __GHOST_BEGIN(focusA, matrix1_ro_focus, "A, 0");
  x += 1;
  x += A[MINDEX1(10, 0)];
  x += 1;
  __GHOST_END(focusA);
  __GHOST_BEGIN(focusA2, matrix1_ro_focus, "A, 0");
  x += 1;
  __GHOST_END(focusA2);
  MFREE1(10, A);

  float* const B = (float* const) MALLOC2(8, 6, sizeof(float));
  // FIXME: __GHOST_BEGIN(focusBi, group_focus, "i := 1, bound_check_start := checked, bound_check_stop := checked, bound_check_step := checked");
  x += 2;
  __GHOST_BEGIN(focusBi, group_focus, "i := 1, bound_check_start := checked, bound_check_stop := checked, bound_check_step := checked");
  __GHOST_BEGIN(focusBj, group_focus, "i := 2, bound_check_start := checked, bound_check_stop := checked, bound_check_step := checked");
  x += 1;
  x += B[MINDEX2(8, 6, 1, 2)];
  x += 1;
  __GHOST_END(focusBj);
  __GHOST_END(focusBi);
  x += 2;
  // FIXME: __GHOST_END(focusBi);
  /* FIXME: RO joins are resource usage
  __GHOST_BEGIN(focusBi, group_ro_focus, "i := 1, bound_check_start := checked, bound_check_stop := checked, bound_check_step := checked");
  x += 2;
  __GHOST_BEGIN(focusBj, group_ro_focus, "i := 2, bound_check_start := checked, bound_check_stop := checked, bound_check_step := checked");
  x += 1;
  x += B[MINDEX2(8, 6, 1, 2)];
  x += 1;
  __GHOST_END(focusBj);
  x += 2;
  __GHOST_END(focusBi);
  */
  MFREE2(8, 6, B);

}
