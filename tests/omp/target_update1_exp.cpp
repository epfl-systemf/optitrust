void init(float*, float*, int);

void init_again(float*, float*, int);

void output(float*, int);

void vec_mult(float* p, float* v1, float* v2, int N) {
  int i;
  init(v1, v2, N);
#pragma omp target data map(to : v1[:N], v2[:N]) map(from : p0 [0:N])
  {
#pragma omp target
#pragma omp parallel for
    for (i = 0; i < N; i++) p[i] = v1[i] * v2[i];
    init_again(v1, v2, N);
#pragma omp target update to(v1[:N], v2[:N])
#pragma omp target
#pragma omp parallel for
    for (i = 0; i < N; i++) p[i] = p[i] + v1[i] * v2[i];
  }
  output(p, N);
}
