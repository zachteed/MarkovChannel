#include "Model.hpp"


namespace Markov {

  int initial_state(Model& m, double vm, double* s)
  {
    double var[] = {1, vm}; int N=m.n_states, P=m.n_prms;
    cblas_dgemv(CblasRowMajor, CblasNoTrans, N, P,
        1.0, m.rs, P, var, 1, 0.0, s, 1);
    cblas_dscal(N, 1.0/cblas_dasum(N, s, 1), s, 1); return 1;
  }

  int transition_matrix(Model& m, double vm, double* Q)
  {

    double var[] = {1, vm};
    int N=m.n_states, E=m.n_edges, P=m.n_prms;

    if(!m.rates) {

      double* B = (double*) malloc(2*E*P*sizeof(double));
      double* rates = (double*) malloc(2*E*P*sizeof(double));
      m.rates = rates; double *x1, *x2, *x3, *x4;

      for (size_t i=0; i<E; i++) {
        int e1=m.edges[i].V1, e2=m.edges[i].V2;
        Math::vdSub(P, &m.rs[e2*P], &m.rs[e1*P], &B[i*P]);
      }
      std::memcpy(&B[P*E], m.rk, P*E*sizeof(double));

      for (size_t i=0; i<E; i++) {
        Math::vdSub(P, &B[P*(i+E)], &B[P*i], &rates[P*2*i]);
        Math::vdAdd(P, &B[P*(i+E)], &B[P*i], &rates[P*(2*i+1)]);
      }
      cblas_dscal(2*E*P, .5, rates, 1); free(B);
    }

    double* expr = (double*) malloc(2*E*sizeof(double));
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 2*E, P, 1.0,
        m.rates, P, var, 1, 0.0, expr, 1);
    Math::vdExp(2*E, expr, expr);

    for(int i = 0; i < 2*E; i++) std::cout<<expr[i]<<std::endl;

    memset(Q, 0, N*N*sizeof(double));
    for (size_t i=0; i<E; i++) {
      int e1 = m.edges[i].V1; int e2 = m.edges[i].V2;
      Q[N*e1+e1] -= expr[2*i]; Q[N*e2+e2] -= expr[2*i+1];
      Q[N*e2+e1] += expr[2*i]; Q[N*e1+e2] += expr[2*i+1];
    }
    free(expr); return 1;
  }
}
