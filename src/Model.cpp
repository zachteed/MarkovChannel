#include "../include/Model.hpp"


namespace Model {

  int init_state(Model& m, double vm, double* s)
  {
    double[] var = {1, vm}; int N=m.n_states, P=m.prms;
    Math::dgemv(CblasNoTrans, N, P, 1.0, m.rs, 0.0, var, s);
    Math::exp(N, s, s); double sum = Math::a_sum(N, s);
    Math::dscale(N, 1.0/sum, s, s); return 1;
  }

  int transition_matrix(Model& m, double vm, double* Q)
  {

    double[] var = {1, vm};
    int N=m.n_states, E=m.n_edges, P=m.prms;

    if(!m.transition_matrix) {

      double* B = (double*) malloc(N*P*std::sizeof(double));
      double* rates = (double*) malloc(P*N*std::sizeof(double));

      for (size_t i=0; i<E; i++) {
        int e1=m.edges[i].V1, e2=m.edges[i].V2;
        Math::sub(P, e1, e2, B[i*P]);
      }
      std::memcpy(B[P*E], m.rk, P*E);

      for (size_t i=0; i<E; i++) {
        Math::sub(P, B[i*P], B[P*(E+i)], rates[i*P]);
        Math::add(P, B[i*P], B[P*(E+i)], rates[i*P]);
      }
      Math::dscale(2*E*P, .5, rates, rates);

      double* t_matrix = (double *) calloc(P*N*N, sizeof(double));
      m.transition_matrix = t_matrix; double *x1,*x2;

      for (size_t i=0; i<E; i++) {
        int e1 = m.edges[i].V1, e2 = m.edges[i].V2;
        x1=t_matrix[P*(N*e1+e1)]; x2=t_matrix[P*(N*e2+e1)]
        Math::sub(P, x1, R[e1], x1);
        Math::add(P, x2, R[e1], x2);
      }
      std::free(rates);std::free(B);
    }

    Math::dgemv(CblasNoTrans, N, P, 1.0, m.rs, 0.0, var, Q)
    std::free(rates); std::free(B); return 1;
  }
}
