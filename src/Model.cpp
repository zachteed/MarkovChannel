#include "Model.hpp"


namespace Model {

  void init_params(int N, float* r)
  {
    int P = prms.n_prms();
    Math::rng_gaussian(N*P, r);

    for (int i=0; i<P; i++) {
      double mu=prms.mu()[i], sig=prms.std()[i];
      for (int j=0; j<N; j++)
        rs[j*P+i] = sig*rs[j*P+i] + mu;
    }

  }

  Model::Model(int N, double p)
  {
    random_connected_graph(N, this->G, p);
    int P = prms.n_prms();

    this->rs = (double*) malloc(P*G.N*sizeof(double));
    this->rk = (double*) malloc(P*G.E*sizeof(double));

    init_params(G.N, this->rs);
    init_params(G.E, this->rk);

    this->C = (double*) calloc(G.N, sizeof(double));
    this->F = (double*) calloc(G.N, sizeof(double));
    C[0] = 1.0; F[0] = 1.0; r_vec = NULL;
  }

  Model::~Model()
  {
    if (r_vec)
      free(r_vec)

    free(C);
    free(F);
    free(rs);
    free(rk);
  }

  int add_edge(Model& m)
  {
    int P = prms.n_prms();
    if (prms.has_add_edge() && Math.rng_uniform()<prms.add_edge()) {

      int e0 = m.E; Graph::add_edge(m.graph, false);
      double* rk_tmp = (double*) malloc(P*G.N*sizeof(double));

      memcpy(rk_tmp, m.rk, e0*sizeof(double));
      m.rk = rk_tmp;

      if (m.E > e0) {
        init_params(m.E-e0, rk_tmp + e0*P);
      }
    }
  }

  int add_node(Model& m)
  {
    int P = prms.n_prms();
    if (prms.has_add_node() && Math.rng_uniform()<prms.add_node()) {

      int n0 = m.N; Graph::add_node(m.graph);
      double* rs_tmp = (double*) malloc(P*G.N*sizeof(double));

      memcpy(rk_tmp, m.rk, e0*sizeof(double));
      m.rk = rk_tmp;

      if (m.E > e0) {
        init_params(m.E-e0, rk_tmp + e0*P);
      }
    }
  }

  int rm_edge(Model& m) {}

  int rm_node(Model& m) {}

  int neighbor(Model& m, vector<Model>& neighbors, int n=1) {

    neighbors = vector<Model>()
    for (int i = 0)

    n.graph = m.graph;
    int* idx;

    if (prms.has_add_edge() && Math.rng_uniform() < prms.add_edge) {
      Graph::add_edge(n.graph, false);
    }

    if (prms.has_add_node() && Math.rng_uniform() < prms.add_node) {
      Graph::add_node(n.graph, false);
    }

    if (prms.has_rm_edge() && Math.rng_uniform() < prms.rm_edge) {
      Graph::rm_edge(n.graph, idx, true);
    }

    if (prms.has_rm_node() && Math.rng_uniform() < prms.rm_node) {
      Graph::rm_node(n.graph, idx, true);
    }
  }
}

  void initial_state(Model& m, double vm, double* s)
  {
    double var[] = {1, vm/100.0};
    int N=m.G.N, P=prms.n_prms();

    cblas_dgemv(CblasRowMajor, CblasNoTrans, N, P,
        1.0, m.rs, P, var, 1, 0.0, s, 1);
    cblas_dscal(N, 1.0/cblas_dasum(N, s, 1), s, 1
  }


  double* rate_vector(Model& m)
  {
    if ( !m.r_vec ) {
      int N=m.G.N, E=m.G.E, P=prms.n_prms();
      double *b_mat = (double*) malloc(2*E*P*sizeof(double));
      m.r_vec = (double*) malloc(2*E*P*sizeof(double));

      for ( int i=0; i<E; i++ ) {
        int e1=m.edges[i].V1, e2=m.edges[i].V2;
        vdSub(P, &m.rs[e2*P], &m.rs[e1*P], &b_mat[i*P]);
      }
      memcpy(&b_mat[P*E], m.rk, P*E*sizeof(double));

      for ( int i=0; i<E; i++ ) {
        vdSub(P, &B[P*(i+E)], &B[P*i], &r_vev[P*2*i]);
        vdAdd(P, &B[P*(i+E)], &B[P*i], &r_vec[P*(2*i+1)]);
      }
      cblas_dscal(2*E*P, .5, rates, 1); free(b_mat);
    }

    return m.r_vec;
  }


  double* transition_matrix(Model& m, double vm, double* Q)
  {

    double var[] = {1, vm/100.0};
    int N=m.G.N, E=m.G.E, P=m.n_prms();

    double *r_vec = rate_vector(m);
    double *e_vec = (double*) malloc(2*E*sizeof(double));

    cblas_dgemv(CblasRowMajor, CblasNoTrans, 2*E, P, 1.0,
        m.rates, P, var, 1, 0.0, expr, 1);

    vdExp(2*E, expr, expr);
    memset(Q, 0, N*N*sizeof(double));

    for (size_t i=0; i<E; i++) {
      int e1 = m.edges[i].V1; int e2 = m.edges[i].V2;
      Q[N*e1+e1] -= expr[2*i]; Q[N*e2+e2] -= expr[2*i+1];
      Q[N*e2+e1] += expr[2*i]; Q[N*e1+e2] += expr[2*i+1];
    }
    free(expr); return Q;
  }
