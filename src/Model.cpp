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

    this->rs = (double*) malloc(P*(G.N+G.E)*sizeof(double));
    this->rk = this->rs + P*G.N;
    init_params(G.N+G.E, this->rs);

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
  }


  Model* neighbor(Model& m, int n=1) {

    Model* n = new Model();
    n->G = m.G;

    std::vector<int> node_idx;
    std::vector<int> edge_idx;

    for (int i=0; i<G.N; i++)
      node_idx.push_back(i);

    for (int i=0; i<G.E; i++) {
      edge_idx.push_back(i);
    }

    double rng[4];
    rng=uniform(4, rng);

    if (prms.has_add_edge() && rng[0] < prms.mutation().add_edge()) {
      Graph::add_edge(n->G, false);

      while(edge_idx.size() < n->G.E)
        edge_idx.push_back(-1);
    }

    if (prms.has_add_node() && rng[1] < prms.mutation().add_node()) {
      Graph::add_node(n->G, false);

      while(edge_idx.size() < n->G.E)
        edge_idx.push_back(-1);

      while(node_idx.size() < n->G.N)
        node_idx.push_back(-1);
    }

    if (prms.has_rm_edge() && rng[2] < prms.mutation().rm_edge()) {
      Graph::rm_edge(n->G, eidx, true);

      while(edge_idx.size() < n->G.E)
        edge_idx.push_back(-1);

      while(node_idx.size() < n->G.N)
        node_idx.push_back(-1);
    }

    if (prms.has_rm_node() && rng[3] < prms.mutation().rm_node()) {
      Graph::rm_node(n->G, nidx, eidx, true);

      while(edge_idx.size() < n->G.E)
        edge_idx.push_back(-1);

      while(node_idx.size() < n->G.N)
        node_idx.push_back(-1);
    }

    int N = n->G.N, E = n->G.E;
    int P = prms.n_prms()

    n->rs = (double*) malloc(P*(N+E)*sizeof(double));
    n->rk = n->rs + P*N;
    init_params(N+E, n->rs);

    int* nidx = &node_idx[0];
    int* eidx = &edge_idx[0];

    for (int i=0; i<N; i++) {
      if (nidx[i] > 0)
        memcpy(n->rs[P*i], m.rs[P*nidx[i]], P*sizeof(double));
    }
    for (int i=0; i<E; i++) {
      if (eidx[i] > 0) {
        memcpy(n->rk[P*i], m.rk[P*eidx[i]], P*sizeof(double));
      }
    }

    int* mask = (int*) malloc(N+E*sizeof(int));
    double* mult = (double*) malloc(N+E*sizeof(double));
    double* r = (double*) malloc(N+E*sizeof(double));

    Math::rng_int(N+E, mask, 0, 2)
    for(int i=0; i<N+E; i++) *(mult++) = *(mask++);

    double sig = prms.mutation().sig();
    Math::rng_gaussian(N+E, r, 0, sig);

    vdMult(N+E, r, mult, r);
    vdAdd(N+E, r, n->rs, n->rs);

    free(mask);
    free(mult);
    free(r);

    return n;
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
