#include "Model.hpp"



namespace Model {

  int Model::Model::count = 0;
  MarkovChannel::ModelParameter prms;

  Model::Model(int N, double p)
  {
    random_connected_graph(N, this->G, p);
    int P = prms.n_prms();

    this->rs = (double*) malloc(P*(G.N+G.E)*sizeof(double));
    this->rk = this->rs + P*G.N;
    std::cout << prms.mu() << '\t' << prms.std() << std::endl;
    Math::rng_gaussian((G.N+G.E)*P, rs, prms.mu(), prms.std());

    this->C = (double*) calloc(G.N, sizeof(double));
    this->F = (double*) calloc(G.N, sizeof(double));
    C[0] = 1.0; F[0] = 1.0; r_vec = NULL;

    id = Model::count++;
  }

  Model::Model()
  {
    r_vec = NULL; rs=NULL; C = NULL; F = NULL;
    id = Model::count++;
  }

  Model::~Model()
  {
    if (r_vec)
      free(r_vec);

    if (C)
      free(C);

    if (F)
      free(F);

    if (rs)
      free(rs);

  }


  Model* neighbor(Model& m, int num) {

    Model* n = new Model();
    n->G = m.G;

    std::vector<int> node_idx;
    std::vector<int> edge_idx;
    int *nidx, *eidx;

    for (int i=0; i<n->G.N; i++)
      node_idx.push_back(i);

    for (int i=0; i<n->G.E; i++) {
      edge_idx.push_back(i);
    }

    double rng[4];
    Math::rng_uniform(4, rng);

    const MarkovChannel::MutationParameter& mut = prms.mutation();

    if (mut.has_add_edge() && rng[0] < mut.add_edge()) {
      Graph::add_edge(n->G, false);

      while(edge_idx.size() < n->G.E)
        edge_idx.push_back(-1);
    }

    if (mut.has_add_node() && rng[1] < mut.add_node()) {
      Graph::add_node(n->G, false);

      while(edge_idx.size() < n->G.E)
        edge_idx.push_back(-1);

      while(node_idx.size() < n->G.N)
        node_idx.push_back(-1);
    }

    if (mut.has_rm_edge() && rng[2] < mut.rm_edge()) {

      nidx = &node_idx[0];
      eidx = &edge_idx[0];
      Graph::rm_edge(n->G, eidx, true);

      while(edge_idx.size() < n->G.E)
        edge_idx.push_back(-1);

      while(node_idx.size() < n->G.N)
        node_idx.push_back(-1);
    }

    if (mut.has_rm_node() && rng[3] < mut.rm_node()) {

      nidx = &node_idx[0];
      eidx = &edge_idx[0];
      Graph::rm_node(n->G, nidx, eidx, true);

      while(edge_idx.size() < n->G.E)
        edge_idx.push_back(-1);

      while(node_idx.size() < n->G.N)
        node_idx.push_back(-1);
    }

    int N = n->G.N, E = n->G.E;
    int P = prms.n_prms();

    n->rs = (double*) malloc(P*(N+E)*sizeof(double));
    n->rk = n->rs + P*N;

    Math::rng_gaussian((N+E)*P, n->rs, prms.mu(), prms.std());

    nidx = &node_idx[0];
    eidx = &edge_idx[0];

    for (int i=0; i<N; i++) {
      if (nidx[i] > 0)
        memcpy(&n->rs[P*i], &m.rs[P*nidx[i]], P*sizeof(double));
    }
    for (int i=0; i<E; i++) {
      if (eidx[i] > 0) {
        memcpy(&n->rk[P*i], &m.rk[P*eidx[i]], P*sizeof(double));
      }
    }

    int* mask = (int*) malloc(N+E*sizeof(int));
    double* mult = (double*) malloc(N+E*sizeof(double));
    double* r = (double*) malloc(N+E*sizeof(double));
    for(int i=0; i<N+E; i++) *(mult++) = *(mask++);

    Math::rng_int(N+E, mask, 0, 2);
    Math::rng_gaussian(N+E, r, 0, mut.std());

    vdMul(N+E, r, mult, r);
    vdAdd(N+E, r, n->rs, n->rs);

    free(mask); free(mult); free(r);
    return n;
  }

  double* initial_state(Model& m, double vm, double* s)
  {
    double var[] = {1, vm/100.0};
    int N=m.G.N, P=prms.n_prms();

    cblas_dgemv(CblasRowMajor, CblasNoTrans, N, P,
        1.0, m.rs, P, var, 1, 0.0, s, 1);
    cblas_dscal(N, 1.0/cblas_dasum(N, s, 1), s, 1);
    return s;
  }


  double* rate_vector(Model& m)
  {
    if ( !m.r_vec ) {
      int N=m.G.N, E=m.G.E, P=prms.n_prms();
      double *b_mat = (double*) malloc(2*E*P*sizeof(double));
      m.r_vec = (double*) malloc(2*E*P*sizeof(double));

      for ( int i=0; i<E; i++ ) {
        int e1=m.G.edges[i].V1, e2=m.G.edges[i].V2;
        vdSub(P, &m.rs[e2*P], &m.rs[e1*P], &b_mat[i*P]);
      }
      memcpy(&b_mat[P*E], m.rk, P*E*sizeof(double));

      for ( int i=0; i<E; i++ ) {
        vdSub(P, &b_mat[P*(i+E)], &b_mat[P*i], &m.r_vec[P*2*i]);
        vdAdd(P, &b_mat[P*(i+E)], &b_mat[P*i], &m.r_vec[P*(2*i+1)]);
      }
      cblas_dscal(2*E*P, .5, m.r_vec, 1); free(b_mat);
    }

    return m.r_vec;
  }


  double* transition_matrix(Model& m, double vm, double* Q)
  {

    double var[] = {1, vm/100.0};
    int N=m.G.N, E=m.G.E, P=prms.n_prms();

    double *r_vec = rate_vector(m);
    double *e_vec = (double*) malloc(2*E*sizeof(double));

    cblas_dgemv(CblasRowMajor, CblasNoTrans, 2*E, P, 1.0,
        r_vec, P, var, 1, 0.0, e_vec, 1);

    vdExp(2*E, e_vec, e_vec);
    memset(Q, 0, N*N*sizeof(double));

    for (size_t i=0; i<E; i++) {
      int e1 = m.G.edges[i].V1; int e2 = m.G.edges[i].V2;
      Q[N*e1+e1] -= e_vec[2*i]; Q[N*e2+e2] -= e_vec[2*i+1];
      Q[N*e2+e1] += e_vec[2*i]; Q[N*e1+e2] += e_vec[2*i+1];
    }
    free(e_vec); return Q;
  }


  std::ostream& operator<<(std::ostream& os, const Model& m)
  {
    int N=m.G.N, E=m.G.E, P=prms.n_prms(), n;
    os << "Model Id:\t" << m.id << "\n" << m.G;

    char buffer[12];

    os << "\n~RS~" << std::endl;
    for (int i=0; i<N; i++) {
      for (int j=0; j<P; j++) {
        n = sprintf(buffer, "%8.4f\t", m.rs[i*P+j]);
        os << std::string(buffer, n);
      }
      os << "\n";
    }

    os << "\n~RK~" << std::endl;
    for (int i=0; i<N; i++) {
      for (int j=0; j<P; j++) {
        sprintf(buffer, "%8.4f\t", m.rk[i*P+j]);
        os << std::string(buffer, n);
      }
      os << "\n";
    }

    os << "\n~G~" << std::endl;
    for (int i=0; i<N; i++) {
      os << m.C[i] << "\t";
    }
    return os << std::endl;
  }

}
