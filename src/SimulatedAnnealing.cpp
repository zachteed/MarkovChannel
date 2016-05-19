#include "SimulatedAnnealing.hpp"

#include <vector>
#include <string>
#include <math.h>

using namespace std;


namespace SimulatedAnnealing
{
  void solve(cost_function cost, MarkovChannel::SAParameter& params)
  {
    int n_chains = params.n_chains();
    int k_max = params.k_max();
    int step = params.step();

    double T0 = params.t0();
    double gamma = params.gamma();

    vector<Model::Model*> models(n_chains, NULL);
    vector<Model::Model*> best(n_chains, NULL);

    vector<double> model_cost(n_chains, 1e12);
    vector<double> best_cost(n_chains, 1e12);

    #pragma omp parallel for
    for (int i=0; i<n_chains; i++) {
      models[i] = new Model::Model(Math::rng_int(3, 10));
      best[i] = new Model::Model(*models[i]);
      model_cost[i] = cost(*models[i]);
      best_cost[i] = model_cost[i];
    }

    Model::Model argmin;
    double min_cost, T=T0;

    for (int i=0; i<k_max; i++) {

      if ( (i > 0) && (i % step == 0) )
        T *= gamma;

      #pragma omp parallel for
      for (int j=0; j<n_chains; j++) {

        Model::Model* neighbor = Model::neighbor(*models[j]);
        double cst = cost(*neighbor);

        if ( Math::rng_uniform() < exp(-(cst-model_cost[j]/T)) ) {
          delete models[j];
          models[j] = neighbor;
          model_cost[j] = cst;
        }

        if ( cst < best_cost[j] ) {
          delete best[j]; best_cost[j] = cst;
          best[j] = new Model::Model(*models[j]);
        }

        if ( Math::rng_uniform() < 0.001 ){
          int idx = Math::rng_int(0, best.size());
          delete models[j];
          models[j] = new Model::Model(*best[idx]);
        }

      }

      for (int j=0; j<n_chains; j++) {
        if ( model_cost[j] < min_cost ) {
          argmin = *models[j];
          min_cost = model_cost[j];
        }
      }

      if ( i % params.display() == 0 ) {
        cout << "Iteration " << i << "\t Loss: " << min_cost << endl;
        cout << argmin << endl;
      }

    }
  }
}
