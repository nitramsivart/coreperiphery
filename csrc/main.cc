#include "main.h"

const int n = 1000;
const int nodes = 2*n;


double ***messages = new double**[nodes];
double **mprob = new double*[nodes];
double *field = new double[2];

double omega[4]; //= {C_cc/n, C_cp/n, C_cp/n, C_pp/n};
double eta[2] = {.5, .5};

// create the graph
std::vector<std::vector<int> > graph(nodes);

int main(int argc, char* argv[]) {
  // Allocate messages, mprob arrays
  clock_t t;
  for(int i = 0; i < nodes; ++i) {
    messages[i] = new double*[nodes];
    for(int j = 0; j < nodes; ++j)
      messages[i][j] = new double[2];
    mprob[i] = new double[2];
  }
  printf("cycles1: %d\n", (int)(clock() - t));
  t = clock();
  srand(time(NULL));
  double C_cc = 22;
  double C_cp = 20;
  double C_pp = 18;
  for(; C_pp < C_cp; ) {
    C_pp++;
    C_cc--;
    omega[0] = C_cc/n;
    omega[1] = C_cp/n;
    omega[2] = C_cp/n;
    omega[3] = C_pp/n;


    graph.clear();
    graph.resize(nodes);

    for(int i = 0; i < nodes; i++) {
      for(int j = i+1; j < n; j++) {
        if(rand() % n < (i < n ? C_cc : C_cp)) {
          graph[i].push_back(j);
          graph[j].push_back(i);
        }
      }
      for(int j = std::max(n, i+1) ; j < nodes; j++) {
        if(rand() % n < (i < n ? C_cp : C_pp)) {
          graph[i].push_back(j);
          graph[j].push_back(i);
        }
      }
    }
    printf("cycles2: %d\n", (int)(clock() - t));
    t = clock();

    // initialize messages and mprob
    for(int i = 0; i < nodes; i++) {
      //printf("(%d)\n",graph[i].size());
      double a = .60;
      mprob[i][0] = a;
      mprob[i][1] = 1-a;
      for(int j = 0; j < nodes; j++) {
        messages[i][j][0] = (i == j ? 0 : a);
        messages[i][j][1] = (i == j ? 0 : 1-a);
      }
    }
    printf("cycles3: %d\n", (int)(clock() - t));
    t = clock();

    // compute field
    calc_field();
    printf("cycles4: %d\n", (int)(clock() - t));
    t = clock();
    for(int loop = 0; loop < 2; loop++) {
      // compute messages, mprob, and update field
      calc_messages_f();
    }
    printf("cycles5: %d\n", (int)(clock() - t));
    t = clock();

    int correct = 0;
    for(int i = 0; i < nodes; i++) {
      correct += (i < n xor mprob[i][0] < .5 ? 1 : 0);
    }
    printf("C_pp: %d, correct: %d\n", (int)C_pp, (int)correct);
  }

  // De-Allocate memory
  for(int i = 0; i < nodes; ++i) {
    for(int j = 0; j < nodes; ++j)
      delete[] messages[i][j];
    delete[] messages[i];
    delete[] mprob[i];
  }
  delete[] messages;
  delete[] mprob;
  delete[] field;
}

void calc_field() {
  for(int r = 0; r < 2; r++) {
    double prod_total = 0.;
    for(int i = 0; i < nodes; i++) {
      double sum_total = 0.;
      for(int s = 0; s < 2; s++)
        sum_total += (1-omega[2*s + r]) * mprob[i][s];
      prod_total += log(sum_total);
    }
    field[r] = exp(prod_total);
  }
} 

void calc_messages_f() {
  for(int u = 0; u < nodes; u++) {
    std::vector<int> nbrs = graph[u];

    //compute the right hand side product, including v
    double rhs_total[2] = {0,0};
    for(int i_w = 0; i_w < nbrs.size(); ++i_w) {
      int w = graph[u][i_w];
      //if(w == v) continue;
      
      double sum_total_top[2] = {0,0};
      double sum_total_bot[2] = {0,0};
      for(int r = 0; r < 2; r++) {
        for(int s = 0; s < 2; s++) {
          sum_total_top[r] += messages[w][u][s] * omega[2*s+r];
          sum_total_bot[r] += mprob[w][s] * (1-omega[2*s+r]);
        }
        rhs_total[r] += log(sum_total_top[r]) - log(sum_total_bot[r]);
      }
    }

    double new_mprob[2] = {0,0};
    for(int v = 0; v < nodes; v++) {
      if(v == u) continue;

      // Calculate whether adjacent
      bool adj = false;
      for(int index = 0; index < nbrs.size(); ++index) {
        if(nbrs[index] == v) {
          adj = true;
          break;
        }
      }

      double total = 0;
      for(int r = 0; r < 2; r++) {
        double sum_total_top = 0.;
        double sum_total_bot = 0.;
        for(int s = 0; s < 2; s++) {
          double prob = (adj ? omega[2*r+s] : 1 - omega[2*r+s]);
          sum_total_top += messages[v][u][s] * prob;
          sum_total_bot += mprob[v][s] * (1-omega[2*s+r]);
        }

        messages[u][v][r] = field[r] * exp(rhs_total[r] - log(sum_total_top)
                                                  + log(sum_total_bot));
        new_mprob[r] += log(sum_total_top);
        total += messages[u][v][r];
      }
      // Normalize!
      messages[u][v][0] /= total;
      messages[u][v][1] /= total;
    }
    double old_mprob[2] = {mprob[u][0], mprob[u][1]};
    // Normalize!
    new_mprob[0] = exp(new_mprob[0]);
    new_mprob[1] = exp(new_mprob[1]);
    for(int r = 0; r < 2; r++) {
      mprob[u][r] = new_mprob[r] / (new_mprob[0] + new_mprob[1]);
    }
    // Adjust field!
    for(int r = 0; r < 2; r++) {
      double old_sum = old_mprob[0] * (1-omega[2*0+r]) 
                     + old_mprob[1] * (1-omega[2*1+r]);
      double new_sum = mprob[u][0] * (1-omega[2*0+r])
                     + mprob[u][1] * (1-omega[2*1+r]);
      //printf("sum_dif: %f\n", new_sum/old_sum);
      field[r] *= new_sum / old_sum;
    }
  }
}
