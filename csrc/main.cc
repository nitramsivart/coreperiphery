#include "main.h"

const int n = 300;
const int nodes = 2*n;


double *messages = new double[nodes*nodes*2];
double *mprob = new double[nodes*2];
double *field = new double[2];

double omega[4]; //= {C_cc/n, C_cp/n, C_cp/n, C_pp/n};
double eta[2] = {.5, .5};

// create the graph
std::vector<std::vector<int> > graph(nodes);

int main(int argc, char* argv[]) {
  srand(time(NULL));
  double C_cc = 30;
  double C_cp = 20;
  double C_pp = 10;
  time_t begin, end;
  time(&begin);
  for(; C_pp < C_cp; ) {
    C_pp++;
    C_cc--;
    omega[0] = C_cc/n;
    omega[1] = C_cp/n;
    omega[2] = C_cp/n;
    omega[3] = C_pp/n;


    std::vector<std::vector<int> > temp_graph(nodes);
    graph = temp_graph;

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

    // initialize messages and mprob
    for(int i = 0; i < nodes; i++) {
      //printf("(%d)\n",graph[i].size());
      double a = .60;
      mprob[A2_TO_I(i,0)] = a;
      mprob[A2_TO_I(i,1)] = 1-a;
      for(int j = 0; j < nodes; j++) {
        messages[A3_TO_I(i,j,0)] = (i == j ? 0 : a);
        messages[A3_TO_I(i,j,1)] = (i == j ? 0 : 1-a);
      }
    }

    // compute field
    calc_field();
    for(int loop = 0; loop < 4; loop++) {
      // compute messages, mprob, and update field
      calc_messages_f();

      // compute mprob
      // calc_mprob();
      //printf("Loop number: %d\n", loop);
      //for(int i = 0; i < nodes; i+=100) {
      //  printf("%f\n", mprob[A2_TO_I(i,0)]);
      //}
    }

    int correct = 0;
    for(int i = 0; i < nodes; i++) {
      correct += (i < n xor mprob[A2_TO_I(i,0)] < .5 ? 1 : 0);
    }
    printf("C_pp: %d, correct: %d\n", (int)C_pp, (int)correct);
  }
  time(&end);
  printf("Time: %f\n", difftime(end, begin));

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
        sum_total += (1-omega[2*s + r]) * mprob[A2_TO_I(i, s)];
      prod_total += log(sum_total);
    }
    field[r] = exp(prod_total);
  }
}

void calc_messages() {
  for(int u = 0; u < nodes; u++) {
    if(u % 10000 == 0)
      printf("node: %d\n", u);
    for(int v = 0; v < nodes; v++) {
      if(u == v) continue;
      double total = 0;
      for(int r = 0; r < 2; r++) {
        double prod_total = 0.;
        //std::vector<int>::iterator it = graph[u].begin();
        for(int index = 0; index < graph[u].size(); ++index) {
          int w = graph[u][index];
          if(w == v) continue;
          double sum_total_top = 0.;
          double sum_total_bot = 0.;
          for(int s = 0; s < 2; s++) {
            sum_total_top += messages[A3_TO_I(w, u, s)] * omega[2*s+r];
            sum_total_bot += mprob[A2_TO_I(w, s)] * (1-omega[2*s+r]);
          }
          prod_total += log(sum_total_top) - log(sum_total_bot);
        }
        messages[A3_TO_I(u,v,r)] = field[r] * exp(prod_total);
        total += messages[A3_TO_I(u,v,r)];
      }
      // Normalize!
      messages[A3_TO_I(u,v,0)] /= total;
      messages[A3_TO_I(u,v,1)] /= total;
    }
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
          sum_total_top[r] += messages[A3_TO_I(w, u, s)] * omega[2*s+r];
          sum_total_bot[r] += mprob[A2_TO_I(w, s)] * (1-omega[2*s+r]);
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
          sum_total_top += messages[A3_TO_I(v, u, s)] * prob;
          sum_total_bot += mprob[A2_TO_I(v, s)] * (1-omega[2*s+r]);
        }

        messages[A3_TO_I(u,v,r)] = field[r] * exp(rhs_total[r] - log(sum_total_top)
                                                  + log(sum_total_bot));
        new_mprob[r] += log(sum_total_top);
        total += messages[A3_TO_I(u,v,r)];
      }
      // Normalize!
      messages[A3_TO_I(u,v,0)] /= total;
      messages[A3_TO_I(u,v,1)] /= total;
    }
    double old_mprob[2] = {mprob[A2_TO_I(u,0)], mprob[A2_TO_I(u,1)]};
    // Normalize!
    new_mprob[0] = exp(new_mprob[0]);
    new_mprob[1] = exp(new_mprob[1]);
    for(int r = 0; r < 2; r++) {
      mprob[A2_TO_I(u,r)] = new_mprob[r] / (new_mprob[0] + new_mprob[1]);
    }
    // Adjust field!
    for(int r = 0; r < 2; r++) {
      double old_sum = old_mprob[0] * (1-omega[2*0+r]) 
                     + old_mprob[1] * (1-omega[2*1+r]);
      double new_sum = mprob[A2_TO_I(u,0)] * (1-omega[2*0+r])
                     + mprob[A2_TO_I(u,1)] * (1-omega[2*1+r]);
      //printf("sum_dif: %f\n", new_sum/old_sum);
      field[r] *= new_sum / old_sum;
    }
  }
}

void calc_mprob() {
  for(int u = 0; u < nodes; u++) { 
    double total = 0;
    for(int r = 0; r < 2; r++) {
      double prod_total = 0;
      for(int w = 0; w < nodes; w++) {
        if(w == u) continue;
        bool adj = false;
        for(int index = 0; index < graph[u].size(); ++index) {
          if(graph[u][index] == w) {
            adj = true;
            break;
          }
        }
        //if(std::find(graph[u].begin(), graph[u].end(), w) != graph[u].end())
          //adj = true;
        double sum_total = 0;
        for(int s = 0; s < 2; s++) {
          double prob = (adj ? omega[2*r+s] : 1 - omega[2*r+s]);
          sum_total += messages[A3_TO_I(w,u,s)] * prob;
        }
        prod_total += log(sum_total);
      }
      mprob[A2_TO_I(u,r)] = exp(prod_total);
      total += mprob[A2_TO_I(u,r)];
    }
    // Normalize!
    mprob[A2_TO_I(u,0)] /= total;
    mprob[A2_TO_I(u,1)] /= total;
  }
}
