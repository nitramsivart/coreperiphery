#include "main.h"

const int n = 50000;
const int nodes = 2*n;


std::unordered_map<long, double> messages;
double *mprob = new double[nodes];
double *field = new double[2];

double omega[4]; //= {C_cc/n, C_cp/n, C_cp/n, C_pp/n};
double eta[2]; //= {.5, .5};

// create the graph
std::vector<std::vector<int> > graph(nodes);
int main(int argc, char* argv[]) {
  // Allocate messages, mprob arrays
  srand(time(NULL));
  double C_cc = 36;
  double C_cp = 32;
  double C_pp = 38;
  for(int loop = 0; loop < 2; loop++) {
    //C_pp+=3;
    C_cc+=1;
    //C_cp-=3;
    omega[0] = C_cc/n;
    omega[1] = C_cp/n;
    omega[2] = C_cp/n;
    omega[3] = C_pp/n;
    eta[0] = .5;
    eta[1] = .5;

    clock_t t = clock();
    messages.clear();
    graph.clear();
    graph.resize(nodes);

    for(int i = 0; i < nodes; i++) {
      for(int j = i+1; j < n; j++) {
        if(rand() % n < (i < n ? C_cc : C_cp)) {
          graph[i].push_back(j);
          graph[j].push_back(i);
        }
      }
      for(int j = std::max(n, i+1); j < nodes; j++) {
        if(rand() % n < (i < n ? C_cp : C_pp)) {
          graph[i].push_back(j);
          graph[j].push_back(i);
        }
      }
    }

    // initialize messages and mprob
    double a = .51;
    for(int i = 0; i < nodes; i++) {
      mprob[i] = a;
      for(int j = 0; j < graph[i].size(); j++) {
        messages[PAIR(i,graph[i][j])]=a;
        messages[PAIR(graph[i][j],i)]=a;
      }
      
    }

    // compute field
    calc_field();
    //printf("field: %f, %f\n", field[0], field[1]);
    double total = 0;
    //for(int outerloop = 0; outerloop < 1; outerloop++) {
      for(int innerloop = 0; innerloop < 10; innerloop++) {
        // compute messages, mprob, and update field
        calc_messages_f();
        double new_total = 0;
        for(int i = 0; i < nodes; i+=10000) {
          new_total += mprob[i];
          printf("mprob %f\n", mprob[i]);
        }
        //if(std::abs(new_total - total) < .0000001) {
        if(std::abs(new_total - .5*(nodes/10000)) < .0000001) {
          printf("breaking: %f,%f,%f\n", new_total, total, std::abs(new_total-total));
          break;
        }
        else {
          printf("diff: %f\n", new_total - total);
          total = new_total;
        }
      }
      //calc_omega();
    //}
    //printf("cycles: %d\n", (int)(clock() - t));

    int correct = 0;
    for(int i = 0; i < nodes; i++)
      correct += (i < n xor mprob[i] < .5 ? 1 : 0);
    
    printf("C_pp: %d, correct: %d\n", (int)C_pp, correct);
  }


  // De-Allocate memory
  delete[] mprob;
  delete[] field;
}

void calc_field() {
  for(int r = 0; r < 2; r++) {
    double prod_total = 0.;
    for(int i = 0; i < nodes; i++) {
      double sum_total = 0.;
      for(int s = 0; s < 2; s++)
        sum_total += (1-omega[2*s + r]) * COND_INV(s,mprob[i]);
      prod_total += log(sum_total);
    }
    field[r] = prod_total;
  }
} 

void calc_messages_f() {
  for(int u = 0; u < nodes; u++) {
    std::vector<int> nbrs = graph[u];

    //compute the right hand side product, including v
    double rhs_total[2] = {0,0};
    for(int i_w = 0; i_w < nbrs.size(); ++i_w) {
      int w = nbrs[i_w];
      
      double sum_total_top[2] = {0,0};
      double sum_total_bot[2] = {0,0};
      //assert(COND_INV(0,messages[PAIR(w,u)]) * omega[2*0+0] + COND_INV(1,messages[PAIR(w,u)]) * omega[2*1+0] ==
             //COND_INV(0,messages[PAIR(w,u)]) * omega[2*0+1] + COND_INV(1,messages[PAIR(w,u)]) * omega[2*1+1]);
             //COND_INV(s,messages[PAIR(w,u)]) * omega[2*s+1]);
      //assert(COND_INV(0,mprob[w]) * (1-omega[2*0+0]) + COND_INV(1,mprob[w]) * (1-omega[2*1+0]) ==
             //COND_INV(0,mprob[w]) * (1-omega[2*0+1]) + COND_INV(1,mprob[w]) * (1-omega[2*1+1]));
             //COND_INV(s,mprob[w]) * (1-omega[2*s+1]));
      double assert_totals[2] = {0,0};
      for(int r = 0; r < 2; r++) {
        for(int s = 0; s < 2; s++) {
          sum_total_top[r] += COND_INV(s,messages[PAIR(w,u)]) * omega[2*s+r];
          assert_totals[r] += COND_INV(s,messages[PAIR(w,u)]) * omega[2*s+r];
          sum_total_bot[r] += COND_INV(s,mprob[w]) * (1-omega[2*s+r]);
        }
        rhs_total[r] += log(sum_total_top[r]) - log(sum_total_bot[r]);
      }
      /*
      if(assert_totals[0] != assert_totals[1]) {
        printf("assert_totals: %f, %f\n", assert_totals[0], assert_totals[1]);
      }

      if(messages[PAIR(w,u)] != .5)
        printf("unequal messages: %f, %d, %d\n", messages[PAIR(w,u)], w, u);
      assert(messages[PAIR(w,u)] == .5);
      assert(COND_INV(0,messages[PAIR(w,u)]) == COND_INV(1,messages[PAIR(w,u)]));
      assert(omega[0+0] == omega[2+1]);
      assert(omega[2+0] == omega[0+1]);
      assert(assert_totals[0] == assert_totals[1]);
      */
    }

    //printf("rhs_total: %f, %f\n", rhs_total[0], rhs_total[1]);

    // compute messages for all neighbors
    for(int i_v = 0; i_v < nbrs.size(); i_v++) {
      int v = nbrs[i_v];

      double tmessages[2] = {0,0};
      for(int r = 0; r < 2; r++) {
        double sum_total_top = 0;
        double sum_total_bot = 0;
        for(int s = 0; s < 2; s++) {
          /*
          assert(COND_INV(s,messages[PAIR(v,u)]) * omega[2*s+0] ==
                 COND_INV(s,messages[PAIR(v,u)]) * omega[2*s+1]);
          assert(COND_INV(s,mprob[v]) * (1-omega[2*s+0]) ==
                 COND_INV(s,mprob[v]) * (1-omega[2*s+1]));
                 */
          sum_total_top += COND_INV(s,messages[PAIR(v,u)]) * omega[2*s+r];
          sum_total_bot += COND_INV(s,mprob[v]) * (1-omega[2*s+r]);
        }

        tmessages[r] = eta[r] * exp(field[r] - field[1] + rhs_total[r] - rhs_total[1]
                                    - log(sum_total_top) + log(sum_total_bot));
      }
      // Normalize messages!
      messages[PAIR(u,v)] = tmessages[0] / (tmessages[0]+tmessages[1]);
      //assert(messages[PAIR(u,v)] == .5);
      //printf("tmessages: %f, %f, %f\n", tmessages[0], tmessages[1], messages[PAIR(u,v)]);
      if(tmessages[0] == 0 && tmessages[1] == 0) {
        printf("both zero, %f, %f\n", rhs_total[0], rhs_total[1]);
        messages[PAIR(u,v)] = (field[0] + rhs_total[0] > field[1] + rhs_total[1]);
        //assert(messages[PAIR(u,v)] == .5);
      }
    }
    // Normalize mprob!
    double old_mprob = mprob[u];
    double mprob0 = eta[0] * exp(field[0]-field[1]+rhs_total[0]-rhs_total[1]);
    double mprob1 = eta[1] * exp(0);
    mprob[u] = mprob0 / (mprob0 + mprob1);
    if(mprob0 == 0 && mprob0 == mprob1) {
      printf("error handling, %f %f\n", rhs_total[0], rhs_total[1]);
      mprob[u] = (field[0] + rhs_total[0] > field[1] + rhs_total[1]);
    }

    //printf("mprobs: %f,%f,%f,%f,%d,%d,%f\n", field[0],rhs_total[0],field[1],rhs_total[1],mprob0>mprob1, mprob1>mprob0,mprob[u]);
    //if(isnan(field[0]))
    //  exit(1);

    // Adjust field!
    for(int r = 0; r < 2; r++) {
      double old_sum = old_mprob * (1-omega[2*0+r]) 
                     + (1-old_mprob) * (1-omega[2*1+r]);
      double new_sum = mprob[u] * (1-omega[2*0+r])
                     + (1-mprob[u]) * (1-omega[2*1+r]);
      //assert(old_sum == new_sum);
      //printf("sum_dif: %f\n", new_sum/old_sum);
      field[r] += log(new_sum) - log(old_sum);
    }
  }
}

void calc_omega() {
  double omega_top[4] = {0,0,0,0};
  eta[0] = 0;
  eta[1] = 0;
  //loop over graph
  for(int u = 0; u < nodes; u++) {
    eta[0] += mprob[u];
    eta[1] += (1-mprob[u]);
    std::vector<int> nbrs = graph[u];

    for(int i_v = 0; i_v < nbrs.size(); i_v++) {
      int v = nbrs[i_v];
      double total = 0;
      double joint_est[4] = {0,0,0,0};
      for(int r = 0; r < 2; r++) {
        for(int s = 0; s < 2; s++) {
          joint_est[2*s+r] = COND_INV(r,messages[PAIR(u,v)]) * 
                             COND_INV(s,messages[PAIR(v,u)]) * omega[2*s+r];
          total += joint_est[2*s+r];
        }
      }
      for(int r = 0; r < 2; r++) {
        for(int s = 0; s < 2; s++) {
          assert(total!=0);
          omega_top[2*s+r] += joint_est[2*s+r] / total;
        }
      }
    }
  }
  // Actually update omega. Eta has not been normalized
  for(int r = 0; r < 2; r++) {
    assert(eta[r] != 0);
    for(int s = 0; s < 2; s++) {
      omega[2*s+r] = omega_top[2*s+r] / (eta[r] * eta[s]);
    }
  }
  // Finalize eta
  eta[0] = eta[0] / nodes;
  eta[1] = 1 - eta[0];
}
