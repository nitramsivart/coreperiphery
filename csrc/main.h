#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <time.h>
#include <unordered_map>
#include <assert.h>

// These are macros to help make a multi-dimensional array
#define A3_TO_I(u,v,r) (u*n*2 + v*2 + r)
#define I_TO_A3_1(i) (i % 2)
#define I_TO_A3_2(i) ((i%(n*2)) / 2)
#define I_TO_A3_3(i) (i / (n*2))

#define A2_TO_I(u,r) (u*2 + r)
#define I_TO_A2_1(i) (i % 2)
#define I_TO_A2_2(i) (i / 2)

// Conditionally invert:
// returns m   if s=0
//         1-m if s=1
#define PAIR(i,j) (i*nodes+j)
#define COND_INV(s,m) (m+s*(1-2*m))

void calc_field();
void calc_messages();
void calc_messages_f();
void calc_mprob();
void calc_omega();
