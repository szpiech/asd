#ifndef __ASD_DIST_H__
#define __ASD_DIST_H__

#include "asd-data.h"

unsigned int *make_thread_partition(int &num_threads, int ncols);

double proportion_shared(short A, short B);
double proportion_shared2(short A1, short A2, short B1, short B2);
double proportion_shared_weighted(short A, short B, short C, short D, double fA, double fB);

void calc_pw_as_dist(void *work_order);
void calc_pw_as_dist2(void *order);

void calc_grm(void *order);

double gen_rel(short A, short B, double mu);

void calc_weighted_asd(void *order);

#endif