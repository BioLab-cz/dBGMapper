/*
 * mapperStruct.h
 *
 *  Created on: May 9, 2020
 *      Author: pluto
 */

#ifndef PRINT_H_
#define PRINT_H_

#include "basic.h"
#include "seeding.h"
#include "FMindex_ExactMatch.h"

void verification(const vector<ms_seed> & vseed, const vector<ms_candidate> & vcand, vector<ms_result> & vrslt,
				const char *read, const char *ref, uint32_t tau);
void printvec_seed(const vector<ms_seed> & vseed);
void printvec_candidate(const vector<ms_candidate> & vcand);
void printvec_result(const vector<ms_result> & vrslt);
void print_extree(const struct TPTnode &node,char *seq);
void print_specificlen(struct TPTnode node,struct seed_extpara ext_set, uint32_t extlen, char *seq, struct seed_segment *ps_segment);

#endif /* PRINT_H_ */
