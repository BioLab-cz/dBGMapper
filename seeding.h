/*
 * seeding.h
 *
 *  Created on: Sep 18, 2019
 *      Author: bio
 */

#ifndef SEEDING_H_
#define SEEDING_H_

#include "basic.h"
#include "FMindex_ExactMatch.h"
#include "load_DBG_full.h"
#include <utility>

extern std::pair<char, const struct sFMindex *> gPair[2];

struct ms_seed{
	uint8_t start_pos;
	uint8_t end_pos;
	char *seed_segment;
};

struct candiPath{
	char *p_candidatepath;
	uint8_t start_seed_id;
	uint8_t end_seed_id;
	uint8_t *sa_ary;

};

struct PH_Node{
	uint8_t start_seed_id;
	uint8_t end_seed_id;
	uint8_t tau;
};

struct ms_candidate
{
	uint8_t start_seed_id;
	uint8_t cur_seed_id;
	uint8_t end_seed_id;
	uint8_t cur_start_pos;
	uint8_t cur_end_pos;
	uint8_t tau[2];
	char *p_path;
	uint32_t *sa_ary[2];
	ms_candidate(uint32_t p) : start_seed_id(p), cur_seed_id(p), end_seed_id(p),cur_start_pos(0), cur_end_pos(0), tau{0,0}, p_path(nullptr){}
};

struct ms_result
{
	uint32_t ref_start;
	uint32_t ref_end;
	uint32_t ed;
};

struct seed_extpara
{
	bool onunipath;
	char dir;
	uint8_t tau;
	char *orignseq;
	char *alignseq;
	const struct sFMindex *pFMidx;
	seed_extpara(std::pair<char, const struct sFMindex *> &qpr) : dir(qpr.first), pFMidx(qpr.second), onunipath(false) {}
};

struct TPTnode
{
	char c;
	uint8_t level;
	uint16_t offset;			//假设unipath最大长度不超过255   统计如果大于255  改为uint16_t
	struct TPTnode * p_parent;
	struct TPTnode * p_child[4];
	uint8_t * edarry;
	uint32_t * saarry;
};

struct seed_segment
{
	char * seedseq;
	uint32_t * saarry;
//	uint32_t *
};

void generate_PHNArray(vector<PH_Node> &PHArray, uint8_t tau);
uint32_t find_PHNAindex(const vector<PH_Node> &PHArray, uint8_t start_id, uint8_t end_id);
void generate_seeds(const char *read, uint32_t tau, vector<ms_seed> & vseed);
void free_seed(vector<ms_seed> & vseed);
void init_candidate(const vector<ms_seed> & vseed, const sFMindex &fmindx, vector<ms_candidate> & vcand);
void candi2candi(const vector<ms_seed> & vseed, const vector<PH_Node> &PHArray, vector<ms_candidate> &vcand, uint32_t i);
void init_rootnode(struct TPTnode *pnode,const struct seed_extpara &ext_para, struct ms_candidate *pcandi = nullptr);
void ext_treenode(struct seed_extpara &ref_para, struct TPTnode *pnode, uint32_t extlen);
void destory_extree(struct TPTnode *pnode);

#endif /* SEEDING_H_ */
