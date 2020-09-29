/*
 * seeding.cpp
 *
 *  Created on: Sep 18, 2019
 *      Author: bio
 */

//不同的seeding的方法
//seeding的操作是将原read序列拆成包含tau+1个seeds的划分
//目标是拆分出的seed的目标
//1）不能精确匹配的seed数量越多越好
//2）精确匹配的seed的长度越长越好

#include "seeding.h"
#include "Hash.h"
#include "Binary_Search.h"
#include "method.h"
#include "print.h"

std::pair<char, const struct sFMindex *> gPair[2];

void generate_seeds(const char *read, uint32_t tau, vector<ms_seed> & vseed)
{
	uint32_t readlen = strlen(read);
	uint32_t cnt = tau + 1;
	struct ms_seed seedtmp;
	uint32_t start = 0;
	while(cnt)
	{
		uint32_t lentmp = (double)readlen / cnt + 0.5;
		seedtmp.seed_segment = new char[lentmp+1]();
		strncpy(seedtmp.seed_segment,read+start,lentmp);
		seedtmp.start_pos = start;
		seedtmp.end_pos = start + lentmp - 1;
		vseed.push_back(seedtmp);
		start += lentmp;
		readlen -= lentmp;
		--cnt;
	}

}

void free_seed(vector<ms_seed> & vseed)
{
	for(uint32_t i = 0; i != vseed.size(); ++i)
	{
		delete [] vseed[i].seed_segment;
	}
}

void generate_PHNArray(vector<PH_Node> &PHArray, uint8_t tau)
{
	uint32_t level = ceil(log(tau+1) / log(2));
	uint32_t nodenum = pow(2,level+1);
	PHArray.resize(nodenum);
	for(uint32_t i = 0; i < nodenum; ++i)
	{
		PHArray[i].start_seed_id = -1;
		PHArray[i].end_seed_id = -1;
		PHArray[i].tau = -1;
	}
	PHArray[1].start_seed_id = 0;
	PHArray[1].end_seed_id = tau;
	PHArray[1].tau = tau;
	uint32_t cnt = 0;
	for(uint32_t i = 2; i < nodenum; ++i)
	{
		if(PHArray[i/2].tau)
		{
			if(i % 2 == 0)
			{
				if(PHArray[i/2].tau == 1 || PHArray[i/2].tau == 3 || PHArray[i/2].tau == 7)
				{
					PHArray[i].tau = PHArray[i/2].tau / 2;
				}
				else
				{
					PHArray[i].tau = (PHArray[i/2].tau + 1) / 2;
				}
				PHArray[i].start_seed_id = PHArray[i/2].start_seed_id;
				PHArray[i].end_seed_id = PHArray[i].start_seed_id + PHArray[i].tau;
			}
			else
			{
				PHArray[i].tau = PHArray[i/2].tau - PHArray[i-1].tau - 1;
				PHArray[i].start_seed_id = PHArray[i-1].end_seed_id + 1;
				PHArray[i].end_seed_id = PHArray[i].start_seed_id + PHArray[i].tau;
			}
			if(PHArray[i].tau == 0)
			{
				++cnt;
				if(cnt == tau + 1)
				{
					break;
				}
			}
		}

	}
	int p = 1;
	for(uint32_t i = 1; i < nodenum; ++i)
	{
		if(i == pow(2,p))
		{
			cout << endl;
			p++;
		}
		cout << + static_cast<uint8_t>(PHArray[i].tau) << " " << + static_cast<uint8_t>(PHArray[i].start_seed_id) << "-" << \
				+ static_cast<uint8_t>(PHArray[i].end_seed_id) << "\t";
	}
}

uint32_t find_PHNAindex(const vector<PH_Node> &PHArray, uint8_t start_id, uint8_t end_id)
{
	for(int i = 1; i < PHArray.size(); ++i) {	//if tau <= 5
		if(PHArray[i].start_seed_id == start_id && PHArray[i].end_seed_id == end_id) {
			return i;
		}
		if(PHArray[i].start_seed_id == -1 && PHArray[i].end_seed_id == -1 && PHArray[i].tau == -1) {
			break;
		}
	}
	return 0; //not found
}

void calc_edarray(struct TPTnode *pnode, char *seq, uint8_t tau)
{
	pnode->edarry = (uint8_t *)malloc(sizeof(uint8_t)*(2*tau+1));
	for(int i = 0; i < 2*tau+1; i++)
	{
		pnode->edarry[i] = tau+1;
	}
	uint8_t flag;
	int upper = pnode->level - tau;
	int seqlen = strlen(seq);
	int start = 0;
	if(upper > 1)
	{
		start = upper - 1;
	}
	uint8_t left, top ,leftop;
	for(int i = 0; i < 2*tau+1; i++)
	{
		flag = 1;
		if(upper + i <= 0)
		{
//			flagr0 = 0;
			continue;
		}
		if(upper + i == seqlen + 1)
		{
			break;
		}
		if(pnode->c == seq[start++])
		{
			flag = 0;
		}
		left = pnode->p_parent->edarry[i+1] + 1;
		top = pnode->edarry[i-1] + 1;
		leftop = pnode->p_parent->edarry[i] + flag;
		if(!i)
		{
			pnode->edarry[i] = min(leftop,left);
			continue;
		}
		if(i == 2*tau)
		{
			pnode->edarry[i] = min(leftop,top);
			break;
		}
		pnode->edarry[i] = min(leftop,top);
		pnode->edarry[i] = min(left,pnode->edarry[i]);
	}
//	for(int i = 0; i < 2*tau+1; i++)
//	{
//		cout << pnode->edarry[i] << " ";
//	}
//	cout << endl;
}

void calc_edarray_iter(char *seq, char c, uint8_t **p2_parent, uint8_t lev, uint8_t tau)
{
	uint8_t *p_cur = (uint8_t*) malloc(sizeof(uint8_t) * (2 * tau + 1));
	for (uint32_t i = 0; i < 2 * tau + 1; i++) {
		p_cur[i] = tau + 1;
	}
	uint8_t flag;
	int upper = lev - tau;
	int seqlen = strlen(seq);
	int start = 0;
	if (upper > 1) {
		start = upper - 1;
	}
	uint8_t left, top, leftop;
	for (int32_t i = 0; i < 2 * tau + 1; i++) {
		flag = 1;
		if (upper + i <= 0) {
			continue;
		}
		if (upper + i == seqlen + 1) {
			break;
		}
		if (c == seq[start++]) {
			flag = 0;
		}
		left = (*p2_parent)[i + 1] + 1;
		top = p_cur[i - 1] + 1;
		leftop = (*p2_parent)[i] + flag;
		if (!i) {
			p_cur[i] = min(leftop, left);
			continue;
		}
		if (i == 2 * tau) {
			p_cur[i] = min(leftop, top);
			break;
		}
		p_cur[i] = min(leftop, top);
		p_cur[i] = min(left, p_cur[i]);
	}
	free(*p2_parent);
	*p2_parent = p_cur;
}

void calc_kmer_edVec(const struct seed_extpara &ext_para, struct TPTnode *pnode, struct ms_candidate *pcandi)
{
	int lev = 0;
	if(ext_para.dir == 'I')
	{
		for(int i = pcandi->cur_start_pos-1; i >= 0; --i)
		{
			calc_edarray_iter(ext_para.alignseq, *(pcandi->p_path+i), &(pnode->edarry), lev++, ext_para.tau);
		}
	}
	else
	{
		int path_len = strlen(pcandi->p_path);
		for(int i = pcandi->cur_end_pos+1; i < path_len; ++i)
		{
			calc_edarray_iter(ext_para.alignseq, *(pcandi->p_path+i), &(pnode->edarry), lev++, ext_para.tau);
		}
	}
}


void init_rootnode(struct TPTnode *pnode,const struct seed_extpara &ext_para, struct ms_candidate *pcandi)   //mod == 0 edvec 0~tau   mod != 0 edvec calc the cur_seed_id + 1 ~ parent's end_seed_id
{
	pnode->c = 'R';
	pnode->level = 0;
	pnode->offset = 0;
	pnode->edarry = (uint8_t*)malloc(sizeof(uint8_t)*(2*ext_para.tau+1));
	int i = 0;
	pnode->p_parent = nullptr;
	for(i = 0 ; i < 4; i++)
	{
		pnode->p_child[i] = nullptr;
	}
	i = 0;
	for(; i < ext_para.tau+1; i++)
	{
		pnode->edarry[i] = 0;
	}
	int tmp = 0;
	if(pcandi == nullptr)
	{
		for(; i < 2*ext_para.tau+1; i++)
		{
			pnode->edarry[i] = ++tmp;
		}
	}
	else
	{
		calc_kmer_edVec(ext_para, pnode, pcandi);
	}
	char *tmpseq = (char *)malloc(sizeof(char)*(strlen(ext_para.orignseq)+1));
	memset(tmpseq,0,strlen(ext_para.orignseq)+1);
	strcpy(tmpseq,ext_para.orignseq);
	if(ext_para.dir == 'O')
	{
		reverseq(tmpseq);
	}
	pnode->saarry = calc_SArangeSeq(*ext_para.pFMidx,tmpseq);
	free(tmpseq);
}

bool init_childnode(const struct seed_extpara &ext_para, struct TPTnode *pnodec, struct TPTnode *pnodep)  //返回真表示继续扩展
{
	pnodec->p_parent = pnodep;
	pnodec->level = pnodep->level + 1;
	calc_edarray(pnodec, ext_para.alignseq, ext_para.tau);
	pnodec->saarry = calc_SArangeChar(*ext_para.pFMidx,pnodep->saarry,pnodec->c);
	if(pnodec->saarry[1] >= pnodec->saarry[0])
	{
		for(int i = 0; i < 2*ext_para.tau+1; i++)
		{
			if(pnodec->edarry[i] <= ext_para.tau)
			{
				for(uint32_t ii = 0 ; ii < 4; ii++)
				{
					pnodec->p_child[ii] = nullptr;
				}
				return true;
			}
		}
	}
	return false;
}

void ext_treenode(struct seed_extpara &ref_para, struct TPTnode *pnode, uint32_t extlen)
{
	if(extlen != 0)
	{
		char *original_save = ref_para.orignseq;
		bool extflag = false;
		if(pnode->offset > 0 && pnode->offset < strlen(ref_para.orignseq) - gbit_para.kmer1Len/2) //offset > 1?
		{
			if(ref_para.dir == 'I')
			{
				pnode->p_child[0] = (struct TPTnode *)malloc(sizeof(struct TPTnode));
				pnode->p_child[0]->offset = pnode->offset - 1;
				pnode->p_child[0]->c = ref_para.orignseq[pnode->p_child[0]->offset];  //
			}
			else
			{
				pnode->p_child[0] = (struct TPTnode *)malloc(sizeof(struct TPTnode));
				pnode->p_child[0]->offset = pnode->offset + 1;
				pnode->p_child[0]->c = ref_para.orignseq[pnode->offset+gbit_para.kmer1Len/2];
			}
			extflag = init_childnode(ref_para, pnode->p_child[0], pnode);
			if(extflag == true)
			{
				ref_para.onunipath = true;
				ext_treenode(ref_para, pnode->p_child[0], extlen-1);
			}
			else
			{
				free(pnode->p_child[0]);
				pnode->p_child[0] = nullptr;
			}
		}
		else
		{
			uint64_t * const hashvalue_tmp = (uint64_t *)malloc(sizeof(uint64_t) * gbit_para.kmer64Len);
			uint32_t seqlen = gbit_para.kmer1Len/2;
			char *seqe = (char *)malloc(sizeof(char) * (seqlen+1));
			memset(seqe,0,seqlen+1);
			if(pnode->c == 'R')
			{
				strcpy(seqe,ref_para.orignseq);
			}
			else
			{
				if(true == ref_para.onunipath)
				{
					strncpy(seqe,ref_para.orignseq+pnode->offset,seqlen);
				}
				else
				{
					if(ref_para.dir == 'I')
					{
						seqe[0] = pnode->c;
						strncpy(seqe+1,ref_para.orignseq,seqlen-1);
					}
					else
					{
						strncpy(seqe,ref_para.orignseq+1,seqlen-1);
						seqe[seqlen-1] = pnode->c;
					}
				}
			}
			cal_hash_value_directly_256bit(seqe,hashvalue_tmp,gbit_para);
			uint64_t bkmerfindret,ukmerfindret;
			bkmerfindret = Tfind_arrindexN(gdBGindex.p2_bkmer, gdBGindex.p_branchedkmer, \
					hashvalue_tmp,gbit_para.kmer64Len);
			if(bkmerfindret != ULLONG_MAX)
			{
				uint8_t adinfo = gdBGindex.p_branchedkmerad[bkmerfindret];
				if(ref_para.dir == 'I')
				{
					adinfo >>= 4;
				}
				for(uint32_t i = 0; i < 4; i++)
				{
					if(adinfo & 0x01)
					{
						pnode->p_child[i] = (struct TPTnode *)malloc(sizeof(struct TPTnode));
						pnode->p_child[i]->offset = 0;
						switch(i)
						{
							case 0 :
								pnode->p_child[i]->c = 'A';break;
							case 1 :
								pnode->p_child[i]->c = 'C';break;
							case 2 :
								pnode->p_child[i]->c = 'G';break;
							case 3 :
								pnode->p_child[i]->c = 'T';break;
							default:
								pnode->p_child[i]->c = 'E';
						}
						extflag = init_childnode(ref_para, pnode->p_child[i], pnode);
						if(extflag == true)
						{
							ref_para.onunipath = false;
//							char *original_save = ref_para.orignseq;
							ref_para.orignseq = seqe;
							ext_treenode(ref_para, pnode->p_child[i], extlen-1);
//							ref_para.orignseq = original_save;
						}
						else
						{
							free(pnode->p_child[i]);
							pnode->p_child[i] = nullptr;
						}
					}
					adinfo >>= 1;
				}
			}
			ukmerfindret = Tfind_arrindexN(gdBGindex.p2_ukmer, gdBGindex.p_unbranchedkmer, \
					hashvalue_tmp,gbit_para.kmer64Len);
			if(ukmerfindret != ULLONG_MAX)//如果是unipath上的kmer  判断在unipath上的offset 根据是向出度方向还是入度方向扩展来找
			{
				uint32_t ukmerid = gdBGindex.p_unbranchedkmerid[ukmerfindret];
				char *pos = nullptr;
				pos = strstr(gdBGindex.upath_arr[ukmerid],seqe);
				pnode->offset = pos - gdBGindex.upath_arr[ukmerid];
				char ch;
				if(ref_para.dir == 'I')
				{
					ch = *(pos-1);
					pnode->p_child[0] = (struct TPTnode *)malloc(sizeof(struct TPTnode));
					pnode->p_child[0]->offset = pnode->offset - 1;
					pnode->p_child[0]->c = ch;
				}
				else
				{
					ch = *(pos + strlen(seqe));
					pnode->p_child[0] = (struct TPTnode *)malloc(sizeof(struct TPTnode));
					pnode->p_child[0]->offset = pnode->offset + 1;
					pnode->p_child[0]->c = ch;
				}
				extflag = init_childnode(ref_para, pnode->p_child[0], pnode);
				if(extflag == true)
				{
					ref_para.onunipath = true;
					ref_para.orignseq = gdBGindex.upath_arr[ukmerid];
					ext_treenode(ref_para, pnode->p_child[0], extlen-1);

				}
				else
				{
					free(pnode->p_child[0]);
					pnode->p_child[0] = nullptr;
				}
			}
			free(seqe);
			free(hashvalue_tmp);
		}
		ref_para.orignseq = original_save;
	}
}

void init_candidate(const vector<ms_seed> & vseed, const sFMindex &fmindx, vector<ms_candidate> & vcand)
{
	for(uint32_t i = 0; i != vseed.size(); ++i)
	{
		ms_candidate candtmp(i);
		candtmp.cur_end_pos = strlen(vseed[i].seed_segment) - 1;
		candtmp.p_path = new char[candtmp.cur_end_pos + 2]();
		strcpy(candtmp.p_path,vseed[i].seed_segment);
		candtmp.sa_ary[0] = calc_SArangeSeq(fmindx,vseed[i].seed_segment);
		candtmp.sa_ary[1] = calc_SArangeSeq(fmindx,vseed[i].seed_segment);
		vcand.push_back(candtmp);
	}
}

void candi2candi(const vector<ms_seed> & vseed, const vector<PH_Node> &PHArray, vector<ms_candidate> &vcand, uint32_t i)
{
	//erase candidate in vector
	vector<ms_candidate>::iterator ite = vcand.begin();
	struct ms_candidate candidate = ite[i];
	printvec_candidate(vcand);
//	free(candidate.sa_ary[0]);
//	free(candidate.sa_ary[1]);
	vcand.erase(ite + i);
	printvec_candidate(vcand);

	//find candidate in PHTree
	int PHTidx = find_PHNAindex(PHArray, candidate.start_seed_id, candidate.end_seed_id);
	if(!PHTidx) {return;}

	struct seed_extpara ext_set(gPair[PHTidx % 2]);
	struct PH_Node PHParent = PHArray[PHTidx/2];
	int aliseqlen = 0;
	int candpathlen = strlen(candidate.p_path);
	ext_set.tau = PHParent.tau - PHArray[PHTidx].tau - 1;
	ext_set.orignseq = new char[gbit_para.kmer1Len/2+1]();
	if(ext_set.dir == 'I')
	{
		strncpy(ext_set.orignseq, candidate.p_path, gbit_para.kmer1Len/2);
		for(int i = candidate.start_seed_id + 1; i >= PHParent.start_seed_id; --i)
		{
			aliseqlen += strlen(vseed[i].seed_segment);
		}
		ext_set.alignseq = new char[aliseqlen+1]();
		for(int i = candidate.start_seed_id + 1; i >= PHParent.start_seed_id; --i)
		{
			strcat(ext_set.alignseq,vseed[i].seed_segment);
		}
		reverseq(ext_set.alignseq);
	}
	else
	{
		strncpy(ext_set.orignseq, candidate.p_path+candpathlen-gbit_para.kmer1Len/2, gbit_para.kmer1Len/2);
		for(int i = candidate.end_seed_id + 1; i <= PHParent.end_seed_id; ++i)
		{
			aliseqlen += strlen(vseed[i].seed_segment);
		}
		ext_set.alignseq = new char[aliseqlen+1]();
		for(int i = candidate.end_seed_id + 1; i <= PHParent.end_seed_id; ++i)
		{
			strcat(ext_set.alignseq,vseed[i].seed_segment);
		}
	}

	struct TPTnode rootnode;
	if( (ext_set.dir == 'I' && candidate.cur_seed_id == candidate.start_seed_id) || (ext_set.dir == 'O' && candidate.cur_seed_id == candidate.end_seed_id))
	{
		init_rootnode(&rootnode, ext_set);
	}
	else
	{
		init_rootnode(&rootnode, ext_set, &candidate);
	}
	ext_treenode(ext_set, &rootnode,strlen(ext_set.alignseq)+ext_set.tau);
	char *ext_seq = nullptr;
	int ext_len = strlen(ext_seq);
	int n_pathlen = ext_len + candpathlen;
	char *new_path = new char[n_pathlen + 1]();
	if(ext_set.dir == 'I')
	{
		candidate.tau[0] = rangeEd(ext_seq, ext_set.alignseq, ext_len, aliseqlen, max(ext_len, aliseqlen));
		candidate.sa_ary[0] = calc_SArangeSeq(*ext_set.pFMidx, ext_seq);
		reverseq(ext_seq);
		strcat(new_path, ext_seq);
		strcat(new_path,candidate.p_path);
		delete candidate.p_path;
		candidate.p_path = new_path;
		candidate.start_seed_id = PHParent.start_seed_id;
		candidate.cur_start_pos += ext_len;
		candidate.cur_end_pos += ext_len;

	}
	else
	{
		candidate.tau[1] = rangeEd(ext_seq, ext_set.alignseq, ext_len, aliseqlen, max(ext_len, aliseqlen));
		candidate.sa_ary[1] = calc_SArangeSeq(*ext_set.pFMidx, ext_seq);
		strcat(new_path,candidate.p_path);
		strcat(new_path, ext_seq);
		delete candidate.p_path;
		candidate.p_path = new_path;
		candidate.end_seed_id = PHParent.end_seed_id;
	}

}

void destory_extree(struct TPTnode *pnode)
{
	if(pnode)
	{
		for(uint32_t i = 0; i < 4; i++)
		{
			if(pnode->p_child[i] != nullptr)
			{
				destory_extree(pnode->p_child[i]);
				free(pnode->p_child[i]);
				pnode->p_child[i] = nullptr;
			}
		}
		if(pnode->saarry)
		{
			free(pnode->saarry);
			pnode->saarry = nullptr;
		}
		if(pnode->edarry)
		{
			free(pnode->edarry);
			pnode->edarry = nullptr;
		}
	}
}

