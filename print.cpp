/*
 * mapperStruct.cpp
 *
 *  Created on: May 9, 2020
 *      Author: pluto
 */

#include "print.h"

void print_extree(const struct TPTnode &node,char *seq)
{
	seq[strlen(seq)] = node.c;
	if(node.p_child[0] == nullptr && node.p_child[1] == nullptr && node.p_child[2] == nullptr && node.p_child[3] == nullptr)
	{
		printf("%s\n",seq);
	}
	else
	{
		for(uint32_t i = 0; i < 4; i++)
		{
			if(node.p_child[i] != nullptr)
			{
				print_extree(*node.p_child[i],seq);
			}
		}
	}
	seq[strlen(seq)-1] = '\0';
}

void print_specificlen(struct TPTnode node,struct seed_extpara ext_set, uint32_t extlen, char *seq, struct seed_segment *ps_segment)
{
	if(node.p_child[0] == nullptr && node.p_child[1] == nullptr && node.p_child[2] == nullptr && node.p_child[3] == nullptr)
	{
		seq[strlen(seq)] = node.c;
		if(node.level >= extlen - ext_set.tau && node.level <= extlen + ext_set.tau)
		{
			printf("%s\n",seq);
			if(ps_segment != nullptr)
			{
				uint32_t seqlen = strlen(seq);
				uint32_t mergelen = seqlen + gbit_para.kmer1Len / 2;
				char *seedmerge = (char *)malloc(sizeof(char) * mergelen);
				memset(seedmerge,0,mergelen);
				uint32_t pos = 0;
				if(ext_set.dir == 'I')
				{
					for(uint32_t i = seqlen-1; i > 0; --i)
					{
						seedmerge[pos++] = seq[i];
					}
					strcpy(seedmerge+pos,ext_set.orignseq);
				}
				else if(ext_set.dir == 'O')
				{
					strcpy(seedmerge+pos,ext_set.orignseq);
					strcat(seedmerge,seq+1);
				}
				ps_segment->seedseq = seedmerge;
			}
		}
		seq[strlen(seq)-1] = '\0';
	}
	else
	{
		for(uint32_t i = 0; i < 4; i++)
		{
			if(node.p_child[i] != nullptr)
			{
				seq[strlen(seq)] = node.c;
				print_specificlen(*node.p_child[i], ext_set, extlen, seq, ps_segment);
				seq[strlen(seq)-1] = '\0';
			}
		}
	}
}

void printvec_seed(const vector<ms_seed> & vseed)
{
	for(uint32_t i = 0; i != vseed.size(); ++i)
	{
		cout << vseed[i].start_pos << ":" << vseed[i].end_pos << ":" << vseed[i].seed_segment << endl;
	}
}

void printvec_candidate(const vector<ms_candidate> & vcand)
{
	for(uint32_t i = 0; i != vcand.size(); ++i)
	{
		cout << (uint32_t)vcand[i].cur_seed_id << " ";
	}
	cout << endl;
}

void printvec_result(const vector<ms_result> & vrslt)
{
	vector<ms_result>::const_iterator ite = vrslt.begin();
	for(; ite != vrslt.end(); ++ite)
	{
		cout << ite->ref_start << "~" << ite->ref_end << " ed:"<< ite->ed << endl;
	}
}

void mapper()
{
}
