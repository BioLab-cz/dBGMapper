/*
 * main.c
 *
 *  Created on: Sep 12, 2019
 *      Author: pluto
 */

#include "analyse_dBG.h"
#include "basic.h"
#include "Binary_Search.h"
#include "read.h"
#include "load_DBG_full.h"
#include "seeding.h"
#include "prealignment.h"
#include "method.h"
#include "print.h"

int main(int argc, char** argv)
{

	char *p_dbg_path;
	char *kmertest,*alignseq;
	uint32_t kmerlen,extlen;
	char dir;
	uint32_t tau;
	uint32_t loop;
	uint32_t read_s;
	for(uint32_t i=1;i<argc;i=i+2)
	{
		if(argv[i][0]=='-'&&argv[i][1]=='r')//
		{
			p_dbg_path=argv[i+1];
		}
		if(argv[i][0]=='-'&&argv[i][1]=='L')//
		{
			kmerlen = atoi(argv[i+1]);
		}
		if(argv[i][0]=='-'&&argv[i][1]=='e')//
		{
			extlen=atoi(argv[i+1]);
		}
		if(argv[i][0]=='-'&&argv[i][1]=='d')//
		{
			dir = *argv[i+1];
		}
		if(argv[i][0]=='-'&&argv[i][1]=='s')//
		{
			kmertest = argv[i+1];
		}
		if(argv[i][0]=='-'&&argv[i][1]=='a')//
		{
			alignseq = argv[i+1];
		}
		if(argv[i][0]=='-'&&argv[i][1]=='t')//
		{
			tau = atoi(argv[i+1]);
		}
		if(argv[i][0]=='-'&&argv[i][1]=='l')//
		{
			loop = atoi(argv[i+1]);
		}
		if(argv[i][0]=='-'&&argv[i][1]=='c')//
		{
			read_s = atoi(argv[i+1]);
		}

	}
//	vector<PH_Node> vPHA;
//	generate_PHNArray(vPHA, tau);
//	cout << loop << "-" << read_s << "index is:" << find_PHNAindex(vPHA, loop, read_s);
//	getchar();

	struct timeval tvs,tve;
//	gettimeofday(&tvs,NULL);
	char nindex[] = "./nindex";
	char rindex[] = "./rindex";
	read_bfile2index(nindex,gnFMidx,0);
	read_bfile2index(rindex,grFMidx,0);
	gPair[0] = make_pair('O',&grFMidx);
	gPair[1] = make_pair('I',&gnFMidx);

	vector<ms_seed> vseed;
	vector<ms_candidate> vcand;

	vector<PH_Node> vPHN;
	generate_PHNArray(vPHN, tau);
	generate_seeds(kmertest, tau, vseed);
	init_candidate(vseed, gnFMidx, vcand);
	candi2candi(vseed,vPHN, vcand, loop);

	cout << "Ee" << endl;
	getchar();



	get_para(&gbit_para, kmerlen);
	char *ref;
	gen_dBG_index(gbit_para, gdBGindex, p_dbg_path,1,&ref);

	gettimeofday(&tve,NULL);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout << "time is: "<< span << endl;
//
	struct TPTnode rootnode;
//	struct seed_extpara ext_set(grFMidx);
//	ext_set.dir = dir;
//	ext_set.orignseq = kmertest;
//	ext_set.alignseq = alignseq;
//	ext_set.tau = tau;
//	ext_set.bit_para = gbit_para;
//	ext_set.sdBGidx = gdBGindex;
//	ext_set.onunipath = false;

//	if(ext_set.dir == 'I')
//	{
//		ext_set.pFMidx = gnFMidx;
//	}
//	if(ext_set.dir == 'O')
//	{
//		ext_set.pFMidx = grFMidx;
//	}

//	init_rootnode(&rootnode, ext_set);
//	ext_treenode(ext_set, &rootnode,extlen+tau);
//	cout << "ext_treenode finished!" << endl;
//	char *extseq = new char[extlen+5]();
//	print_extree(rootnode, extseq);
//	cout << "print_extree done!\n";
	/*
	gettimeofday(&tve,NULL);
	free_FMindex(&nFMidx);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout << "time is: "<< span << endl;
	cout << timetmp << endl;
	cout << "verification cost " << timetmp / span * 100 << "%" << endl;

	printf("%s is over!\n",funname);
	*/
	free_dBGindex(gdBGindex);
	return 0;
}
