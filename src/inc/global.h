/*
 * inc/Global.h
 * 
 * Copyright (c) 2011-2013 BGI-Shenzhen <soap at genomics dot org dot cn>. 
 *
 * This file is part of SOAPdenovo-Trans.
 *
 * SOAPdenovo-Trans is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SOAPdenovo-Trans is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SOAPdenovo-Trans.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

int * contig_index_array = NULL;
int scaffNum=0;
int gapNum=1;
boolean fill=0;
int overlaplen = 23;
int inGraph;
long long n_ban;
long long n_solexa = 0;
long long prevNum = 0;
int ins_size_var = 20;
PE_INFO * pes = NULL;
MEM_MANAGER * rv_mem_manager = NULL;
MEM_MANAGER * cn_mem_manager = NULL;
MEM_MANAGER * arc_mem_manager = NULL;
unsigned int num_vt = 0;
unsigned long long new_num_vt=0;
unsigned int ** found_routes = NULL;
unsigned int * so_far = NULL;
int max_n_routes = 10;
int num_trace;
Kmer WORDFILTER;
unsigned int num_ed = 0;
unsigned int num_ctg = 0;
unsigned int num_ed_limit;
unsigned int extraEdgeNum;
EDGE * edge_array = NULL;
VERTEX * vt_array = NULL;
unsigned int * index_array = NULL;
CONTIG * contig_array = NULL;
int lineLen;
int len_bar = 100;
int weakPE = 3;
int fillGap = 0;
boolean globalFlag;
long long arcCounter;
MEM_MANAGER * prearc_mem_manager = NULL;
MEM_MANAGER ** preArc_mem_managers = NULL;
int maxReadLen = 0;
int maxReadLen4all = 0;
int minReadLen = 0;
int maxNameLen = 0;
ARC ** arcLookupTable = NULL;
long long * markersArray = NULL;
boolean deLowKmer = 0;
boolean deLowEdge = 1;
long long newCntCounter;
boolean repsTie = 0;
CONNECT ** cntLookupTable = NULL;
int num_libs = 0;
LIB_INFO * lib_array = NULL;
int libNo = 0;
long long readNumBack;
int gradsCounter;
unsigned int ctg_short = 0;
int thrd_num = 8;
int cvgAvg = 0;
KmerSet ** KmerSets = NULL;
KmerSet ** KmerSetsPatch = NULL;
DARRAY * readSeqInGap = NULL;
DARRAY * gapSeqDarray = NULL;
DARRAY ** darrayBuf;
boolean orig2new;
int maxSteps;
boolean maskRep = 1;
int GLDiff = 50;
int initKmerSetSize = 0;
int RPKM = 0;  //Tang 2013-3-30
long known_genome_size = 0;
int smallKmer = 0;
int deltaKmer = 0;
//mao
CONTIG_PATH *contigpath=NULL;
READ_PATH *readpath=NULL;
int *pool=NULL;
int num_pool=0;
unsigned int* flag_array =NULL;	//lzy 0602
CONTIG_PATHID * contig_path_array=NULL;
PATH_CONTIGID * path_contig_array=NULL;
LOCUS * component=NULL;
int cut_length=48;
int max_num=5;
int max_cnt=0;
int ctg_mask=100;

//mao 9-21
int dk = 5;
int de = 20;
int dd = 5 ;//mao 2011 10 21
int da = 5;
int dA =2;
int delowArc =200;
int SupportNum=2;
int N_kmer=0;//mao 2011-10-13
int read_trace=0;//mao 2011-10-21

#ifdef MER127
Kmer kmerZero={0,0,0,0};
#endif
#ifdef MER63
Kmer kmerZero={0,0};
#endif
#ifdef MER31
Kmer kmerZero={0};
#endif

