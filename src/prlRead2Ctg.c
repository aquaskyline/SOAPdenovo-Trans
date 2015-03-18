/*
 * prlRead2Ctg.c
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

#include "stdinc.h"
#include "newhash.h"
#include "extfunc.h"
#include "extvab.h"

static long long readsInGap = 0;

static int buffer_size = 100000000;
static long long readCounter;
static long long mapCounter;
static int ALIGNLEN = 0;

//buffer related varibles for chop kmer
static int read_c;
static char **rcSeq;
static char **seqBuffer;
static int *lenBuffer;
static unsigned int *ctgIdArray;
static int *posArray;
static char *orienArray;
static char *footprint;		// flag indicates whether the read shoulld leave markers on contigs

// kmer related variables
static int kmer_c;
static Kmer *kmerBuffer;
static ubyte8 *hashBanBuffer;
static kmer_t **nodeBuffer;
static boolean *smallerBuffer;
static unsigned int *indexArray;
static int *insSizeArray; //20120104

//contig -> read
READSET ** ctg2read;
unsigned int *countContig;

static int *deletion;
static void parse1read (int t);
static void threadRoutine (void *thrdID);
static void searchKmer (int t, KmerSet * kset);
static void chopKmer4read (int t, int threadID);
static void thread_wait (pthread_t * threads);

static void creatThrds (pthread_t * threads, PARAMETER * paras)
{
	unsigned char i;
	int temp;

	for (i = 0; i < thrd_num; i++)
	{
		//printf("to create %dth thread\n",(*(char *)&(threadID[i])));
		if ((temp = pthread_create (&threads[i], NULL, (void *) threadRoutine, &(paras[i]))) != 0)
		{
			printf ("create threads failed\n");
			exit (1);
		}
	}

	printf ("%d thread created in prlRead2Ctg\n", thrd_num);
}
// 2->1->3
static void threadRoutine (void *para)
{
	PARAMETER *prm;
	int i, t;
	unsigned char id;

	prm = (PARAMETER *) para;
	id = prm->threadID;

	//printf("%dth thread with task %d, hash_table %p\n",id,prm.task,prm.hash_table);
	while (1)
	{
		if (*(prm->selfSignal) == 1)
		{
			for (i = 0; i < kmer_c; i++)
			{
				//if((hashBanBuffer[i]&taskMask)!=prm.threadID)
				if ((hashBanBuffer[i] % thrd_num) != id)
					continue;
				searchKmer (i, KmerSets[id]);
			}
			*(prm->selfSignal) = 0;
		}
		else if (*(prm->selfSignal) == 2)
		{
			for (i = 0; i < read_c; i++)
			{
				if (i % thrd_num != id)
					continue;
				chopKmer4read (i, id + 1);
			}
			*(prm->selfSignal) = 0;
		}
		else if (*(prm->selfSignal) == 3)
		{
			// parse reads
			for (t = 0; t < read_c; t++)
			{
				if (t % thrd_num != id)
					continue;
				parse1read (t);
			}
			*(prm->selfSignal) = 0;
		}
		else if (*(prm->selfSignal) == 5)
		{
			*(prm->selfSignal) = 0;
			break;
		}
		usleep (1);
	}
}

static void chopKmer4read (int t, int threadID)
{
	int len_seq = lenBuffer[t];

	if (len_seq < overlaplen + 1)
	{
		return;
	}

	char *src_seq = seqBuffer[t];
	char *bal_seq = rcSeq[threadID];
	int j, bal_j;
	ubyte8 hash_ban, bal_hash_ban;
	Kmer word, bal_word;
	int index;

	word = kmerZero;
	for (index = 0; index < overlaplen; index++)
	{
#ifdef MER127
		word = KmerLeftBitMoveBy2 (word);
		word.low2 |= src_seq[index];
#endif
#ifdef MER63
		word = KmerLeftBitMoveBy2 (word);
		word.low |= src_seq[index];
#endif
#ifdef MER31
		word = KmerLeftBitMoveBy2 (word);
		word += src_seq[index];
#endif
	}

	reverseComplementSeq (src_seq, len_seq, bal_seq);
	// complementary node
	bal_word = reverseComplement (word, overlaplen);
	bal_j = len_seq - 0 - overlaplen;	//  0;
	index = indexArray[t];

	if (KmerSmaller (word, bal_word))
	{
		hash_ban = hash_kmer (word);
		kmerBuffer[index] = word;
		smallerBuffer[index] = 1;
		hashBanBuffer[index++] = hash_ban;
	}
	else
	{
		bal_hash_ban = hash_kmer (bal_word);
		kmerBuffer[index] = bal_word;
		smallerBuffer[index] = 0;
		hashBanBuffer[index++] = bal_hash_ban;
	}

	//printf("%dth: %p with %p\n",kmer_c-1,bal_word,bal_hash_ban);
	for (j = 1; j <= len_seq - overlaplen; j++)
	{
		word = nextKmer (word, src_seq[j - 1 + overlaplen]);
		bal_j = len_seq - j - overlaplen;	//  j;
		bal_word = prevKmer (bal_word, bal_seq[bal_j]);

		if (KmerSmaller (word, bal_word))
		{
			hash_ban = hash_kmer (word);
			kmerBuffer[index] = word;
			smallerBuffer[index] = 1;
			hashBanBuffer[index++] = hash_ban;
			//printf("%dth: %p with %p\n",kmer_c-1,word,hashBanBuffer[kmer_c-1]);
		}
		else
		{
			// complementary node
			bal_hash_ban = hash_kmer (bal_word);
			kmerBuffer[index] = bal_word;
			smallerBuffer[index] = 0;
			hashBanBuffer[index++] = bal_hash_ban;
			//printf("%dth: %p with %p\n",kmer_c-1,bal_word,hashBanBuffer[kmer_c-1]);
		}
	}
}

//splay for one kmer in buffer and save the node to nodeBuffer
static void searchKmer (int t, KmerSet * kset)
{
	kmer_t *node;
	boolean found = search_kmerset (kset, kmerBuffer[t], &node);

	if (found && !node->deleted)
	{
		nodeBuffer[t] = node;
	}
	else
	{
		nodeBuffer[t] = NULL;
	}
}

static void parse1read (int t)
{
	unsigned int j, i, s;
	unsigned int contigID;
	int counter = 0, counter2 = 0;
	unsigned int ctgLen, pos;
	kmer_t *node;
	boolean isSmaller;
	int flag, maxOcc = 0;
	kmer_t *maxNode = NULL;
	int alldgnLen = lenBuffer[t] > ALIGNLEN ? ALIGNLEN : lenBuffer[t];
	int multi = alldgnLen - overlaplen + 1 < 5 ? 5 : alldgnLen - overlaplen + 1;
	unsigned int start, finish;
	int count_Contig=0;
	unsigned int  pos_temp[20];
	int  alignLen[20];

	footprint[t] = 0;
	start = indexArray[t];
	finish = indexArray[t + 1];

	if (finish == start)	//too short
	{
		ctgIdArray[t] = 0;
		return;
	}
	for (j = start; j < finish; j++)
	{
		node = nodeBuffer[j];
		if (!node)	//same as previous
			continue;
		flag = 1;
		for (s = j + 1; s < finish; s++)
		{
			if (!nodeBuffer[s])
				continue;
			if (nodeBuffer[s]->l_links == node->l_links)
			{
				flag++;
				nodeBuffer[s] = NULL;
			}
		}

		if (flag >= 2)//mao 127kmer->31kmer
			counter2++;
		if (flag >= multi)
			counter++;
		else
			continue;
		pos_temp[count_Contig]=j;
		alignLen[count_Contig]=flag;
		count_Contig++;
		
		if (flag > maxOcc)
		{
			pos = j;
			maxOcc = flag;
			maxNode = node;
		}
	}

	if (!counter)		//no match
	{
		ctgIdArray[t] = 0;
		return;
	}

	if (counter2 > 1)
	{
		footprint[t] = 1;
	}			//use as a flag

	j = pos;
	i = pos - start + 1;
	node = nodeBuffer[j];
	isSmaller = smallerBuffer[j];
	contigID = node->l_links;
	ctgLen = contig_array[contigID].length;
	pos = node->r_links;

	if (node->twin == isSmaller)
	{
		orienArray[t] = '-';
		ctgIdArray[t] = getTwinCtg (contigID);
		posArray[t] = ctgLen - pos - overlaplen - i + 1;
	}
	else
	{
		orienArray[t] = '+';
		ctgIdArray[t] = contigID;
		posArray[t] = pos - i + 1;
	}
	int m,n;
	ctg2read[t]=(READSET*)ckalloc(count_Contig * sizeof(READSET));
	countContig[t]=count_Contig;
	for(m=0;m<count_Contig;m++)
	{
		j=pos_temp[m];
		i=pos_temp[m]-start+1;
		node=nodeBuffer[j];
		isSmaller=smallerBuffer[j];
		contigID=node->l_links;
		ctgLen=contig_array[contigID].length;
		pos = node->r_links;
		if(node->twin == isSmaller)
		{
			ctg2read[t][m].orien='-';
			ctg2read[t][m].contigID=getTwinCtg(contigID);
			ctg2read[t][m].contigOffset=ctgLen - pos - overlaplen;
			ctg2read[t][m].readOffset=i;
			ctg2read[t][m].alignLength=alignLen[m];
		}
		else
		{
			ctg2read[t][m].orien='+';
			ctg2read[t][m].contigID=contigID;
			ctg2read[t][m].contigOffset=pos;
			ctg2read[t][m].readOffset=i;
			ctg2read[t][m].alignLength=alignLen[m];
		}
	}
}

static void sendWorkSignal (unsigned char SIG, unsigned char *thrdSignals)
{
	int t;

	for (t = 0; t < thrd_num; t++)
	{
		thrdSignals[t + 1] = SIG;
	}

	while (1)
	{
		usleep (10);

		for (t = 0; t < thrd_num; t++)
			if (thrdSignals[t + 1])
			{
				break;
			}

		if (t == thrd_num)
		{
			break;
		}
	}
}

static void locate1read (int t)
{
	int i, j, start, finish;
	kmer_t *node;
	unsigned int contigID;
	int pos, ctgLen;
	boolean isSmaller;

	start = indexArray[t];
	finish = indexArray[t + 1];

	for (j = start; j < finish; j++)
	{
		node = nodeBuffer[j];

		if (!node)	//same as previous
		{
			continue;
		}

		i = j - start + 1;
		isSmaller = smallerBuffer[j];
		contigID = node->l_links;
		ctgLen = contig_array[contigID].length;
		pos = node->r_links;

		if (node->twin == isSmaller)
		{
			ctgIdArray[t] = getTwinCtg (contigID);
			posArray[t] = ctgLen - pos - overlaplen - i + 1;
		}
		else
		{
			ctgIdArray[t] = contigID;
			posArray[t] = pos - i + 1;
		}
	}
}

//static void output1read (int t, FILE * outfp)
static void output1read (int t, FILE * outfp, gzFile *outfp2,char orien, int dhflag)
{
	int len = lenBuffer[t];
	int index;

	readsInGap++;

	for (index = 0; index < len; index++)
	{
		writeChar2tightString (seqBuffer[t][index], rcSeq[1], index);
	}

	fwrite (&len, sizeof (int), 1, outfp);
	fwrite (&ctgIdArray[t], sizeof (int), 1, outfp);
	fwrite (&posArray[t], sizeof (int), 1, outfp);
	fwrite (rcSeq[1], sizeof (char), len / 4 + 1, outfp);

	if(fill && insSizeArray[t] < 2000 && len > 0){ 
		gzprintf(outfp2,">%d\t%d\t%d\t%c\t%d\t%d\n",len,ctgIdArray[t],posArray[t],orien,insSizeArray[t],dhflag);
		for(index=0;index<len;index++)
			gzprintf(outfp2,"%c",int2base((int)seqBuffer[t][index]));
		gzprintf(outfp2,"\n");
	}
}

static void getReadIngap (int t, int insSize, FILE * outfp,gzFile * outfp5,boolean readOne)
{
	int read1, read2;
	char orientation;

	if (readOne)
	{
		read1 = t;
		read2 = t + 1;
		if(orienArray[read2]=='+'){
			orientation = '-';
		}else{
			orientation = '+';
		}	
		ctgIdArray[read1] = ctgIdArray[read2];
		posArray[read1] = posArray[read2] + insSize - lenBuffer[read1];	//   --> R2       <-- R1
		output1read (read1, outfp,outfp5,orientation, 1);
	}
	else
	{
		read2 = t;
		read1 = t - 1;
		if(orienArray[read1]=='+'){
			orientation = '-';
		}else{
			orientation = '+';
		}
		ctgIdArray[read2] = ctgIdArray[read1];
		posArray[read2] = posArray[read1] + insSize - lenBuffer[read2];	// --> R1     <-- R2
		output1read (read2,outfp,outfp5,orientation, 1);
	}
}

static void recordLongRead (FILE * outfp)
{
	int t;

	for (t = 0; t < read_c; t++)
	{
		readCounter++;

		if (footprint[t])
		{
//			output1read (t, outfp);
		}
	}
}
static void getPEreadOnContig(int t, gzFile* outfp) //by cyx and lbh.
{
	int len1,len2,index;
	char orien1,orien2;
	len1 = lenBuffer[t-1];
	len2 = lenBuffer[t];
	orien1 = orienArray[t-1];
	orien2 = orienArray[t];

	if(insSizeArray[t]<2000&&insSizeArray[t]==insSizeArray[t-1]){
		gzwrite (outfp, &len1, sizeof (int));
		gzwrite (outfp, &ctgIdArray[t-1], sizeof (int));
		gzwrite (outfp, &posArray[t-1], sizeof (int));
		gzwrite (outfp, &orien1, sizeof (char));
		gzwrite (outfp, &insSizeArray[t-1], sizeof (int));
		
		for (index = 0; index < len1; index++)
		{
			writeChar2tightString (seqBuffer[t-1][index], rcSeq[1], index);
		}
          	gzwrite (outfp, rcSeq[1], (unsigned)(len1 / 4 + 1));

		gzwrite (outfp, &len2, sizeof (int));
		gzwrite (outfp, &ctgIdArray[t], sizeof (int));
		gzwrite (outfp, &posArray[t], sizeof (int));
		gzwrite (outfp, &orien2, sizeof (char));
		gzwrite (outfp, &insSizeArray[t], sizeof (int));
		
		for (index = 0; index < len2; index++)
		{
			writeChar2tightString (seqBuffer[t][index], rcSeq[1], index);
		}
          	gzwrite (outfp, rcSeq[1], (unsigned)(len2 / 4 + 1));
	}
}

static void recordAlldgn (FILE * outfp, int insSize, FILE * outfp2,FILE * outfp3,FILE * outfp4,gzFile * outfp5,gzFile * outfp6)
{
	int t, ctgId;
	boolean rd1gap, rd2gap;
	int m;
	char orientation;
	
	for (t = 0; t < read_c; t++)
	{
		readCounter++;
		rd1gap = rd2gap = 0;
		ctgId = ctgIdArray[t];

		if (outfp2 && t % 2 == 1)	//make sure this is read2 in a pair
		{
			if (ctgIdArray[t] < 1 && ctgIdArray[t - 1] > 0)
			{
				getReadIngap (t, insSizeArray[t], outfp2,outfp5, 0);	//read 2 in gap
				rd2gap = 1;
			}
			else if (ctgIdArray[t] > 0 && ctgIdArray[t - 1] < 1)
			{
				getReadIngap (t - 1,  insSizeArray[t-1], outfp2,outfp5, 1);	//read 1 in gap
				rd1gap = 1;
			}
			else if(ctgIdArray[t]>0&&ctgIdArray[t-1]>0){ //PE read on contig
				if(fill)
					getPEreadOnContig(t,outfp6);
			}
		}
		if (ctgId < 1)
			continue;
		mapCounter++;
//		fprintf (outfp, "%lld\t%u\t%d\t%c\n", readCounter, ctgIdArray[t], posArray[t], orienArray[t]);
		
		if(readCounter% 2 == 1)
			fprintf (outfp, "%lld\t%u\t%d\t%c\n", readCounter, ctg2read[t][countContig[t]-1].contigID, ctg2read[t][countContig[t]-1].contigOffset-ctg2read[t][countContig[t]-1].readOffset+1, ctg2read[t][countContig[t]-1].orien);
		else
			fprintf (outfp, "%lld\t%u\t%d\t%c\n", readCounter, ctg2read[t][0].contigID, ctg2read[t][0].contigOffset-ctg2read[t][0].readOffset+1, ctg2read[t][0].orien);
			
		for(m=0;m<countContig[t];m++)
		{
			if(ctg2read[t][m].alignLength >=5)
				fprintf (outfp3, "%lld\t%u\t%d\t%c\n", readCounter, ctg2read[t][m].contigID, ctg2read[t][m].readOffset - ctg2read[t][m].contigOffset, ctg2read[t][m].orien);
			if(read_trace && (ctg2read[t][m].alignLength >=5))
			{
				if(ctg2read[t][m].orien == '+')
					fprintf (outfp4,"%lld\t%d\t%llu\t%d\t%d\t%c\n",readCounter,ctg2read[t][m].readOffset-1,ctg2read[t][m].contigID,
						ctg2read[t][m].contigOffset,ctg2read[t][m].alignLength+overlaplen-1,ctg2read[t][m].orien);
				else
					fprintf (outfp4,"%lld\t%d\t%llu\t%d\t%d\t%c\n",readCounter,ctg2read[t][m].readOffset-1,getTwinCtg(ctg2read[t][m].contigID),
						contig_array[ctg2read[t][m].contigID].length-ctg2read[t][m].contigOffset-(ctg2read[t][m].alignLength+overlaplen-1),ctg2read[t][m].alignLength+overlaplen-1,ctg2read[t][m].orien);
			}
		}
		free(ctg2read[t]);
		if (t % 2 == 0)
			continue;


		// reads are not located by pe info but across edges
		if (outfp2 && footprint[t - 1] && !rd1gap)
		{
			if (ctgIdArray[t - 1] < 1)
				locate1read (t - 1);
			if(orienArray[t]=='+'){
					orientation = '-';
			}else{
					orientation = '+';
			}
			output1read (t - 1, outfp2,outfp5,orientation,1); //read1 in gap. by lbh
		}

		if (outfp2 && footprint[t] && !rd2gap)
		{
			if (ctgIdArray[t] < 1)
				locate1read (t);
			if(orienArray[t-1]=='+'){
					orientation = '-';
			}else{
					orientation = '+';
			}
			output1read (t, outfp2,outfp5,orientation,2); //read1 in gap. by lbh
		}
	}
}

//load contig index and length
void basicContigInfo (char *infile)
{
	char name[256], lldne[1024];
	FILE *fp;
	int length, bal_ed, num_all, num_long, index;

	sprintf (name, "%s.ContigIndex", infile);
	fp = ckopen (name, "r");
	fgets (lldne, sizeof (lldne), fp);
	sscanf (lldne + 8, "%d %d", &num_all, &num_long);
	printf ("%d edges in graph\n", num_all);
	num_ctg = num_all;
	contig_array = (CONTIG *) ckalloc ((num_all + 1) * sizeof (CONTIG));
	fgets (lldne, sizeof (lldne), fp);
	num_long = 0;

	while (fgets (lldne, sizeof (lldne), fp) != NULL)
	{
		sscanf (lldne, "%d %d %d", &index, &length, &bal_ed);
		contig_array[++num_long].length = length;
		contig_array[num_long].bal_edge = bal_ed + 1;

		if (index != num_long)
		{
			printf ("basicContigInfo: %d vs %d\n", index, num_long);
		}

		if (bal_ed == 0)
		{
			continue;
		}

		contig_array[++num_long].length = length;
		contig_array[num_long].bal_edge = -bal_ed + 1;
	}

	fclose (fp);
}

void prlRead2Ctg (char *libfile, char *outfile)
{
	long long i;
	char *src_name, *next_name, name[256];
	FILE *fo, *outfp2 = NULL,*outfp3=NULL,*outfp4=NULL;
	int maxReadNum, libNo, prevLibNo, insSize;
	boolean flag, pairs = 1;
	pthread_t threads[thrd_num];
	unsigned char thrdSignal[thrd_num + 1];
	PARAMETER paras[thrd_num];

	maxReadLen = 0;
	maxNameLen = 256;
	scan_libInfo (libfile);
	alloc_pe_mem (num_libs);

	if (!maxReadLen)
	{
		maxReadLen = 100;
	}

	printf ("In file: %s, max seq len %d, max name len %d\n\n", libfile, maxReadLen, maxNameLen);

	if (maxReadLen > maxReadLen4all)
	{
		maxReadLen4all = maxReadLen;
	}

	src_name = (char *) ckalloc ((maxNameLen + 1) * sizeof (char));
	next_name = (char *) ckalloc ((10*maxNameLen + 1) * sizeof (char));
	kmerBuffer = (Kmer *) ckalloc (buffer_size * sizeof (Kmer));
	hashBanBuffer = (ubyte8 *) ckalloc (buffer_size * sizeof (ubyte8));
	nodeBuffer = (kmer_t **) ckalloc (buffer_size * sizeof (kmer_t *));
	smallerBuffer = (boolean *) ckalloc (buffer_size * sizeof (boolean));
	maxReadNum = buffer_size / (maxReadLen - overlaplen + 1);
	maxReadNum = maxReadNum % 2 == 0 ? maxReadNum : maxReadNum - 1;	//make sure paired reads are processed at the same batch
	seqBuffer = (char **) ckalloc (maxReadNum * sizeof (char *));
	lenBuffer = (int *) ckalloc (maxReadNum * sizeof (int));
	indexArray = (unsigned int *) ckalloc ((maxReadNum + 1) * sizeof (unsigned int));
	ctgIdArray = (unsigned int *) ckalloc ((maxReadNum + 1) * sizeof (unsigned int));
	posArray = (int *) ckalloc ((maxReadNum + 1) * sizeof (int));
	orienArray = (char *) ckalloc ((maxReadNum + 1) * sizeof (char));
	footprint = (char *) ckalloc ((maxReadNum + 1) * sizeof (char));
	insSizeArray = (int*)ckalloc((maxReadNum+1)*sizeof(int));
	ctg2read=(READSET**)ckalloc((maxReadNum+1)*sizeof(READSET*));
	countContig=(unsigned int *)ckalloc((maxReadNum+1)*sizeof(unsigned int ));
	for (i = 0; i < maxReadNum; i++)
	{
		seqBuffer[i] = (char *) ckalloc (maxReadLen * sizeof (char));
	}

	rcSeq = (char **) ckalloc ((thrd_num + 1) * sizeof (char *));
	thrdSignal[0] = 0;

	if (1)
	{
		for (i = 0; i < thrd_num; i++)
		{
			rcSeq[i + 1] = (char *) ckalloc (maxReadLen * sizeof (char));
			thrdSignal[i + 1] = 0;
			paras[i].threadID = i;
			paras[i].mainSignal = &thrdSignal[0];
			paras[i].selfSignal = &thrdSignal[i + 1];
		}

		creatThrds (threads, paras);
	}

	if (!contig_array)
	{
		basicContigInfo (outfile);
	}

	sprintf (name, "%s.readInGap", outfile);
	outfp2 = ckopen (name, "wb");
	
	sprintf (name, "%s.readOnContig", outfile);
	fo = ckopen (name, "w");
	fprintf (fo, "read\tcontig\tpos\n");
	
	sprintf (name, "%s.ctg2Read", outfile);
	outfp3=ckopen(name,"w");
//	fprintf(outfp3,"readID\tcontigID\treadOffset\tcontigOffset\tAlignLength\tdirection\n");
	fprintf(outfp3,"read\tcontig\tpos\n");

	if(read_trace)
	{
		sprintf(name, "%s.readInformation", outfile);
		outfp4=ckopen(name,"w");
	}
	gzFile * outfp5, *outfp6;
	if (fill)
	{
		sprintf(name,"%s.shortreadInGap.gz",outfile);
   		outfp5 = gzopen(name,"w");
	
		sprintf(name,"%s.PEreadOnContig.gz",outfile);
		outfp6 = gzopen(name,"wb");	
	}
	
	readCounter = mapCounter = readsInGap = 0;
	kmer_c = n_solexa = read_c = i = libNo = readNumBack = gradsCounter = 0;
	prevLibNo = -1;

	int type = 0;		//deside the PE reads is good or bad

	while ((flag = read1seqInLib (seqBuffer[read_c], next_name, &(lenBuffer[read_c]), &libNo, pairs, 0, &type)) != 0)
	{
		if (type == -1)	//if the reads is bad, go back.
		{
			i--;
			if (lenBuffer[read_c - 1] >= overlaplen + 1)
			{
				kmer_c -= lenBuffer[read_c - 1] - overlaplen + 1;
			//	read_c--;
			}
			read_c--;
			n_solexa -= 2;
			continue;
		}
		if (libNo != prevLibNo)
		{
			prevLibNo = libNo;
			insSize = lib_array[libNo].avg_ins;
			ALIGNLEN = lib_array[libNo].map_len;
			printf ("current insert size %d, map_len %d\n", insSize, ALIGNLEN);
			if (insSize > 1000)
				ALIGNLEN = ALIGNLEN < 35 ? 35 : ALIGNLEN;
			else
				ALIGNLEN = ALIGNLEN < 32 ? 32 : ALIGNLEN;
		}
		insSizeArray[read_c] = insSize;
		
		if (insSize > 1000)
			ALIGNLEN = ALIGNLEN < (lenBuffer[read_c] / 2 + 1) ? (lenBuffer[read_c] / 2 + 1) : ALIGNLEN;

//		printf ("--- %lldth reads\n", i);
		if ((++i) % 100000000 == 0)
			printf ("--- %lldth reads\n", i);

		indexArray[read_c] = kmer_c;
		if (lenBuffer[read_c] >= overlaplen + 1)
			kmer_c += lenBuffer[read_c] - overlaplen + 1;
		read_c++;

		if (read_c == maxReadNum)
		{
			indexArray[read_c] = kmer_c;
			sendWorkSignal (2, thrdSignal);
			sendWorkSignal (1, thrdSignal);
			sendWorkSignal (3, thrdSignal);
			recordAlldgn (fo, insSize, outfp2,outfp3,outfp4, outfp5, outfp6);
			kmer_c = 0;
			read_c = 0;
		}
	}

	if (read_c)
	{
		indexArray[read_c] = kmer_c;
		sendWorkSignal (2, thrdSignal);
		sendWorkSignal (1, thrdSignal);
		sendWorkSignal (3, thrdSignal);
		recordAlldgn (fo, insSize, outfp2,outfp3,outfp4,outfp5, outfp6);
		printf ("Output %lld out of %lld (%.1f)%% reads in gaps\n", readsInGap, readCounter, (float) readsInGap / readCounter * 100);
	}

	printf ("%lld out of %lld (%.1f)%% reads mapped to contigs\n", mapCounter, readCounter, (float) mapCounter / readCounter * 100);
	sendWorkSignal (5, thrdSignal);
	thread_wait (threads);
	fclose (fo);
	sprintf (name, "%s.peGrads", outfile);
	fo = ckopen (name, "w");
	fprintf (fo, "grads&num: %d\t%lld\t%d\n", gradsCounter, n_solexa, maxReadLen4all);

	if (pairs)
	{
		if (gradsCounter)
			printf ("%d pe insert size, the largest boundary is %lld\n\n", gradsCounter, pes[gradsCounter - 1].PE_bound);
		else
		{
			printf ("no paired reads found\n");
		}

		for (i = 0; i < gradsCounter; i++)
		{
			fprintf (fo, "%d\t%lld\t%d\t%d\n", pes[i].insertS, pes[i].PE_bound, pes[i].rank, pes[i].pair_num_cut);
		}

		fclose (fo);
	}

	fclose (outfp2);
	fclose(outfp3);
	if(read_trace)
		fclose(outfp4);
	if(fill)
	{
		gzclose(outfp5);
		gzclose(outfp6);
	}
	free_pe_mem ();
	free_libs ();

	if (1)			// multi-threads
	{
		for (i = 0; i < thrd_num; i++)
		{
			free ((void *) rcSeq[i + 1]);
		}
	}

	free ((void *) rcSeq);

	for (i = 0; i < maxReadNum; i++)
	{
		free ((void *) seqBuffer[i]);
	}

	free ((void *) seqBuffer);
	free ((void *) lenBuffer);
	free ((void *) indexArray);
	free ((void *) kmerBuffer);
	free ((void *) smallerBuffer);
	free ((void *) hashBanBuffer);
	free ((void *) nodeBuffer);
	free ((void *) ctgIdArray);
	free ((void *) posArray);
	free ((void *) orienArray);
	free ((void *) footprint);
	free ((void *) src_name);
	free ((void *) next_name);

	if (contig_array)
	{
		free ((void *) contig_array);
		contig_array = NULL;
	}
}

static void thread_wait (pthread_t * threads)
{
	int i;

	for (i = 0; i < thrd_num; i++)
		if (threads[i] != 0)
		{
			pthread_join (threads[i], NULL);
		}
}

/********************* map long reads for gap filling ************************/
void prlLongRead2Ctg (char *libfile, char *outfile)
{
	long long i;
	char *src_name, *next_name, name[256];
	FILE *outfp2;
	int maxReadNum, libNo, prevLibNo;
	boolean flag, pairs = 0;
	pthread_t threads[thrd_num];
	unsigned char thrdSignal[thrd_num + 1];
	PARAMETER paras[thrd_num];

	maxReadLen = 0;
	maxNameLen = 256;
	scan_libInfo (libfile);

	if (!maxReadLen)
	{
		maxReadLen = 100;
	}

	int longReadLen = getMaxLongReadLen (num_libs);

	if (longReadLen < 1)	// no long reads
	{
		return;
	}

	maxReadLen4all = maxReadLen < longReadLen ? longReadLen : maxReadLen;
	printf ("In file: %s, long read len %d, max name len %d\n\n", libfile, longReadLen, maxNameLen);
	maxReadLen = longReadLen;
	src_name = (char *) ckalloc ((maxNameLen + 1) * sizeof (char));
	next_name = (char *) ckalloc ((maxNameLen + 1) * sizeof (char));
	kmerBuffer = (Kmer *) ckalloc (buffer_size * sizeof (Kmer));
	hashBanBuffer = (ubyte8 *) ckalloc (buffer_size * sizeof (ubyte8));
	nodeBuffer = (kmer_t **) ckalloc (buffer_size * sizeof (kmer_t *));
	smallerBuffer = (boolean *) ckalloc (buffer_size * sizeof (boolean));
	maxReadNum = buffer_size / (maxReadLen - overlaplen + 1);
	maxReadNum = maxReadNum % 2 == 0 ? maxReadNum : maxReadNum - 1;	//make sure paired reads are processed at the same batch
	seqBuffer = (char **) ckalloc (maxReadNum * sizeof (char *));
	lenBuffer = (int *) ckalloc (maxReadNum * sizeof (int));
	indexArray = (unsigned int *) ckalloc ((maxReadNum + 1) * sizeof (unsigned int));
	ctgIdArray = (unsigned int *) ckalloc ((maxReadNum + 1) * sizeof (unsigned int));
	posArray = (int *) ckalloc ((maxReadNum + 1) * sizeof (int));
	orienArray = (char *) ckalloc ((maxReadNum + 1) * sizeof (char));
	footprint = (char *) ckalloc ((maxReadNum + 1) * sizeof (char));

	for (i = 0; i < maxReadNum; i++)
	{
		seqBuffer[i] = (char *) ckalloc (maxReadLen * sizeof (char));
	}

	rcSeq = (char **) ckalloc ((thrd_num + 1) * sizeof (char *));
	deletion = (int *) ckalloc ((thrd_num + 1) * sizeof (int));
	thrdSignal[0] = 0;
	deletion[0] = 0;

	if (1)
	{
		for (i = 0; i < thrd_num; i++)
		{
			rcSeq[i + 1] = (char *) ckalloc (maxReadLen * sizeof (char));
			deletion[i + 1] = 0;
			thrdSignal[i + 1] = 0;
			paras[i].threadID = i;
			paras[i].mainSignal = &thrdSignal[0];
			paras[i].selfSignal = &thrdSignal[i + 1];
		}

		creatThrds (threads, paras);
	}

	if (!contig_array)
	{
		basicContigInfo (outfile);
	}

	sprintf (name, "%s.longReadInGap", outfile);
	outfp2 = ckopen (name, "wb");
	readCounter = 0;
	kmer_c = n_solexa = read_c = i = libNo = 0;
	prevLibNo = -1;

	int type = 0;		//deside the PE reads is good or bad

	while ((flag = read1seqInLib (seqBuffer[read_c], next_name, &(lenBuffer[read_c]), &libNo, pairs, 4, &type)) != 0)
	{
		if (type == -1)	//if the reads is bad, go back.
		{
			i--;
			if (lenBuffer[read_c - 1] >= overlaplen + 1)
			{
				kmer_c -= lenBuffer[read_c - 1] - overlaplen + 1;
			//	read_c--;
			}
			read_c--;
			n_solexa -= 2;
			continue;
		}
		if (libNo != prevLibNo)
		{
			prevLibNo = libNo;
			ALIGNLEN = lib_array[libNo].map_len;
			ALIGNLEN = ALIGNLEN < 35 ? 35 : ALIGNLEN;
			printf ("Map_len %d\n", ALIGNLEN);
		}

		if ((++i) % 100000000 == 0)
		{
			printf ("--- %lldth reads\n", i);
		}

		indexArray[read_c] = kmer_c;

		if (lenBuffer[read_c] >= overlaplen + 1)
		{
			kmer_c += lenBuffer[read_c] - overlaplen + 1;
		}

		read_c++;

		if (read_c == maxReadNum)
		{
			indexArray[read_c] = kmer_c;
			sendWorkSignal (2, thrdSignal);
			sendWorkSignal (1, thrdSignal);
			sendWorkSignal (3, thrdSignal);
			recordLongRead (outfp2);
			kmer_c = 0;
			read_c = 0;
		}
	}

	if (read_c)
	{
		indexArray[read_c] = kmer_c;
		sendWorkSignal (2, thrdSignal);
		sendWorkSignal (1, thrdSignal);
		sendWorkSignal (3, thrdSignal);
		recordLongRead (outfp2);
		printf ("Output %lld out of %lld (%.1f)%% reads in gaps\n", readsInGap, readCounter, (float) readsInGap / readCounter * 100);
	}

	sendWorkSignal (5, thrdSignal);
	thread_wait (threads);
	fclose (outfp2);
	free_libs ();

	if (1)			// multi-threads
	{
		for (i = 0; i < thrd_num; i++)
		{
			deletion[0] += deletion[i + 1];
			free ((void *) rcSeq[i + 1]);
		}
	}

	printf ("%d reads deleted\n", deletion[0]);
	free ((void *) rcSeq);
	free ((void *) deletion);

	for (i = 0; i < maxReadNum; i++)
	{
		free ((void *) seqBuffer[i]);
	}

	free ((void *) seqBuffer);
	free ((void *) lenBuffer);
	free ((void *) indexArray);
	free ((void *) kmerBuffer);
	free ((void *) smallerBuffer);
	free ((void *) hashBanBuffer);
	free ((void *) nodeBuffer);
	free ((void *) ctgIdArray);
	free ((void *) posArray);
	free ((void *) orienArray);
	free ((void *) footprint);
	free ((void *) src_name);
	free ((void *) next_name);
}
