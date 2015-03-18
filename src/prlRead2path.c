/*
 * prlRead2path.c
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
#include <stdinc.h>
#include "newhash.h"
#include <extfunc.h>
#include <extvab.h>

#define preARCBLOCKSIZE 100000
//static const Kmer kmerZero = { 0, 0, 0, 0 };

static unsigned int *arcCounters;
static int buffer_size = 100000000;
static long long markCounter = 0;
static unsigned int *fwriteBuf;
static unsigned char *markerOnEdge;

//buffer related varibles for chop kmer
static int read_c;
static char **rcSeq;
static char **seqBuffer;
static int *lenBuffer;

//edge and (K+1)mer related variables
static preARC **preArc_array;
static Kmer *mixBuffer;
static boolean *flagArray;	//indicate each item in mixBuffer where it's a (K+1)mer

// kmer related variables
static char **flags;
static int kmer_c;
static Kmer *kmerBuffer;
static ubyte8 *hashBanBuffer;
static kmer_t **nodeBuffer;
static boolean *smallerBuffer;
static int *indexArray;

static int *deletion;

static struct aiocb aio1;
static struct aiocb aio2;
static char *aioBuffer1;
static char *aioBuffer2;
static char *readBuffer1;
static char *readBuffer2;

static void parse1read (int t, int threadID);
static void search1kmerPlus (int j, unsigned char thrdID);
static void threadRoutine (void *thrdID);
static void searchKmer (int t, KmerSet * kset);
static void chopKmer4read (int t, int threadID);
static void thread_wait (pthread_t * threads);
static void thread_add1preArc (unsigned int from_ed, unsigned int to_ed, unsigned int thrdID);

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

	printf ("%d thread created prlRead2path\n", thrd_num);
}
// 2 -> 1 -> 3 -> 4 -> 6
static void threadRoutine (void *para)
{
	PARAMETER *prm;
	int i, t, j, start, finish;
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
				{
					continue;
				}

				searchKmer (i, KmerSets[id]);
			}

			*(prm->selfSignal) = 0;
		}
		else if (*(prm->selfSignal) == 2)
		{
			for (i = 0; i < read_c; i++)
			{
				if (i % thrd_num != id)
				{
					continue;
				}

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
				{
					continue;
				}

				parse1read (t, id + 1);
			}

			*(prm->selfSignal) = 0;
		}
		else if (*(prm->selfSignal) == 4)
		{
			//printf("thread %d, reads %d splay kmerplus\n",id,read_c);
			for (t = 0; t < read_c; t++)
			{
				start = indexArray[t];
				finish = indexArray[t + 1];

				for (j = start; j < finish; j++)
				{
					if (flagArray[j] == 0)
					{
#ifdef MER127
						if(mixBuffer[j].low2==0)
#endif
#ifdef MER63
						if(mixBuffer[j].low==0)
#endif
#ifdef MER31
						if(mixBuffer[j]==0)
#endif
						{
							break;
						}
					}
					else if (hashBanBuffer[j] % thrd_num == id)
					{
						//fprintf(stderr,"thread %d search for ban %lld\n",id,hashBanBuffer[j]);
						search1kmerPlus (j, id);
					}

					/*
					   if(flagArray[j]==0&&mixBuffer[j]==0)
					   break;
					   if(!flagArray[j]||(hashBanBuffer[j]%thrd_num)!=id)
					   continue;
					   search1kmerPlus(j,id);
					 */
				}
			}

			*(prm->selfSignal) = 0;
		}
		else if (*(prm->selfSignal) == 6)
		{
			for (t = 0; t < read_c; t++)
			{
				start = indexArray[t];
				finish = indexArray[t + 1];

				for (j = start; j < finish - 1; j++)
				{
#ifdef MER127
					if (mixBuffer[j].low2 == 0 || mixBuffer[j + 1].low2 == 0)
					{
						break;
					}

					if (mixBuffer[j].low2 % thrd_num != id)
					{
						continue;
					}

					thread_add1preArc (mixBuffer[j].low2, mixBuffer[j + 1].low2, id);
#endif
#ifdef MER63
					if (mixBuffer[j].low == 0 || mixBuffer[j + 1].low == 0)
					{
						break;
					}

					if (mixBuffer[j].low % thrd_num != id)
					{
						continue;
					}

					thread_add1preArc (mixBuffer[j].low, mixBuffer[j + 1].low, id);
#endif
#ifdef MER31
					if (mixBuffer[j] == 0 || mixBuffer[j + 1] == 0)
					{
						break;
					}

					if (mixBuffer[j] % thrd_num != id)
					{
						continue;
					}

					thread_add1preArc (mixBuffer[j], mixBuffer[j + 1], id);
#endif				
				}
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
	char *src_seq = seqBuffer[t];
	char *bal_seq = rcSeq[threadID];
	int len_seq = lenBuffer[t];
	int j, bal_j;
	ubyte8 hash_ban, bal_hash_ban;
	Kmer word, bal_word;
	int index;

	//mao 2011 10 8
	Kmer InvalidKmer ;
	InvalidKmer=kmerZero;
	int n_num;
	word = kmerZero;


	for (index = 0; index < overlaplen; index++)
	{
#ifdef MER127
		word = KmerLeftBitMoveBy2(word);
		word.low2 |= src_seq[index];
#endif
#ifdef MER63
		word = KmerLeftBitMoveBy2(word);
		word.low |= src_seq[index];
#endif
#ifdef MER31
		word = KmerLeftBitMoveBy2(word);
		word += src_seq[index];
#endif

		//mao 2011 10 8
		if(src_seq[index] == 4)
		        n_num = overlaplen;
		else if(n_num >0)
		        n_num--;
	}

	reverseComplementSeq (src_seq, len_seq, bal_seq);
	// complementary node
	bal_word = reverseComplement (word, overlaplen);
	bal_j = len_seq - 0 - overlaplen;	//  0;
	index = indexArray[t];

	//mao 2011 10 8
	if(n_num > 0 && N_kmer)
	{
		hash_ban = hash_kmer (InvalidKmer);
		hashBanBuffer[index] = hash_ban;
		kmerBuffer[index] = InvalidKmer;
		smallerBuffer[index] = 1;
	}
	else if (KmerSmaller (word, bal_word))
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
		bal_word = reverseComplement (word, overlaplen);

		//mao 2011 10 8
		if(src_seq[j - 1 + overlaplen] == 4)
		        n_num = overlaplen;
		else if(n_num >0)
		        n_num--;

		//mao 2011 10 8
		if(n_num > 0 && N_kmer)
		{
			hash_ban = hash_kmer (InvalidKmer);
			hashBanBuffer[index] = hash_ban;
			kmerBuffer[index] = InvalidKmer;
			smallerBuffer[index] = 1;
		}
		else if (KmerSmaller (word, bal_word))
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
#ifdef MER127
	if( kmerBuffer[t].low2== 0 && N_kmer)
#endif
#ifdef MER63
	if( kmerBuffer[t].low== 0 && N_kmer)
#endif
#ifdef MER31
	if( kmerBuffer[t]== 0 && N_kmer)
#endif
		return ;
	kmer_t *node;
	boolean found = search_kmerset (kset, kmerBuffer[t], &node);

	if (!found)
	{
#ifdef MER127
		printf ("searchKmer: kmer %llx %llx %llx %llx is not found\n",
			kmerBuffer[t].high1, kmerBuffer[t].low1, kmerBuffer[t].high2, kmerBuffer[t].low2);
#endif
#ifdef MER63
		printf ("searchKmer: kmer %llx %llx is not found\n",
			kmerBuffer[t].high, kmerBuffer[t].low);
#endif
#ifdef MER31
		printf ("searchKmer: kmer %llx is not found\n", kmerBuffer[t]);
#endif
	}

	nodeBuffer[t] = node;
}

static preARC *getPreArcBetween (unsigned int from_ed, unsigned int to_ed)
{
	preARC *parc;

	parc = preArc_array[from_ed];

	while (parc)
	{
		if (parc->to_ed == to_ed)
		{
			return parc;
		}

		parc = parc->next;
	}

	return parc;
}

static void thread_add1preArc (unsigned int from_ed, unsigned int to_ed, unsigned int thrdID)
{
	preARC *parc = getPreArcBetween (from_ed, to_ed);

	if (parc)
	{
		parc->multiplicity++;
	}
	else
	{
		parc = prlAllocatePreArc (to_ed, preArc_mem_managers[thrdID]);
		arcCounters[thrdID]++;
		parc->next = preArc_array[from_ed];
		preArc_array[from_ed] = parc;
	}
}

static void memoAlloc4preArc ()
{
	unsigned int i;

	preArc_array = (preARC **) ckalloc ((num_ed + 1) * sizeof (preARC *));

	for (i = 0; i <= num_ed; i++)
	{
		preArc_array[i] = NULL;
	}
}

static void memoFree4preArc ()
{
	prlDestroyPreArcMem ();

	if (preArc_array)
	{
		free ((void *) preArc_array);
	}
}

static void output_arcs (char *outfile)
{
	unsigned int i;
	char name[256];
	FILE *outfp, *outfp2 = NULL;
	preARC *parc;

	sprintf (name, "%s.preArc", outfile);
	outfp = ckopen (name, "w");

	if (repsTie)
	{
		sprintf (name, "%s.markOnEdge", outfile);
		outfp2 = ckopen (name, "w");
	}

	markCounter = 0;

	for (i = 1; i <= num_ed; i++)
	{
		if (repsTie)
		{
			markCounter += markerOnEdge[i];
			fprintf (outfp2, "%d\n", markerOnEdge[i]);
		}

		parc = preArc_array[i];

		if (!parc)
		{
			continue;
		}

		fprintf (outfp, "%u", i);

		while (parc)
		{
			fprintf (outfp, " %u %u", parc->to_ed, parc->multiplicity);
			parc = parc->next;
		}

		fprintf (outfp, "\n");
	}

	fclose (outfp);

	if (repsTie)
	{
		fclose (outfp2);
		printf ("%lld markers counted\n", markCounter);
	}
}

static void recordPathBin (FILE * outfp)
{
	int t, j, start, finish;
	unsigned char counter;

	for (t = 0; t < read_c; t++)
	{
		start = indexArray[t];
		finish = indexArray[t + 1];

#ifdef MER127
		if (finish - start < 3 || mixBuffer[start].low2 == 0 || mixBuffer[start + 1].low2 == 0 || mixBuffer[start + 2].low2 == 0)
		{
			continue;
		}
#endif
#ifdef MER63
		if (finish - start < 3 || mixBuffer[start].low == 0 || mixBuffer[start + 1].low == 0 || mixBuffer[start + 2].low== 0)
		{
			continue;
		}
#endif
#ifdef MER31
		if (finish - start < 3 || mixBuffer[start] == 0 || mixBuffer[start + 1] == 0 || mixBuffer[start + 2] == 0)
		{
			continue;
		}
#endif
		counter = 0;

#ifdef MER127
		for (j = start; j < finish; j++)
		{
			if (mixBuffer[j].low2 == 0)
				break;
			fwriteBuf[counter++] = (unsigned int) mixBuffer[j].low2;
			if (markerOnEdge[mixBuffer[j].low2] <  255)
				markerOnEdge[mixBuffer[j].low2]++;
			markCounter++;
		}
#endif
#ifdef MER63
		for (j = start; j < finish; j++)
		{
			if (mixBuffer[j].low == 0 )
				break;
			fwriteBuf[counter++] = (unsigned int) mixBuffer[j].low;
			if (markerOnEdge[mixBuffer[j].low] <  255)
				markerOnEdge[mixBuffer[j].low]++;
			markCounter++;
		}
#endif
#ifdef MER31
		for (j = start; j < finish; j++)
		{
			if (mixBuffer[j] == 0)
				break;
			fwriteBuf[counter++] = (unsigned int) mixBuffer[j];
			if (markerOnEdge[mixBuffer[j]] <  255)
				markerOnEdge[mixBuffer[j]]++;
			markCounter++;
		}
#endif
		fwrite (&counter, sizeof (char), 1, outfp);
		fwrite (fwriteBuf, sizeof (unsigned int), (int) counter, outfp);
	}
}

static void search1kmerPlus (int j, unsigned char thrdID)
{
	kmer_t *node;
	boolean found = search_kmerset (KmerSetsPatch[thrdID], mixBuffer[j], &node);

	if (!found)
	{
		mixBuffer[j] =kmerZero;
		return;
	}
#ifdef MER127
	if (smallerBuffer[j])
	{
		mixBuffer[j].low2 = node->l_links;
	}
	else
	{
		mixBuffer[j].low2 = node->l_links + node->twin - 1;
	}
#endif
#ifdef MER63
	if (smallerBuffer[j])
	{
		mixBuffer[j].low= node->l_links;
	}
	else
	{
		mixBuffer[j].low= node->l_links + node->twin - 1;
	}
#endif
#ifdef MER31
	if (smallerBuffer[j])
	{
		mixBuffer[j] = node->l_links;
	}
	else
	{
		mixBuffer[j] = node->l_links + node->twin - 1;
	}
#endif
}

static void parse1read (int t, int threadID)
{
	unsigned int j, retain = 0;
	unsigned int edge_index = 0;
	kmer_t *node;
	boolean isSmaller;
	Kmer wordplus, bal_wordplus;
	unsigned int start, finish, pos;
	Kmer prevKmer, currentKmer;
	boolean IsPrevKmer = 0;

	start = indexArray[t];
	finish = indexArray[t + 1];
	pos = start;

	for (j = start; j < finish; j++)
	{
#ifdef MER127
		if( kmerBuffer[j].low2== 0 && N_kmer)
#endif
#ifdef MER63
		if( kmerBuffer[j].low== 0 && N_kmer)
#endif
#ifdef MER31
		if( kmerBuffer[j]== 0 && N_kmer)
#endif
		{
			IsPrevKmer=0;
			continue;
		}
		node = nodeBuffer[j];

		//extract edges or keep kmers
		if ((node->deleted) || (node->linear && !node->inEdge))	// deleted or in a floating loop
		{
			if (retain < 2)
			{
				retain = 0;
				pos = start;
			}
			else
			{
				break;
			}

			continue;
		}

		isSmaller = smallerBuffer[j];

		if (node->linear)
		{
			if (isSmaller)
			{
				edge_index = node->l_links;
			}
			else
			{
				edge_index = node->l_links + node->twin - 1;
			}
#ifdef MER127
			if (retain == 0 || IsPrevKmer)
			{
				retain++;
				mixBuffer[pos].low2 = edge_index;
				flagArray[pos++] = 0;
				IsPrevKmer = 0;
			}
			else if (edge_index != mixBuffer[pos - 1].low2)
			{
				retain++;
				mixBuffer[pos].low2 = edge_index;
				flagArray[pos++] = 0;
			}
#endif
#ifdef MER63
			if (retain == 0 || IsPrevKmer)
			{
				retain++;
				mixBuffer[pos].low= edge_index;
				flagArray[pos++] = 0;
				IsPrevKmer = 0;
			}
			else if (edge_index != mixBuffer[pos - 1].low)
			{
				retain++;
				mixBuffer[pos].low= edge_index;
				flagArray[pos++] = 0;
			}
#endif
#ifdef MER31
			if (retain == 0 || IsPrevKmer)
			{
				retain++;
				mixBuffer[pos] = edge_index;
				flagArray[pos++] = 0;
				IsPrevKmer = 0;
			}
			else if (edge_index != mixBuffer[pos - 1])
			{
				retain++;
				mixBuffer[pos] = edge_index;
				flagArray[pos++] = 0;
			}
#endif
		}
		else
		{
			if (isSmaller)
			{
				currentKmer = node->seq;
			}
			else
			{
				currentKmer = reverseComplement (node->seq, overlaplen);
			}

			if (IsPrevKmer)
			{
				retain++;
				wordplus = KmerPlus (prevKmer, lastCharInKmer (currentKmer));
				bal_wordplus = reverseComplement (wordplus, overlaplen + 1);

				if (KmerSmaller (wordplus, bal_wordplus))
				{
					smallerBuffer[pos] = 1;
					hashBanBuffer[pos] = hash_kmer (wordplus);
					mixBuffer[pos] = wordplus;
				}
				else
				{
					smallerBuffer[pos] = 0;
					hashBanBuffer[pos] = hash_kmer (bal_wordplus);
					mixBuffer[pos] = bal_wordplus;
				}

				//  fprintf(stderr,"%lld\n",hashBanBuffer[pos]);
				flagArray[pos++] = 1;
			}

			IsPrevKmer = 1;
			prevKmer = currentKmer;
		}
	}

	/*
	   for(j=start;j<pos;j++)
	   fprintf(stderr,"%d ",flagArray[j]);
	   fprintf(stderr,"\n");
	 */
	if (retain < 1)
	{
		deletion[threadID]++;
	}

	if (retain < 2)
	{
		flagArray[start] = 0;
		mixBuffer[start] = kmerZero;
		return;
	}

	if ((pos - start) != retain)
	{
		printf ("read %d, %d vs %d\n", t, retain, edge_index - start);
	}

	if (pos < finish)
	{
		flagArray[pos] = 0;
		mixBuffer[pos] = kmerZero;
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

void prlRead2edge (char *libfile, char *outfile)
{
	char *cach1;
	char *cach2;
	unsigned char asm_ctg = 1;

	long long i;
	char name[256], *src_name, *next_name;
	FILE *outfp = NULL;
	int maxReadNum, libNo;
	boolean flag, pairs = 0;
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

	maxReadLen4all = maxReadLen;
	printf ("In file: %s, max seq len %d, max name len %d\n\n", libfile, maxReadLen, maxNameLen);

	if (repsTie)
	{
		sprintf (name, "%s.path", outfile);
		outfp = ckopen (name, "wb");
	}

	src_name = (char *) ckalloc ((maxNameLen + 1) * sizeof (char));
	next_name = (char *) ckalloc ((10*maxNameLen + 1) * sizeof (char));
	kmerBuffer = (Kmer *) ckalloc (buffer_size * sizeof (Kmer));
	mixBuffer = (Kmer *) ckalloc (buffer_size * sizeof (Kmer));
	hashBanBuffer = (ubyte8 *) ckalloc (buffer_size * sizeof (ubyte8));
	nodeBuffer = (kmer_t **) ckalloc (buffer_size * sizeof (kmer_t *));
	smallerBuffer = (boolean *) ckalloc (buffer_size * sizeof (boolean));
	flagArray = (boolean *) ckalloc (buffer_size * sizeof (boolean));
	maxReadNum = buffer_size / (maxReadLen - overlaplen + 1);
	//printf("buffer for at most %d reads\n",maxReadNum);
	
	int maxAIOSize = 32768;/*
	aioBuffer1 = (char *) ckalloc ((maxAIOSize) * sizeof (char));
	aioBuffer2 = (char *) ckalloc ((maxAIOSize) * sizeof (char));
	readBuffer1 = (char *) ckalloc ((maxAIOSize + 1024) * sizeof (char));	//(char *)ckalloc(maxAIOSize*sizeof(char));
	readBuffer2 = (char *) ckalloc ((maxAIOSize + 1024) * sizeof (char));
	cach1 = (char *) ckalloc (1024 * sizeof (char));
	cach2 = (char *) ckalloc (1024 * sizeof (char));
	memset(cach1,'\0',1024);
	memset(cach2,'\0',1024);*/
        aioBuffer1 = (char *) ckalloc ((maxAIOSize) * sizeof (char));
        aioBuffer2 = (char *) ckalloc ((maxAIOSize) * sizeof (char));
        readBuffer1 = (char *) ckalloc ((maxAIOSize + (maxReadLen+1024)) * sizeof (char));      //(char *)ckalloc(maxAIOSize*sizeof(char));     //1024
        readBuffer2 = (char *) ckalloc ((maxAIOSize + (maxReadLen+1024)) * sizeof (char));      //1024
        cach1 = (char *) ckalloc ((maxReadLen+1024) * sizeof (char));   //1024
        cach2 = (char *) ckalloc ((maxReadLen+1024) * sizeof (char));   //1024
        memset(cach1,'\0',(maxReadLen+1024));   //1024
        memset(cach2,'\0',(maxReadLen+1024));   //1024


	seqBuffer = (char **) ckalloc (maxReadNum * sizeof (char *));
	lenBuffer = (int *) ckalloc (maxReadNum * sizeof (int));
	indexArray = (int *) ckalloc ((maxReadNum + 1) * sizeof (int));

	for (i = 0; i < maxReadNum; i++)
	{
		seqBuffer[i] = (char *) ckalloc (maxReadLen * sizeof (char));
	}

	memoAlloc4preArc ();
	flags = (char **) ckalloc ((thrd_num + 1) * sizeof (char *));
	deletion = (int *) ckalloc ((thrd_num + 1) * sizeof (int));
	rcSeq = (char **) ckalloc ((thrd_num + 1) * sizeof (char *));

	if (repsTie)
	{
		markerOnEdge = (unsigned char *) ckalloc ((num_ed + 1) * sizeof (unsigned char));

		for (i = 1; i <= num_ed; i++)
		{
			markerOnEdge[i] = 0;
		}

		fwriteBuf = (unsigned int *) ckalloc ((maxReadLen - overlaplen + 1) * sizeof (unsigned int));
	}

	thrdSignal[0] = 0;

	if (1)
	{
		preArc_mem_managers = (MEM_MANAGER **) ckalloc (thrd_num * sizeof (MEM_MANAGER *));
		arcCounters = (unsigned int *) ckalloc (thrd_num * sizeof (unsigned int));

		for (i = 0; i < thrd_num; i++)
		{
			arcCounters[i] = 0;
			preArc_mem_managers[i] = createMem_manager (preARCBLOCKSIZE, sizeof (preARC));
			deletion[i + 1] = 0;
			flags[i + 1] = (char *) ckalloc (2 * maxReadLen * sizeof (char));
			rcSeq[i + 1] = (char *) ckalloc (maxReadLen * sizeof (char));
			thrdSignal[i + 1] = 0;
			paras[i].threadID = i;
			paras[i].mainSignal = &thrdSignal[0];
			paras[i].selfSignal = &thrdSignal[i + 1];
		}

		creatThrds (threads, paras);
	}

	if (1)
	{
		deletion[0] = 0;
		flags[0] = (char *) ckalloc (2 * maxReadLen * sizeof (char));
		rcSeq[0] = (char *) ckalloc (maxReadLen * sizeof (char));
	}

	kmer_c = n_solexa = read_c = i = libNo = readNumBack = gradsCounter = 0;
	int t0, t1, t2, t3, t4, t5, t6;

	t0 = t1 = t2 = t3 = t4 = t5 = t6 = 0;
	time_t read_start, read_end, time_bef, time_aft;

	time (&read_start);
	
		while (openNextFile (&libNo, pairs, asm_ctg))
	{
		if (lib_array[libNo].curr_type == 4)
		{
			int type = 0;	//deside the PE reads is good or bad

			while ((flag = read1seqInLibBam (seqBuffer[read_c], next_name, &(lenBuffer[read_c]), &libNo, pairs, 1, &type)) != 0)
			{
				if (type == -1)	//if the reads is bad, go back.
				{
					i--;
					if (lenBuffer[read_c - 1] >= overlaplen + 1)
					{
						kmer_c -= lenBuffer[read_c - 1] - overlaplen + 1;
						read_c--;
					}
					n_solexa -= 2;
					continue;
				}
				if ((++i) % 1000000 == 0)
				{
					printf ("--- %lldth reads\n", i);
				}

				if (lenBuffer[read_c] < overlaplen + 1)
				{
					continue;
				}

				//if(lenBuffer[read_c]>70)
				//    lenBuffer[read_c] = 70;
				//else if(lenBuffer[read_c]>40)
				//    lenBuffer[read_c] = 40;

				indexArray[read_c] = kmer_c;
				kmer_c += lenBuffer[read_c] - overlaplen + 1;
				read_c++;

				if (read_c == maxReadNum)
				{
					indexArray[read_c] = kmer_c;
					time (&read_end);
					t0 += read_end - read_start;
					time (&time_bef);
					sendWorkSignal (2, thrdSignal);
					time (&time_aft);
					t1 += time_aft - time_bef;
					time (&time_bef);
					sendWorkSignal (1, thrdSignal);
					time (&time_aft);
					t2 += time_aft - time_bef;
					time (&time_bef);
					sendWorkSignal (3, thrdSignal);
					time (&time_aft);
					t3 += time_aft - time_bef;
					time (&time_bef);
					sendWorkSignal (4, thrdSignal);
					time (&time_aft);
					t4 += time_aft - time_bef;
					time (&time_bef);
					sendWorkSignal (6, thrdSignal);
					time (&time_aft);
					t5 += time_aft - time_bef;
					time (&time_bef);

					//recordPreArc();
					if (repsTie)
					{
						recordPathBin (outfp);
					}

					time (&time_aft);
					t6 += time_aft - time_bef;
					//output_path(read_c,edge_no,flags,outfp);
					kmer_c = 0;
					read_c = 0;
					time (&read_start);
				}
			}
		}
		else if (lib_array[libNo].curr_type == 1 || lib_array[libNo].curr_type == 2)
		{
			initAIO (&aio1, aioBuffer1, fileno (lib_array[libNo].fp1), maxAIOSize);
			initAIO (&aio2, aioBuffer2, fileno (lib_array[libNo].fp2), maxAIOSize);
			int offset1, offset2, flag1, flag2, rt1, rt2;

			offset1 = offset2 = 0;
			rt1 = aio_read (&aio1);
			rt2 = aio_read (&aio2);
			flag1 = AIORead (&aio1, &offset1, readBuffer1, cach1, &rt1, lib_array[libNo].curr_type);
			flag2 = AIORead (&aio2, &offset2, readBuffer2, cach2, &rt2, lib_array[libNo].curr_type);
			if(flag1 && flag2)
			{
				int start1, start2, turn;

				start1 = start2 = 0;
				turn = 1;
				while (start1 < offset1 || start2 < offset2)
				{
					if (turn == 1)
					{
						turn = 2;
						readseqInLib (seqBuffer[read_c], next_name, &(lenBuffer[read_c]), readBuffer1, &start1, offset1, libNo);
						if ((++i) % 1000000 == 0)
							printf ("--- %lldth reads\n", i);
/*						if (lenBuffer[read_c] < overlaplen + 1)
							continue;*/
						if (lenBuffer[read_c] < overlaplen + 1)
						{
							if(start1>=offset1)
							{
								start1=0;
								flag1=AIORead (&aio1, &offset1, readBuffer1, cach1, &rt1, lib_array[libNo].curr_type);
							}
							continue;
						}

						indexArray[read_c] = kmer_c;
						kmer_c += lenBuffer[read_c] - overlaplen + 1;
						read_c++;
						if(start1>=offset1){
							start1=0;
							flag1=AIORead (&aio1, &offset1, readBuffer1, cach1, &rt1, lib_array[libNo].curr_type);
						}
						if (read_c == maxReadNum) {
							indexArray[read_c] = kmer_c;

							time (&read_end);
							t0 += read_end - read_start;
							time (&time_bef);
							sendWorkSignal (2, thrdSignal);
							time (&time_aft);
							t1 += time_aft - time_bef;
							time (&time_bef);
							sendWorkSignal (1, thrdSignal);
							time (&time_aft);
							t2 += time_aft - time_bef;
							time (&time_bef);
							sendWorkSignal (3, thrdSignal);
							time (&time_aft);
							t3 += time_aft - time_bef;
							time (&time_bef);
							sendWorkSignal (4, thrdSignal);
							time (&time_aft);
							t4 += time_aft - time_bef;
							time (&time_bef);
							sendWorkSignal (6, thrdSignal);
							time (&time_aft);
							t5 += time_aft - time_bef;
							time (&time_bef);

							//recordPreArc();
							if (repsTie)
								recordPathBin (outfp);
							time (&time_aft);
							t6 += time_aft - time_bef;
							//output_path(read_c,edge_no,flags,outfp);
							kmer_c = 0;
							read_c = 0;
							time (&read_start);
						}
						continue;
					}
					if (turn == 2)
					{
						turn = 1;
						readseqInLib (seqBuffer[read_c], next_name, &(lenBuffer[read_c]), readBuffer2, &start2, offset2, libNo);
						if ((++i) % 1000000 == 0)
							printf ("--- %lldth reads\n", i);
/*						if (lenBuffer[read_c] < overlaplen + 1)
							continue;*/
						if (lenBuffer[read_c] < overlaplen + 1)
						{
							if((flag2 == 2) && (start2 >= offset2))
								break;

							if(start2 >= offset2)
							{
								start2=0;
								flag2 = AIORead (&aio2, &offset2, readBuffer2, cach2, &rt2, lib_array[libNo].curr_type);
							}
							continue;
						}

						indexArray[read_c] = kmer_c;
						kmer_c += lenBuffer[read_c] - overlaplen + 1;
						read_c++;
						if((flag2 == 2) && (start2 >= offset2))
							break;
						if(start2 >= offset2){
			                        	start2=0;
							flag2 = AIORead (&aio2, &offset2, readBuffer2, cach2, &rt2, lib_array[libNo].curr_type);
						}
						if (read_c == maxReadNum){
							indexArray[read_c] = kmer_c;

							time (&read_end);
							t0 += read_end - read_start;
							time (&time_bef);
							sendWorkSignal (2, thrdSignal);
							time (&time_aft);
							t1 += time_aft - time_bef;
							time (&time_bef);
							sendWorkSignal (1, thrdSignal);
							time (&time_aft);
							t2 += time_aft - time_bef;
							time (&time_bef);
							sendWorkSignal (3, thrdSignal);
							time (&time_aft);
							t3 += time_aft - time_bef;
							time (&time_bef);
							sendWorkSignal (4, thrdSignal);
							time (&time_aft);
							t4 += time_aft - time_bef;
							time (&time_bef);
							sendWorkSignal (6, thrdSignal);
							time (&time_aft);
							t5 += time_aft - time_bef;
							time (&time_bef);

							//recordPreArc();
							if (repsTie)
								recordPathBin (outfp);
							time (&time_aft);
							t6 += time_aft - time_bef;
							//output_path(read_c,edge_no,flags,outfp);
							kmer_c = 0;
							read_c = 0;
							time (&read_start);
						}
						continue;
					}
				}
			}
		}
		else
		{
			initAIO (&aio1, aioBuffer1, fileno (lib_array[libNo].fp1), maxAIOSize);
			int offset, flag1, rt;

			offset = 0;
			rt = aio_read (&aio1);
			while ((flag1 = AIORead (&aio1, &offset, readBuffer1, cach1, &rt, lib_array[libNo].curr_type)))
			{
				int start = 0;

				while (start < offset)
				{
					readseqInLib (seqBuffer[read_c], next_name, &(lenBuffer[read_c]), readBuffer1, &start, offset, libNo);
					if ((++i) % 1000000 == 0)
						printf ("--- %lld reads\n", i);
					if (lenBuffer[read_c] < overlaplen + 1)
						continue;
					indexArray[read_c] = kmer_c;
					kmer_c += lenBuffer[read_c] - overlaplen + 1;
					read_c++;

					if (read_c > maxReadNum - 1024)
					{
						indexArray[read_c] = kmer_c;

						time (&read_end);
						t0 += read_end - read_start;
						time (&time_bef);
						sendWorkSignal (2, thrdSignal);
						time (&time_aft);
						t1 += time_aft - time_bef;
						time (&time_bef);
						sendWorkSignal (1, thrdSignal);
						time (&time_aft);
						t2 += time_aft - time_bef;
						time (&time_bef);
						sendWorkSignal (3, thrdSignal);
						time (&time_aft);
						t3 += time_aft - time_bef;
						time (&time_bef);
						sendWorkSignal (4, thrdSignal);
						time (&time_aft);
						t4 += time_aft - time_bef;
						time (&time_bef);
						sendWorkSignal (6, thrdSignal);
						time (&time_aft);
						t5 += time_aft - time_bef;
						time (&time_bef);

						//recordPreArc();
						if (repsTie)
							recordPathBin (outfp);
						time (&time_aft);
						t6 += time_aft - time_bef;
						//output_path(read_c,edge_no,flags,outfp);
						kmer_c = 0;
						read_c = 0;
						time (&read_start);
					}
				}
				if (flag1 == 2)
					break;

			}
		}
	}

	printf ("%lld reads processed\n", i);
	printf ("time %d,%d,%d,%d,%d,%d,%d\n", t0, t1, t2, t3, t4, t5, t6);

	if (read_c)
	{
		indexArray[read_c] = kmer_c;
		sendWorkSignal (2, thrdSignal);
		sendWorkSignal (1, thrdSignal);
		sendWorkSignal (3, thrdSignal);
		sendWorkSignal (4, thrdSignal);
		sendWorkSignal (6, thrdSignal);

		//recordPreArc();
		if (repsTie)
		{
			recordPathBin (outfp);
		}
	}

	printf ("%lld markers outputed\n", markCounter);
	sendWorkSignal (5, thrdSignal);
	thread_wait (threads);
	output_arcs (outfile);
	memoFree4preArc ();

	if (1)			// multi-threads
	{
		arcCounter = 0;

		for (i = 0; i < thrd_num; i++)
		{
			arcCounter += arcCounters[i];
			free ((void *) flags[i + 1]);
			deletion[0] += deletion[i + 1];
			free ((void *) rcSeq[i + 1]);
		}
	}

	if (1)
	{
		free ((void *) flags[0]);
		free ((void *) rcSeq[0]);
	}

	printf ("done mapping reads, %d reads deleted, %lld arcs created\n", deletion[0], arcCounter);

	if (repsTie)
	{
		free ((void *) markerOnEdge);
		free ((void *) fwriteBuf);
	}

	free ((void *) arcCounters);
	free ((void *) rcSeq);

	for (i = 0; i < maxReadNum; i++)
	{
		free ((void *) seqBuffer[i]);
	}

	free ((void *) seqBuffer);
	free ((void *) lenBuffer);
	free ((void *) indexArray);
	free ((void *) flags);
	free ((void *) deletion);
	free ((void *) kmerBuffer);
	free ((void *) mixBuffer);
	free ((void *) smallerBuffer);
	free ((void *) flagArray);
	free ((void *) hashBanBuffer);
	free ((void *) nodeBuffer);
	free ((void *) src_name);
	free ((void *) next_name);
	free ((void *) aioBuffer1);
    free ((void *) aioBuffer2);
    free ((void *) readBuffer1);
    free ((void *) readBuffer2);
    free ((void *) cach1);
    free ((void *) cach2);

	if (repsTie)
	{
		fclose (outfp);
	}

	free_pe_mem ();
	free_libs ();
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
