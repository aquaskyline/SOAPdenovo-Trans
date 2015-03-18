/*
 * attachPEinfo.c
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
#include "stack.h"
#include "zlib.h"

#define CNBLOCKSIZE 10000
static STACK *isStack;
static int ignorePE1, ignorePE2, ignorePE3, static_flag;
static int onsameCtgPE;
static unsigned long long peSUM;

static boolean staticF;

static int existCounter;

static int calcuIS (STACK * intStack);

static int cmp_pe (const void *a, const void *b)
{
	PE_INFO *A, *B;

	A = (PE_INFO *) a;
	B = (PE_INFO *) b;

	if (A->rank > B->rank)
	{
		return 1;
	}
	else if (A->rank == B->rank)
	{
		return 0;
	}
	else
	{
		return -1;
	}
}

void loadPEgrads (char *infile)
{
	FILE *fp;
	char name[256], line[1024];
	int i;
	boolean rankSet = 1;

	sprintf (name, "%s.peGrads", infile);
	fp = fopen (name, "r");

	if (!fp)
	{
		printf ("can not open file %s\n", name);
		gradsCounter = 0;
		return;
	}

	while (fgets (line, sizeof (line), fp) != NULL)
	{
		if (line[0] == 'g')
		{
			sscanf (line + 10, "%d %lld %d", &gradsCounter, &n_solexa, &maxReadLen);
			printf ("there're %d grads, %lld reads, max read len %d\n", gradsCounter, n_solexa, maxReadLen);
			break;
		}
	}

	alloc_pe_mem (gradsCounter);

	for (i = 0; i < gradsCounter; i++)
	{
		fgets (line, sizeof (line), fp);
		pes[i].rank = 0;
		sscanf (line, "%d %lld %d %d", &(pes[i].insertS), &(pes[i].PE_bound), &(pes[i].rank), &(pes[i].pair_num_cut));

		if (pes[i].rank < 1)
		{
			rankSet = 0;
		}
	}

	fclose (fp);

	if (rankSet)
	{
		qsort (&pes[0], gradsCounter, sizeof (PE_INFO), cmp_pe);
		return;
	}

	int lastRank = 0;

	for (i = 0; i < gradsCounter; i++)
	{
		if (i == 0)
		{
			pes[i].rank = ++lastRank;
		}
		else if (pes[i].insertS < 300)
		{
			pes[i].rank = lastRank;
		}
		else if (pes[i].insertS < 800)
		{
			if (pes[i - 1].insertS < 300)
			{
				pes[i].rank = ++lastRank;
			}
			else
			{
				pes[i].rank = lastRank;
			}
		}
		else if (pes[i].insertS < 3000)
		{
			if (pes[i - 1].insertS < 800)
			{
				pes[i].rank = ++lastRank;
			}
			else
			{
				pes[i].rank = lastRank;
			}
		}
		else if (pes[i].insertS < 7000)
		{
			if (pes[i - 1].insertS < 3000)
			{
				pes[i].rank = ++lastRank;
			}
			else
			{
				pes[i].rank = lastRank;
			}
		}
		else
		{
			if (pes[i - 1].insertS < 7000)
			{
				pes[i].rank = ++lastRank;
			}
			else
			{
				pes[i].rank = lastRank;
			}
		}
	}
}

CONNECT *add1Connect (unsigned int e1, unsigned int e2, int gap, int weight, boolean inherit)
{
	if (e1 == e2 || e1 == getTwinCtg (e2))
	{
		return NULL;
	}

	CONNECT *connect = NULL;
	long long sum;

	if (weight > 255)
	{
		weight = 255;
	}

	connect = getCntBetween (e1, e2);

	if (connect)
	{
		if (!weight)
		{
			return connect;
		}

		existCounter++;

		if (!inherit)
		{
			sum = connect->weightNotInherit * connect->gapLen + gap * weight;
			connect->gapLen = sum / (connect->weightNotInherit + weight);

			if (connect->weightNotInherit + weight <= 255)
			{
				connect->weightNotInherit += weight;
			}
			else if (connect->weightNotInherit < 255)
			{
				connect->weightNotInherit = 255;
			}
		}
		else
		{
			sum = connect->weight * connect->gapLen + gap * weight;
			connect->gapLen = sum / (connect->weight + weight);

			if (!connect->inherit)
			{
				connect->maxSingleWeight = connect->weightNotInherit;
			}

			connect->inherit = 1;
			connect->maxSingleWeight = connect->maxSingleWeight > weight ? connect->maxSingleWeight : weight;
		}

		if (connect->weight + weight <= 255)
		{
			connect->weight += weight;
		}
		else if (connect->weight < 255)
		{
			connect->weight = 255;
		}
	}
	else
	{
		newCntCounter++;
		connect = allocateCN (e2, gap);

		if (cntLookupTable)
		{
			putCnt2LookupTable (e1, connect);
		}

		connect->weight = weight;

		if (contig_array[e1].mask || contig_array[e2].mask)
		{
			connect->mask = 1;
		}

		connect->next = contig_array[e1].downwardConnect;
		contig_array[e1].downwardConnect = connect;

		if (!inherit)
		{
			connect->weightNotInherit = weight;
		}
		else
		{
			connect->weightNotInherit = 0;
			connect->inherit = 1;
			connect->maxSingleWeight = weight;
		}
	}

	return connect;
}

int attach1PE (unsigned int e1, int pre_pos, unsigned int bal_e2, int pos, int insert_size)
{
	int gap, realpeSize;
	unsigned int bal_e1, e2;

	if (e1 == bal_e2)
	{
		ignorePE1++;
		return -1;	//orientation wrong
	}

	bal_e1 = getTwinCtg (e1);
	e2 = getTwinCtg (bal_e2);

	if (e1 == e2)
	{
		realpeSize = contig_array[e1].length + overlaplen - pre_pos - pos;

		if (realpeSize > 0)
		{
			peSUM += realpeSize;
			onsameCtgPE++;

			if ((int) contig_array[e1].length > insert_size)
			{
				int *item = (int *) stackPush (isStack);

				(*item) = realpeSize;
			}
		}

		return 2;
	}

	gap = insert_size - overlaplen + pre_pos + pos - contig_array[e1].length - contig_array[e2].length;

	if (gap < -(insert_size / 10))
	{
		ignorePE2++;
		return 0;
	}

	if (gap > insert_size)
	{
		ignorePE3++;
		return 0;
	}

	add1Connect (e1, e2, gap, 1, 0);
	add1Connect (bal_e2, bal_e1, gap, 1, 0);
	return 1;
}

int connectByPE_grad (gzFile * fp, int peGrad, char *line)
{
	long long pre_readno, readno, minno, maxno;
	int pre_pos, pos, flag, PE, count = 0, Total_PE=0;
	unsigned int pre_contigno, contigno, newIndex;

	if (peGrad < 0 || peGrad > gradsCounter)
	{
		printf ("specified pe grad is out of bound\n");
		return 0;
	}

	maxno = pes[peGrad].PE_bound;

	if (peGrad == 0)
	{
		minno = 0;
	}
	else
	{
		minno = pes[peGrad - 1].PE_bound;
	}

	onsameCtgPE = peSUM = 0;
	PE = pes[peGrad].insertS;

	if (strlen (line))
	{
		sscanf (line, "%lld %d %d", &pre_readno, &pre_contigno, &pre_pos);

		//printf("first record %d %d %d\n",pre_readno,pre_contigno,pre_pos);
		if (pre_readno <= minno)
		{
			pre_readno = -1;
		}
	}
	else
	{
		pre_readno = -1;
	}

	ignorePE1 = ignorePE2 = ignorePE3 = 0;
	static_flag = 1;
	isStack = (STACK *) createStack (CNBLOCKSIZE, sizeof (int));

	while (gzgets (fp, line, lineLen) != NULL)
	{
		sscanf (line, "%lld %d %d", &readno, &contigno, &pos);

		if (readno > maxno)
		{
			break;
		}

		if (readno <= minno)
		{
			continue;
		}

		newIndex = index_array[contigno];

		//if(contig_array[newIndex].bal_edge==0)
		if (isSameAsTwin (newIndex))
		{
			continue;
		}

		if (PE && (readno % 2 == 0) && (pre_readno == readno - 1))	// they are a pair of reads
		{
			Total_PE++;
			flag = attach1PE (pre_contigno, pre_pos, newIndex, pos, PE);

			if (flag == 1)
			{
				count++;
			}
		}

		pre_readno = readno;
		pre_contigno = newIndex;
		pre_pos = pos;
	}
/*
	printf ("%d PEs with insert size %d attached, %d + %d + %d ignored\n", count, PE, ignorePE1, ignorePE2, ignorePE3);
	printf("for insert size: %d\nTotal_PE_link\tSame_contig_right\tSame_contig_wrong\tMinus_dis_PE\tPlus_dis_PE\tCorrect_PE\tAccumulate_connect\n%d\t%d\t%d\t%d\t%d\t%d\t%lld\n",PE,Total_PE,onsameCtgPE,ignorePE1,ignorePE2,ignorePE3,count,newCntCounter);

	if (onsameCtgPE > 0)
	{
		printf ("estimated PE size %lli, by %d pairs\n", peSUM / onsameCtgPE, onsameCtgPE);
	}

	printf ("on contigs longer than %d, %d pairs found,", PE, isStack->item_c);
	printf ("insert_size estimated: %d\n", calcuIS (isStack));
*/
	printf("for insert size: %d\nTotal_PE_link\tSame_contig_right\tSame_contig_wrong\tMinus_dis_PE\tPlus_dis_PE\tCorrect_PE\tAccumulate_connect\n%d\t%d\t%d\t%d\t%d\t%d\t%lld\n",PE,Total_PE,onsameCtgPE,ignorePE1,ignorePE2,ignorePE3,count,newCntCounter);
	printf("Using contigs longer than %d to estimate insert size: \n",PE);
	printf("Pair_num\tSD\tinsert_size\n%d\t",isStack->item_c);
	calcuIS(isStack);
	
	freeStack (isStack);
	return count;
}

static int calcuIS (STACK * intStack)
{
	long long sum = 0;
	int avg = 0;
	int *item;
	int num = intStack->item_c;

	if (num < 100)
	{
		return avg;
	}

	stackBackup (intStack);

	while ((item = (int *) stackPop (intStack)) != NULL)
	{
		sum += *item;
	}

	stackRecover (intStack);
	num = intStack->item_c;
	avg = sum / num;
	sum = 0;
	stackBackup (intStack);

	while ((item = (int *) stackPop (intStack)) != NULL)
	{
		sum += (*item - avg) * (*item - avg);
	}

	int SD = sqrt (sum / (num - 1));

	if (SD == 0)
	{
		printf ("SD=%d, ", SD);
		return avg;
	}

	stackRecover (intStack);
	sum = num = 0;

	while ((item = (int *) stackPop (intStack)) != NULL)
		if (abs (*item - avg) < 3 * SD)
		{
			sum += *item;
			num++;
		}

	if(num == 0) avg = 0;
	else avg = sum / num;
	printf ("SD=%d, ", SD);
	return avg;
}

unsigned int getTwinCtg (unsigned int ctg)
{
	return ctg + contig_array[ctg].bal_edge - 1;
}

boolean isSmallerThanTwin (unsigned int ctg)
{
	return contig_array[ctg].bal_edge > 1;
}

boolean isLargerThanTwin (unsigned int ctg)
{
	return contig_array[ctg].bal_edge < 1;
}

boolean isSameAsTwin (unsigned int ctg)
{
	return contig_array[ctg].bal_edge == 1;
}
