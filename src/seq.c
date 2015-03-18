/*
 * Seq.c
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
/*
void print_kmer(FILE *fp,Kmer kmer,char c)
{
    fprintf(fp,"%llx %llx %llx %llx",kmer.high1,kmer.low1,kmer.high2,kmer.low2);
    fprintf(fp,"%c",c);

}*/
void printTightString (char *tightSeq, int len)
{
	int i;

	for (i = 0; i < len; i++)
	{
		printf ("%c", int2base ((int) getCharInTightString (tightSeq, i)));
		if ((i + 1) % 100 == 0)
		{
			printf ("\n");
		}
	}
	printf ("\n");
}

void writeChar2tightString (char nt, char *tightSeq, int pos)
{
	char *byte = tightSeq + pos / 4;

	switch (pos % 4)
	{
	case 0:
		*byte &= 63;
		*byte += nt << 6;
		return;
	case 1:
		*byte &= 207;
		*byte += nt << 4;
		return;
	case 2:
		*byte &= 243;
		*byte += nt << 2;
		return;
	case 3:
		*byte &= 252;
		*byte += nt;
		return;
	}
}

char getCharInTightString (char *tightSeq, int pos)
{
	char *byte = tightSeq + pos / 4;

	switch (pos % 4)
	{
	case 3:
		return (*byte & 3);
	case 2:
		return (*byte & 12) >> 2;
	case 1:
		return (*byte & 48) >> 4;
	case 0:
		return (*byte & 192) >> 6;
	}
	return 0;
}

// complement of sequence denoted 0, 1, 2, 3
void reverseComplementSeq (char *seq, int len, char *bal_seq)
{
	int i, index = 0;

	if (len < 1)
	{
		return;
	}
	for (i = len - 1; i >= 0; i--)
	{	
		//mao 2011 10 8
	//	if(seq[i] >= 4)
	//		bal_seq[index++]=seq[i];
	//	else
			bal_seq[index++] = int_comp (seq[i]);
	}
	return;
}

// complement of sequence denoted 0, 1, 2, 3
char *compl_int_seq (char *seq, int len)
{
	char *bal_seq = NULL, c, bal_c;
	int i, index;

	if (len < 1)
	{
		return bal_seq;
	}
	bal_seq = (char *) ckalloc (len * sizeof (char));
	index = 0;
	for (i = len - 1; i >= 0; i--)
	{
		c = seq[i];
		if (c < 4)
		{
			bal_c = int_comp (c);
		}		//3-c;
		else
		{
			bal_c = c;
		}
		bal_seq[index++] = bal_c;
	}
	return bal_seq;
}

long long trans_seq (char *seq, int len)
{
	int i;
	long long res;

	res = 0;
	for (i = 0; i < len; i++)
	{
		res = res * 4 + seq[i];
	}
	return (res);
}

/*
char *kmer2seq(Kmer word)
{
    int i;
    char *seq;
    Kmer charMask = 3;

    seq = (char *)ckalloc(overlaplen*sizeof(char));
    for(i=overlaplen-1;i>=0;i--){
        seq[i] = charMask&word;
        word >>= 2;
    }
    return seq;
}
*/
