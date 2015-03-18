/*
 * cutTip_graph.c
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

static int caseA, caseB, caseC, caseD, caseE;

static int  removeArc()
{
	int count=0;
	unsigned int index;
	ARC * arc,*arc_temp;
	for(index=1;index<=num_ed;index++)
	{
		arc=edge_array[index].arcs;
		while(arc)
		{
			arc_temp=arc;
			arc=arc->next;

			if(arc_temp->multiplicity == 0)
			{
				if(arc_temp->prev)
				{
					arc_temp->prev->next=arc_temp->next;					
				}
				else
				{
					edge_array[index].arcs=arc_temp->next;	
				}
				if(arc_temp->next)
				{
					arc_temp->next->prev = arc_temp->prev;
				}
				dismissArc (arc_temp);
				count++;				
			}
		}		
	}	
	return count;
}
void destroyEdge (unsigned int edgeid)
{
	unsigned int bal_ed = getTwinEdge (edgeid);
	ARC *arc;

	if (bal_ed == edgeid)
	{
		edge_array[edgeid].length = 0;
		return;
	}

	arc = edge_array[edgeid].arcs;

	while (arc)
	{
		arc->bal_arc->to_ed = 0;
		arc = arc->next;
	}

	arc = edge_array[bal_ed].arcs;

	while (arc)
	{
		arc->bal_arc->to_ed = 0;
		arc = arc->next;
	}

	edge_array[edgeid].arcs = NULL;
	edge_array[bal_ed].arcs = NULL;
	edge_array[edgeid].length = 0;
	edge_array[bal_ed].length = 0;
	edge_array[edgeid].deleted = 1;
	edge_array[bal_ed].deleted = 1;
	//printf("Destroyed %d and %d\n",edgeid,bal_ed);
}

ARC *arcCount (unsigned int edgeid, unsigned int *num)
{
	ARC *arc;
	ARC *firstValidArc = NULL;
	unsigned int count = 0;

	arc = edge_array[edgeid].arcs;

	while (arc)
	{
		if (arc->to_ed > 0)
		{
			count++;

			if (count == 1)
			{
				firstValidArc = arc;
			}
			else if (count > 1)
			{
				*num = count;
				return firstValidArc;
			}
		}

		arc = arc->next;
	}

	*num = count;
	return firstValidArc;
}

void removeWeakEdges (int lenCutoff, unsigned int multiCutoff)
{
	unsigned int bal_ed;
	unsigned int arcRight_n, arcLeft_n;
	ARC *arcLeft, *arcRight;
	unsigned int i;
	int counter = 0;

	for (i = 1; i <= num_ed; i++)
	{
		if (edge_array[i].deleted || edge_array[i].length == 0 || edge_array[i].length > lenCutoff || EdSameAsTwin (i))
		{
			continue;
		}

		bal_ed = getTwinEdge (i);
		arcRight = arcCount (i, &arcRight_n);

		if (arcRight_n > 1 || !arcRight || arcRight->multiplicity > multiCutoff)
		{
			continue;
		}

		arcLeft = arcCount (bal_ed, &arcLeft_n);

		if (arcLeft_n > 1 || !arcLeft || arcLeft->multiplicity > multiCutoff)
		{
			continue;
		}

		destroyEdge (i);
		counter++;
	}

	printf ("%d weak inner edges destroyed\n", counter);
	removeDeadArcs ();
	/*
	   linearConcatenate();
	   compactEdgeArray();
	 */
}

void removeLowCovEdges (int lenCutoff, unsigned short covCutoff)
{
	unsigned int bal_ed;
	unsigned int arcRight_n, arcLeft_n;
	ARC *arcLeft, *arcRight;
	unsigned int i;
	int counter = 0;

	for (i = 1; i <= num_ed; i++)
	{
		if (edge_array[i].deleted || edge_array[i].cvg == 0 || edge_array[i].cvg > covCutoff * 10 || edge_array[i].length >= lenCutoff || EdSameAsTwin (i) || edge_array[i].length == 0)
		{
			continue;
		}

		bal_ed = getTwinEdge (i);
		arcRight = arcCount (i, &arcRight_n);
		arcLeft = arcCount (bal_ed, &arcLeft_n);

		if (arcLeft_n < 1 || arcRight_n < 1)
		{
			continue;
		}

		destroyEdge (i);
		counter++;
	}

	printf ("Remove low coverage(%d): %d inner edges destroyed\n", covCutoff, counter);
	removeDeadArcs ();
	linearConcatenate ();
	compactEdgeArray ();
}

boolean isUnreliableTip (unsigned int edgeid, int cutLen, boolean strict)
{
	unsigned int arcRight_n, arcLeft_n;
	unsigned int bal_ed;
	unsigned int currentEd = edgeid;
	int length = 0;
	unsigned int mult = 0;
	ARC *arc, *activeArc = NULL, *tempArc;

	if (edgeid == 0)
	{
		return 0;
	}

	bal_ed = getTwinEdge (edgeid);

	if (bal_ed == edgeid)
	{
		return 0;
	}

	arcCount (bal_ed, &arcLeft_n);

	if (arcLeft_n > 0)
	{
		return 0;
	}

	while (currentEd)
	{
		arcCount (bal_ed, &arcLeft_n);
		tempArc = arcCount (currentEd, &arcRight_n);

		if (arcLeft_n > 1 || arcRight_n > 1)
		{
			break;
		}

		length += edge_array[currentEd].length;

		if (tempArc)
		{
			activeArc = tempArc;
			currentEd = activeArc->to_ed;
			bal_ed = getTwinEdge (currentEd);
		}
		else
		{
			currentEd = 0;
		}
	}

	if (length >= cutLen)
	{
		return 0;
	}

	if (currentEd == 0)
	{
		caseB++;
		return 1;
	}

	if (!strict)
	{
		if (arcLeft_n < 2)
		{
			length += edge_array[currentEd].length;
		}

		if (length >= cutLen)
		{
			return 0;
		}
		else
		{
			caseC++;
			return 1;
		}
	}

	if (arcLeft_n < 2)
	{
		return 0;
	}

	if (!activeArc)
	{
		printf ("no activeArc while checking edge %d\n", edgeid);
	}

	if (activeArc->multiplicity == 1)
	{
		caseD++;
		return 1;
	}

	for (arc = edge_array[bal_ed].arcs; arc != NULL; arc = arc->next)
		if (arc->multiplicity > mult)
		{
			mult = arc->multiplicity;
		}

	if (mult > activeArc->multiplicity)
	{
		caseE++;
	}

	return mult > activeArc->multiplicity;
}

boolean isUnreliableTip_strict (unsigned int edgeid, int cutLen)
{
	unsigned int arcRight_n, arcLeft_n;
	unsigned int bal_ed;
	unsigned int currentEd = edgeid;
	int length = 0;
	unsigned int mult = 0;
	ARC *arc, *activeArc = NULL, *tempArc;

	if (edgeid == 0)
	{
		return 0;
	}

	bal_ed = getTwinEdge (edgeid);

	if (bal_ed == edgeid)
	{
		return 0;
	}

	arcCount (bal_ed, &arcLeft_n);

	if (arcLeft_n > 0)
	{
		return 0;
	}

	while (currentEd)
	{
		arcCount (bal_ed, &arcLeft_n);
		tempArc = arcCount (currentEd, &arcRight_n);

		if (arcLeft_n > 1 || arcRight_n > 1)
		{
			if (arcLeft_n == 0 || length == 0)
			{
				return 0;
			}
			else
			{
				break;
			}
		}

		length += edge_array[currentEd].length;

		if (length >= cutLen)
		{
			return 0;
		}

		if (tempArc)
		{
			activeArc = tempArc;
			currentEd = activeArc->to_ed;
			bal_ed = getTwinEdge (currentEd);
		}
		else
		{
			currentEd = 0;
		}
	}

	if (currentEd == 0)
	{
		caseA++;
		return 1;
	}

	if (!activeArc)
	{
		printf ("no activeArc while checking edge %d\n", edgeid);
	}

	if (activeArc->multiplicity == 1)
	{
		caseB++;
		return 1;
	}

	for (arc = edge_array[bal_ed].arcs; arc != NULL; arc = arc->next)
		if (arc->multiplicity > mult)
		{
			mult = arc->multiplicity;
		}

	if (mult > activeArc->multiplicity)
	{
		caseC++;
	}

	return mult > activeArc->multiplicity;
}

void removeDeadArcs ()
{
	unsigned int i, count = 0;
	ARC *arc, *arc_temp;

	for (i = 1; i <= num_ed; i++)
	{
		arc = edge_array[i].arcs;

		while (arc)
		{
			arc_temp = arc;
			arc = arc->next;

			if (arc_temp->to_ed == 0)
			{
				count++;
				edge_array[i].arcs = deleteArc (edge_array[i].arcs, arc_temp);
			}
		}
	}

	printf ("%d dead arcs removed\n", count);
}

void cutTipsInGraph (int cutLen, boolean strict)
{
	int flag = 1;
	unsigned int i;

	if (!cutLen)
	{
		cutLen = 2 * overlaplen;
	}

	//if(cutLen > 100) cutLen = 100;
	printf ("strict %d, cutLen %d\n", strict, cutLen);

	if (strict)
	{
		linearConcatenate ();
	}

	caseA = caseB = caseC = caseD = caseE = 0;

	while (flag)
	{
		flag = 0;

		for (i = 1; i <= num_ed; i++)
		{
			if (edge_array[i].deleted)
			{
				continue;
			}

			if (isUnreliableTip (i, cutLen, strict))
			{
				destroyEdge (i);
				flag++;
			}
		}

		printf ("a cutTipsInGraph lap, %d tips cut\n", flag);
	}

	removeDeadArcs ();

	if (strict)
	{
		printf ("case A %d, B %d C %d D %d E %d\n", caseA, caseB, caseC, caseD, caseE);
	}

	linearConcatenate ();
	compactEdgeArray ();
	printf("\n");
}
void delowHighArc(int multi)
{
	unsigned int i, twin,to_edge,count = 0;
	ARC *arc, *arc_temp;
	unsigned int in_weight,out_weight,curr_weight;

	for (i = 1; i <= num_ed; i++)
	{
		in_weight=0;		
		curr_weight=0;

		//获取i的in_flow权重
		twin=getTwinEdge(i);
		arc=edge_array[twin].arcs;
		while(arc)
		{
			in_weight += arc->multiplicity;
			arc=arc->next;
		}
		
		arc = edge_array[i].arcs;		
		while (arc)
		{				
			curr_weight = arc->multiplicity;

			to_edge = arc->to_ed;
			arc_temp = edge_array[to_edge].arcs;
			out_weight=0;
			while(arc_temp)
			{
				out_weight += arc_temp->multiplicity;
				arc_temp=arc_temp->next;
			}
			if( in_weight != 0 && curr_weight !=0 && curr_weight > in_weight*multi && curr_weight > out_weight*multi)
			{
				count++;
				arc->multiplicity= in_weight > out_weight ? in_weight : out_weight;
			}
			arc=arc->next;
			
		}		
	}	
//	printf("delow arc : %d\n",count);
}
static unsigned int deleteLightOutArc(double min_arc_rate)
{
	unsigned int index,total_weight,count=0;
	ARC * arc,*next_arc,*twin_arc;
	unsigned int to_ed,twin_te;
	for(index=1;index<=num_ed;index++)
	{
		total_weight=0;
		arc = edge_array[index].arcs;
		while(arc)
		{
			total_weight += arc->multiplicity;
			arc=arc->next;
		}

		if(total_weight > 0)
		{
			arc = edge_array[index].arcs;
			while(arc)
			{
				next_arc = arc->next;
				to_ed = arc ->to_ed;
				if(arc->multiplicity != 0 && arc->multiplicity <= total_weight*min_arc_rate)
				{					
					twin_arc=arc->bal_arc;
					arc->multiplicity=0;
					twin_arc->multiplicity=0;
//					twin_te=getTwinEdge(to_ed);
//					edge_array[index].arcs=deleteArc (edge_array[index].arcs, arc);
//					edge_array[twin_te].arcs=deleteArc(edge_array[twin_te].arcs,twin_arc);
					count++;
				}
				arc=next_arc;
			}		
		}
//		printf("%d\t%d\n",index,edge_array[5998].arcs->to_ed);
		
	}
	return count;	
}
static unsigned int deleteLightFlowArc(double min_arc_rate)
{
	unsigned int index,twin,count=0;
	unsigned int total_in_weight,total_out_weight,coverage;
	ARC * arc,*next_arc,*twin_arc;
	unsigned int to_ed,twin_te;
	for(index=1;index<=num_ed;index++)
	{
		total_in_weight=0;
		total_out_weight=0;
		twin = getTwinEdge(index);
		coverage = (double)edge_array[index].cvg/10;

		arc = edge_array[index].arcs;
		while(arc)
		{
			total_out_weight += arc->multiplicity;
			arc=arc->next;
		}

		arc = edge_array[twin].arcs;
		while(arc)
		{
			total_in_weight += arc->multiplicity;
			arc=arc->next;
		}

		arc = edge_array[index].arcs;
		while(arc)
		{
			next_arc = arc->next;
			to_ed = arc ->to_ed;
			if(arc->multiplicity != 0 && arc->multiplicity <= (double)total_in_weight*min_arc_rate || arc->multiplicity <= (double)coverage*min_arc_rate)
			{
				twin_arc=arc->bal_arc;
				arc->multiplicity=0;
				twin_arc->multiplicity=0;
				count++;
			}
			arc=next_arc;
		}
		arc = edge_array[twin].arcs;
		while(arc)
		{
			next_arc = arc->next;
			to_ed = arc ->to_ed;
			if(arc->multiplicity != 0 && arc->multiplicity <= (double)total_out_weight*min_arc_rate || arc->multiplicity <= (double)coverage*min_arc_rate)
			{
				twin_arc=arc->bal_arc;
				arc->multiplicity=0;
				twin_arc->multiplicity=0;
				count++;
			}
			arc=next_arc;
		}
		if(twin != index)
			index++;
	}
	return count;	
}
int deleteLightArc()
{
	printf("Start to remove light arc using the rate %f and %f\n",(double)da/100,(double)dA/100);
	unsigned int changed = deleteLightOutArc((double)da/100);
	unsigned int flow_changed = deleteLightFlowArc((double)dA/100);
/*
	unsigned int index,total_changed=0;	
	ARC * arc ,*arc_temp;
	for(index=1;index<=num_ed;index++)
	{
		arc=edge_array[index].arcs;
		while(arc)
		{
			arc_temp=arc;
			arc=arc->next;

			if(arc_temp->multiplicity == 0)
			{
				if(arc_temp->prev)
				{
					arc_temp->prev->next=arc_temp->next;					
				}
				else
				{
					edge_array[index].arcs=arc_temp->next;	
				}
				if(arc_temp->next)
				{
					arc_temp->next->prev = arc_temp->prev;
				}
				dismissArc (arc_temp);
				total_changed++;
			}
		}		
	}*/
	int total_changed = removeArc();
	printf("delete light arc :%d/%d\n",total_changed/2,changed+flow_changed);
	return (changed+flow_changed)>0 ? 1 : 0;
}
void deleteUnlikeArc()
{
	unsigned int index;
	ARC * arc ,*arc_temp;
	unsigned short source_multi,target_multi,max_multi;
	for(index=1;index<=num_ed;index++)
	{
		arc=edge_array[index].arcs;
		source_multi=edge_array[index].cvg;		
		while(arc)
		{
			target_multi=edge_array[arc->to_ed].cvg;
			max_multi=source_multi>target_multi?source_multi:target_multi;
//			printf("test:\t%4.1f",arc->multiplicity/max_multi);
			if(arc->multiplicity <  (double)max_multi/25 ||arc->multiplicity<3 )
			{
				arc->multiplicity=0;
			}			
			arc=arc->next;
		}
	}
	/*
	for(index=1;index<=num_ed;index++)
	{
		arc=edge_array[index].arcs;
		while(arc)
		{
			arc_temp=arc;
			arc=arc->next;

			if(arc_temp->multiplicity == 0)
			{
				if(arc_temp->prev)
				{
					arc_temp->prev->next=arc_temp->next;					
				}
				else
				{
					edge_array[index].arcs=arc_temp->next;	
				}
				if(arc_temp->next)
				{
					arc_temp->next->prev = arc_temp->prev;
				}
				dismissArc (arc_temp);
				total_changed++;
			}
		}		
	}*/
	int total_changed = removeArc();
	printf("delete unlike arc :%d\n",total_changed);
}
static void computeNextCov(unsigned int contigID,double * cov )
{
	*cov=0;
	int length=0;
	ARC *arc = edge_array[contigID].arcs;
	while(arc)
	{
		if(!(edge_array[arc->to_ed].deleted))
		{
			*cov += (double)(edge_array[arc->to_ed].cvg * edge_array[arc->to_ed].length);
			length += edge_array[arc->to_ed].length;
		}
		arc=arc->next;
	}
	if(length >0 )
		*cov /= length;
	else
		*cov =0;
}
static void delete1contig(unsigned int edgeid)
{
	edge_array[edgeid].cvg=0;
	edge_array[edgeid].deleted=1;
	edge_array[edgeid].length=0;

	ARC *arc=edge_array[edgeid].arcs;
	while(arc)
	{
		arc->multiplicity=0;
		arc->bal_arc->multiplicity=0;
		arc=arc->next;
	}

	if(EdSameAsTwin(edgeid))
		return;
	edge_array[getTwinEdge(edgeid)].cvg=0;
	edge_array[getTwinEdge(edgeid)].deleted=1;
	edge_array[getTwinEdge(edgeid)].length=0;
	arc = edge_array[getTwinEdge(edgeid)].arcs;
	while(arc)
	{
		arc->multiplicity=0;
		arc->bal_arc->multiplicity=0;
		arc=arc->next;
	}
}
int deleteLightContig()
{
	double prev_cov,next_cov,max,min,curr_cov;
	unsigned int index;
	int change=0;
	ARC * arc,*arc_temp;

	for(index=1;index<=num_ed;index++)
	{		
		if(EdSameAsTwin(index))
			continue;
		computeNextCov(index,&next_cov);
		computeNextCov(getTwinEdge(index),&prev_cov);
		if(next_cov ==0 || prev_cov ==0)
			continue;
		if(next_cov > prev_cov)
		{
			max=next_cov;
			min=prev_cov;
		}
		else
		{
			max=prev_cov;
			min=next_cov;
		}
		curr_cov = (double)edge_array[index].cvg;
		printf("contig_cov:\t%0.1f\t%0.1f\t%0.1f\n",curr_cov,max,min);


		if(min / max <0.1)
		{
			if(curr_cov /min < 0.5)
			{				
				delete1contig(index);
			}
		}
		else
		{
			if(curr_cov / max <0.05)
			{				
				delete1contig(index);
			}
		}	
		index++;
	}
	/*
	for(index=1;index<=num_ed;index++)
	{
		arc=edge_array[index].arcs;
		while(arc)
		{
			arc_temp=arc;
			arc=arc->next;

			if(arc_temp->multiplicity == 0)
			{
				if(arc_temp->prev)
				{
					arc_temp->prev->next=arc_temp->next;					
				}
				else
				{
					edge_array[index].arcs=arc_temp->next;	
				}
				if(arc_temp->next)
				{
					arc_temp->next->prev = arc_temp->prev;
				}
				dismissArc (arc_temp);
				change++;
			}
		}		
	}*/
	change = removeArc();
	return change>0?1:0;;
}

static int extern_contig(unsigned int edgeid,int pool_index)
{
	if(pool[edgeid]!=0)
		return 0;
	pool[edgeid]=pool_index;	
	pool[ getTwinEdge(edgeid)]=pool_index;
	int length=0;
	length += edge_array[edgeid].length;
	ARC *arc;
	unsigned int best_id;
	int max_arc;
	unsigned int curr_edge= edgeid;

	while(curr_edge)
	{
		max_arc=0;
		arc = edge_array[edgeid].arcs;
		while(arc)
		{
			if(pool[arc->to_ed] ==0)
			{
				if(arc->multiplicity > max_arc)
				{
					max_arc=arc->multiplicity;
					best_id=arc->to_ed;
				}
			}
			arc=arc->next;
		}
		if(max_arc>0)
		{			
			pool[best_id]=pool_index;
			pool[getTwinEdge(best_id)]=pool_index;
			length += edge_array[best_id].length;
			curr_edge=best_id;
		}
		else
			curr_edge=0;
	}

	curr_edge= getTwinEdge(edgeid);
	

	while(curr_edge)
	{
		max_arc=0;
		arc = edge_array[edgeid].arcs;
		while(arc)
		{
			if(pool[arc->to_ed] ==0)
			{
				if(arc->multiplicity > max_arc)
				{
					max_arc=arc->multiplicity;
					best_id=arc->to_ed;
				}
			}
			arc=arc->next;
		}
		if(max_arc>0)
		{			
			pool[best_id]=pool_index;
			pool[getTwinEdge(best_id)]=pool_index;
			length += edge_array[best_id].length;
			curr_edge=best_id;
		}
		else
			curr_edge=0;
	}
	return length;
}
typedef struct cov_list
{
	unsigned int contig;
	unsigned short cov;
}COV_LIST;

static int cmp_cov (const void *a, const void *b) 
{
	COV_LIST *A, *B;
	
	A = (COV_LIST *) a;
	B = (COV_LIST *) b;
	if (A->cov > B->cov)		
	{
		return -1;
	}
	
	else if (A->cov == B->cov)		
	{
		return 0;
	}
	
	else		
	{
		return 1;
	}
}
void deleteShortContig(int cutLength)
{
	unsigned int index;
	if(pool== NULL)
		pool= (int*)ckalloc (sizeof(int)*(num_ed+1));
	int * poolid_length=(int*)ckalloc(sizeof(int)*(num_ed+1));
	for(index=0;index<=num_ed;index++)
	{
		pool[index]=0;
		poolid_length[index]=0;
	}
	int poolid_index=1;

	COV_LIST * cov = (COV_LIST * ) ckalloc (sizeof(COV_LIST)*(num_ed+1));
	for(index=1;index<=num_ed;index++)
	{
		cov[index].contig=index;
		cov[index].cov=edge_array[index].cvg;
	}
	qsort(&cov[1], num_ed, sizeof(COV_LIST), cmp_cov);	
	for(index=1;index<=num_ed;index++)
	{		
		poolid_length[poolid_index]=extern_contig(cov[index].contig,poolid_index);
		if(poolid_length[poolid_index]!=0)
			poolid_index++;		
	}

	int num_delelte=0;
	for(index=1;index<=num_ed;index++)
	{
		if(poolid_length[pool[index]]<cutLength)
		{	
			delete1contig(index);
			num_delelte++;
		}
		if(!EdSameAsTwin(index))
			index++;
	}
	
	free(poolid_length);
	free(pool);
	free(cov);
	printf("%d short contig(<%d) removed \n",num_delelte,cutLength);
	
	removeArc();
}
void deleteWeakEdge(unsigned short cutoff)
{
	if(cutoff > 30)
		cutoff=30;
	printf("Start to remove the low coverage edge < %d\n",cutoff/10);
	unsigned int index;
	int total=0;
	for(index=1;index<=num_ed;index++)
	{
		if(edge_array[index].cvg < cutoff)
		{
			delete1contig(index);
			total++;
		}
		if(!EdSameAsTwin(index))
			index++;
	}
	printf("%d edges removed\n\n",total);
	removeArc();
}
void resetCov()
{
	unsigned int index;
	if(pool== NULL)
		pool= (int*)ckalloc (sizeof(int)*(num_ed+1));
	int * poolid_length=(int*)ckalloc(sizeof(int)*(num_ed+1));
	for(index=0;index<=num_ed;index++)
	{
		pool[index]=0;
		poolid_length[index]=0;
	}
	int poolid_index=1;

	COV_LIST * cov = (COV_LIST * ) ckalloc (sizeof(COV_LIST)*(num_ed+1));
	for(index=1;index<=num_ed;index++)
	{
		cov[index].contig=index;
		cov[index].cov=edge_array[index].cvg;
	}
	qsort(&cov[1], num_ed, sizeof(COV_LIST), cmp_cov);
/*
	for(index=1;index<=num_ed;index++)
	{
		printf("cov:\t%d\n",cov[index].cov);
	}
*/
	
	for(index=1;index<=num_ed;index++)
	{		
		poolid_length[poolid_index]=extern_contig(cov[index].contig,poolid_index);
		if(poolid_length[poolid_index]!=0)
			poolid_index++;		
	}
	int i;
	unsigned int *contig_cov = (unsigned int *) ckalloc ( sizeof(unsigned int ) * (poolid_index + 1));
	unsigned int *contig_length = (unsigned int *) ckalloc ( sizeof(unsigned int ) * (poolid_index + 1));
	for(i=1;i<poolid_index;i++)
	{
		contig_cov[i]=0;
		contig_length[i]=0;
	}
	for(index=1;index<=num_ed;index++)
	{
		contig_cov[pool[index]]+=edge_array[index].cvg *edge_array[index].length;
		contig_length[pool[index]]+=edge_array[index].length;
		if(!(EdSameAsTwin(index)))		
			index++;		
	}
	for(i=1;i<poolid_index;i++)
	{		
		if(contig_length[i]>0)
			contig_cov[i] /= contig_length[i];
		else
			printf("pool length == 0\n");
	}
	for(index=1;index<=num_ed;index++)
	{
		edge_array[index].cvg = contig_cov[pool[index]];
	}
	free(cov);
	free(pool);
	pool=NULL;
	free(poolid_length);
	free(contig_cov);
	free(contig_length);
}
void deleteSimpleLoop()
{
	unsigned int index;
	int loop1=0,loop2=0;
	for(index=1;index<=num_ed;index++)
	{
		ARC * arc1 = getArcBetween(index,index);
		if(arc1)
		{
			arc1->multiplicity=0;
			arc1->bal_arc->multiplicity=0;
			loop1++;
		}
		ARC * arc = edge_array[index].arcs;
		
		while(arc)
		{
			unsigned int to_ed = arc->to_ed;			
			ARC * arc2 = getArcBetween(to_ed,index);
			if(arc2)
			{
				arc2->multiplicity=0;
				arc2->bal_arc->multiplicity=0;
				arc->multiplicity=0;
				arc->bal_arc->multiplicity=0;
				loop2++;
			}
			arc=arc->next;				
		}
	}
	printf("%d loops removed in Graph\n",loop1+loop2);
	removeArc();
}
