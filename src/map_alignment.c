// This file contains all functions for the fine tuned aligning of preplaced reads
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "referenceAssembly.h"
#include "ref_io.h"
#include "ref_slidingWindow.h"
#include "ref_math.h"
#include "map_alignment.h"


// void locallyAlignSequencesGtoh(char * seq1, unsigned int length1, char * seq2, unsigned int length2)
// {
//  int mismatch_penalty = 1;
//  int gap_penalty = 1;
//     
//     
//     
//     
//     
//     
//     
//     
// }

void printQuadraticMatrix(double** matrix, int dim){
    
    int i,j;
    
    for(i=0; i<dim; i++)
    {
        for(j=0; j<dim; j++)
        {
        print_selective("%.2f\t",matrix[i][j]);
        }        
    print_selective("\n");
    }
    
    
}

void printQuadraticMatrixInt(int** matrix, int dim){
    
    int i,j;
    
    for(i=0; i<dim; i++)
    {
        for(j=0; j<dim; j++)
        {
        print_selective("%d\t",matrix[i][j]);
        }        
    print_selective("\n");
    }
    
    
}

void printCirseqMatrix(int* jump_to, char* seq, int seq_length, int starting_i){
 int i,j;
 
 for(i = 0; i < seq_length; i++){
 for(j = 0; j < seq_length; j++){
     printf("\t%d", jump_to[i+j*seq_length]);
 }
     printf("\n");
 }
 
 i=starting_i;
 for(j = seq_length-1; j>=0; j--){
  printf("%c", seq[j]);
  
  if(jump_to[i+j*seq_length] >0){
   i = jump_to[i+j*seq_length];
   printf("\n");
  }else{
      i--;
      
}
 }
}

void printReplicates(char * seq, int** jump_to, int dim, int start){
    
    int i;
    int pos = start;
    for(i = dim-1 ; i >= 0 ; i--){
        print_selective("%c",seq[i]);
        
//         print_selective("%d %d\n", pos, i);
        if(jump_to[pos][i] > 0){
         print_selective("\n");
         pos = jump_to[pos][i];
        }else{
         pos--;   
        }
    }
     print_selective("\n");
    
}

double unadjustedScore(char * seq, int** jump_to, int dim, int start, double match, double miss){
    
    int i;
    int pos = start;
    
    double score = 0.0;
    
    for(i = dim-1 ; i >= 0 ; i--){
        
        if(seq[i]==seq[pos])
            score+=match;
        else
            score+=miss;
        
        
        if(jump_to[pos][i] > 0){
         pos = jump_to[pos][i];
        }else{
         pos--;   
        }
    }
//      print_selective("\n");
    return score;
}


int findNumberOfReplicates(char * seq, int** jump_to, int dim, int start, int estimated_length){
    
    assert(estimated_length > 0);
    
    int i;
    int pos = start;
    int number = 0;
    
    int * num_hits = (int*)calloc(dim, sizeof(int));
    
    for(i = dim-1 ; i >= 0 ; i--){
        
        num_hits[pos]++;
        if(num_hits[pos] > number)
        {
            number = num_hits[pos];
        }
        
        if(jump_to[pos][i] > 0){
         pos = jump_to[pos][i];
        }else{
         pos--;   
        }
    }

    int sum;
    
    while(number > 1)
    {
        sum = 0;
        for(i = 0 ; i < dim ; i++ )
        {
            if(num_hits[i] > number)
            {
            sum++;   
            }
        
        }
        
        if(sum/(double)estimated_length >= 0.8)
            break;
        
        number--;
    }
    
    free(num_hits);
    
//     if(number <= 1 || estimated_length == 0 ||  sum/(double)estimated_length < 0.8)
    if(number <= 1)
    {
        return 1;
    }else{
        return number;
    }
}


int findStartReplicates(char * seq, int** jump_to, int dim, int start){
    
    int i;
    int pos = start;
    int start_pos = 0;
    
    int * num_hits = (int*)calloc(dim, sizeof(int));
    
    for(i = dim-1 ; i >= 0 ; i--){
        
        num_hits[pos]++;
        
        if(jump_to[pos][i] > 0){
         pos = jump_to[pos][i];
        }else{
         pos--;   
        }
    }
//      print_selective("\n");
// for(i = 0 ; i < dim ; i++)
// {
//     print_selective("%d ", num_hits[i]);
// }
// print_selective("\n");

    
    for(i = 0 ; i < dim ; i++ )
    {
        start_pos = i;
        
        if(num_hits[i] > 1)
        {
         break;   
        }
        
        
    }
    
//     print_selective("start_pos %d\n", start_pos);
    
    free(num_hits);
    
    if(start_pos >= dim-1)
    {
        return 0;
    }else{
        return start_pos;
    }
}



void calculateJumpMatrix(char *seq, unsigned int seq_length, double match, double miss, double jump, double** score, int** jump_to, double* max, int* which_max)
{
     int lane, i, j;
        
        double multiplier = 1.0;
        
        double add = miss;
//     Top row initialization
    i = 0;
    for(j=0; j<seq_length; j++)
    {
         multiplier = 1.0 + 2*j/(double)seq_length;
        if(mapIgnoreDegenerate(seq[i]) == mapIgnoreDegenerate('N') || mapIgnoreDegenerate(seq[j]) == mapIgnoreDegenerate('N'))
        {
         add = 0;   
        }else if(seq[i] == seq[j])
        {
        add = match*multiplier; 
        }else{
        add = miss*multiplier;
        }
        
        score[i][j] = add;
    }
    
// Calculating diagonal movement only     
    for(i=1; i<seq_length; i++)
    {
//         Start with diagonal
  
        for(j=i; j<seq_length; j++)
        {
         lane = j-i;
         assert(lane >= 0);
         
         multiplier = 1.0 + 2*lane/(double)seq_length;
         
        if(mapIgnoreDegenerate(seq[i]) == mapIgnoreDegenerate('N') || mapIgnoreDegenerate(seq[j]) == mapIgnoreDegenerate('N'))
        {
         add = 0;   
        }else if(seq[i] == seq[j])
            {
                add = match*multiplier; 
            }else{
                add = miss*multiplier;
            }
        
          score[i][j] = score[i-1][j-1]+add;
        }  
    }
 
//      print_selective("%s\n\n",seq);
//     printQuadraticMatrix(score, seq_length);
    
    
//     Postprocessing for vertical jumps by looking at each column seperately
// Diagonal: 
    for(i=0; i<seq_length; i++)
    {
    max[i] = score[i][i];    
    which_max[i] = i;
    }
    
// The rest:    
    for(lane=1; lane<seq_length; lane++)
    {
        
        multiplier = 1.0 + 2*lane/(double)seq_length;
                 
        for(i=0; i<seq_length-lane; i++)
        {

         j = i + lane;   
         
         

        if(mapIgnoreDegenerate(seq[i]) == mapIgnoreDegenerate('N') || mapIgnoreDegenerate(seq[j]) == mapIgnoreDegenerate('N'))
        {
         add = 0;   
        }else if(seq[i] == seq[j])
            {
                add = match*multiplier; 
            }else{
                add = miss*multiplier;
            }
         

        if(i>0)
        {          
            if(score[i-1][j-1] > max[j-1] + jump){
                score[i][j] = score[i-1][j-1] + add;   
            }else{
                score[i][j] =   max[j-1] + add + jump;
                jump_to[i][j] = which_max[j-1];
            }

         }else{//if i == 0
            
            if(0 > max[j-1] + jump){
                score[i][j] = 0 + add;   
            }else{
                score[i][j] =   max[j-1] + add + jump;
                jump_to[i][j] = which_max[j-1];
            }
             
         }
         
         
         if(score[i][j] > max[j])
         {
             max[j] = score[i][j];
             which_max[j] = i;
         }

        }
//     print_selective("%d $d\n",i,j);
//     fflush(stderr);
    }
}

// TODO must be much refined to actually take qualityvalues. This is just a hack to have something to work with
// should be in ref_slidingWindow instead
char combine(char c, char b){
//     this works only for 4 base DNA encoding. 0010 & 0011 becomes 0010
    int temp = map(c) & map(b);
    
    if(temp)
    {
        return reMap(temp);
    }else{
//         or if there is no overlap, take the union instead 0100 | 0011 becomes 0111
        temp = map(c) | map(b);
        return reMap(temp);
    }
    
}


void findConsensusReplicates(char* result_seq, double* result_phred, char * seq, char* seqQ, int** jump_to, int dim, int start, int start_limit, int end_limit)
{
            int i;
    int pos = start;
    
    
    for(i = dim-1 ; i >= start_limit ; i--){
        
        if(i < end_limit)
        {
//          if we already encountered a hit here
            if(result_seq[pos]){
//                 double temp_phred;
            result_seq[pos] = combine(result_seq[pos], seq[i]);
//             result_seq[pos] = combine_exact(&temp_phred, result_seq[pos], result_phred[pos], seq[i], cQ2P(seqQ[i]));
//             result_phred[pos] = temp_phred;
            
            }else{//If this is the first hit for the position
                result_seq[pos] = seq[i];
                result_phred[pos] = cQ2P(seqQ[i]);
            }
        }
        
        if(jump_to[pos][i] > 0){
         pos = jump_to[pos][i];
        }else{
         pos--;   
        }
    }
        
        
        
    char * temp = strdup(seq);
                

    int count = 0;
    
        for(i = 0 ; i < dim; i++)
        {
         if(result_seq[i])
         {
          temp[count] = result_seq[i];
          count++;
         }
        }
        temp[count] = '\0';
        
        
        
        for(i = 0 ; i < dim; i++)
        {
            result_seq[i] = temp[i];
        }
        
        
        free(temp);
        
//                 print_selective("\n%s\n", result_seq);
}
    
    

int findEndReplicates(char *seq, unsigned int seq_length, double match, double miss, double jump){
    int i;
    int result;
    
        double** score = (double**)calloc(seq_length, sizeof(double*));
    for(i=0; i<seq_length; i++)
    {
    score[i] = (double*)calloc(seq_length, sizeof(double));
    }
    
    int** jump_to = (int**)calloc(seq_length, sizeof(int*));
    for(i=0; i<seq_length; i++)
    {
    jump_to[i] = (int*)calloc(seq_length, sizeof(int));
    }
    
    double* max = (double*)calloc(seq_length, sizeof(double));
    
    int* which_max = (int*)calloc(seq_length, sizeof(int));
    
    char* reverse_seq = strdup(seq);
    
    reverseSequence(reverse_seq, seq);
    calculateJumpMatrix(reverse_seq, seq_length, match, miss, jump, score, jump_to, max, which_max);
    result = findStartReplicates(seq, jump_to, seq_length, which_max[seq_length-1]);
    
//     TODO should there be a '-1' here?
    result = seq_length - result;    
    
    free(reverse_seq);
    for(i=0; i<seq_length; i++)
    {
        free(score[i]);
        free(jump_to[i]);
    }
    free(score);
    free(jump_to);
    free(max);
    free(which_max);
    
    return result;
}

// Finds the consensus of cirseq data reads. Returns the number of replicates found.
int alignReplicates(setting arg, char* result, double * result_phred, char *seq, char* seqQ, unsigned int seq_length){
    
    int i;
    
    double match = 3;
    double miss =  -1;
    
//     double jump = 0;
    double jump = -1* 4* 2 * match;
    
    
    

    
    double** score = (double**)calloc(seq_length, sizeof(double*));
    for(i=0; i<seq_length; i++)
    {
    score[i] = (double*)calloc(seq_length, sizeof(double));
    }
    
    int** jump_to = (int**)calloc(seq_length, sizeof(int*));
    for(i=0; i<seq_length; i++)
    {
    jump_to[i] = (int*)calloc(seq_length, sizeof(int));
    }
    
    double* max = (double*)calloc(seq_length, sizeof(double));
    
    int* which_max = (int*)calloc(seq_length, sizeof(int));
    
//     --------------------------------------
    calculateJumpMatrix(seq, seq_length, match, miss, jump, score, jump_to, max, which_max);
 
    
    
//     printReplicates(seq, jump_to, seq_length, which_max[seq_length-1]);
    
    int start = findStartReplicates(seq, jump_to, seq_length, which_max[seq_length-1]);
    int end   = findEndReplicates(seq, seq_length,  match,  miss,  jump);
    

    
    
     char* result_seq = (char*)calloc(seq_length+1, sizeof(char));
    findConsensusReplicates(result_seq, result_phred, seq, seqQ, jump_to, seq_length, which_max[seq_length-1], start, end);
    
      
//     print_selective("Start: %d End: %d\n", start, end);
    
    if(start<end)
    {
        for(i=0; i<seq_length; i++){
         result[i] = result_seq[i];   
        }
    }else{
        result[i] = seq[i];
    }
    
    int replication_number = findNumberOfReplicates(seq, jump_to, seq_length, which_max[seq_length-1], strlen(result));
    
//     TODO return this value
    free(result_seq);
    
    for(i=0; i<seq_length; i++)
    {
        free(score[i]);
        free(jump_to[i]);
    }
    free(score);
    free(jump_to);
    free(max);
    free(which_max);
    
    return replication_number;
}

int alignReplicates_old(char *seq, unsigned int seq_length){
    
    int lane,i,j;
    
    double match = 1;
    double miss =  -4;
    double add = miss;
    
    double* score = (double*)calloc(seq_length * seq_length, sizeof(double));
    
    int* jump_to = (int*)calloc(seq_length * seq_length, sizeof(int));
    
    double* max = (double*)calloc(seq_length, sizeof(double));
    
    int* which_max = (int*)calloc(seq_length, sizeof(int));
    
    double multiplier = 1.0;
    
    for(lane = 0; lane < seq_length ; lane++){
    
    multiplier = 1.0 + 2*lane/(double)seq_length;
    
        for(i = 0; i < seq_length - lane; i++){
            j = lane + i;
            
            if(seq[i] == seq[j])
                add = match;
            else
                add = miss;
            
            add = add*multiplier;
            
            
            if(i == 0 && j == 0)//fist element (no jump possible and no diagonal move)
            {
                score[i +j*seq_length] = add;
                max[j] = add;
                which_max[j] = i;
                
            }else if(i == 0 && j != 0){// top row (no diagonal parent)
                score[i +j*seq_length] = add;
                
                if(max[j-1]>0){
                    score[i +j*seq_length] += max[j-1];
                    jump_to[i +j*seq_length] = which_max[j-1];                    
                }
                
                if(max[j] < score[i +j*seq_length]){
                    max[j] = score[i +j*seq_length];
                    which_max[j] = i;
                }
            }else{// Regular recursion in the matrix
                
                
                if(score[(i-1) +(j-1)*seq_length] > max[j-1]){
                score[i +j*seq_length] = (score[(i-1) +(j-1)*seq_length] ) + add;
                }else{
                    score[i +j*seq_length] = max[j-1] + add;
                    jump_to[i +j*seq_length] = which_max[j-1];
                }
                
                
                if(max[j] < score[i +j*seq_length]){
                    max[j] = score[i +j*seq_length];
                    which_max[j] = i;
                }
            }
            
            
        }
    }
    
    
//     printCirseqMatrix(jump_to, seq, seq_length, which_max[seq_length-1]);
    
    printf("last row max: %f at: %i\n", max[seq_length-1], which_max[seq_length-1]);
 
    free(score);
    free(max);
 
    return 0;
}

int alignSingleDeletion(unsigned int length, char* ref1, char* ref2, char* read)
{
    
    assert(length > 0);
    
    int i;
    double temp, add;
    
    int begin_score1 = -1;
    int begin_score1_before_indel = 0;
    
    int indel_point = 0;
    double max_score2 = -1.0;
    int max_position_score2 = -1;
    
    double score1 = 0, score2 = 0;
    
    double match_score = 2.0;
    double mismatch_penalty = -4.0;
    
    int verbose = 1;
    
        if(mapUnsafe(ref1[0]) == mapUnsafe(read[0]))
            temp = score1 + match_score;
        else
            temp = score1 + mismatch_penalty;
            
        score1 = 0;
        if(temp > 0)
            score1 = temp;    
 
           
        if(mapUnsafe(ref2[0]) == mapUnsafe(read[0]))
            score2 = match_score;
        else    
            score2 = mismatch_penalty;
 
    
        for(i = 1; i < length; i++)
        {
//             Calculate value for tail part of the read 
            if(mapUnsafe(ref2[i]) == mapUnsafe(read[i]))
                add = match_score;
            else
                add = mismatch_penalty;
          
//             either insert indel here (from score1) or extend alignment
           if(score1 > score2)
           {
            score2 = score1 + add;
            begin_score1_before_indel = begin_score1;
            indel_point = i;
           }else{
            score2 = score2 + add;   
           }
           
           if(score2 > max_score2)
               max_position_score2 = i;
           
            
//             Calcualte score for head of squence.
//            The oder is important, as score2 must be calculated from Score1[i-1][i-1]
            if(mapUnsafe(ref1[i]) == mapUnsafe(read[i]))
                temp = score1 + match_score;
            else
                temp = score1 + mismatch_penalty;
            
            score1 = 0;
            if(temp > 0)
            {
                score1 = temp;
                if(begin_score1 < 0)
                    begin_score1 = i;
            }else{
             begin_score1 = -1;   
            }
            
            
        }
        
        
    if(begin_score1 > indel_point || begin_score1 < 0)
            begin_score1 = begin_score1_before_indel;
            
            
    if(verbose > 1)
        {
        print_selective("ref \tread\n");   
        
        for(i = 0; i < begin_score1; i++)
        {
         print_selective(" %c  X %c\n", ref1[i], read[i]);   
        }
         print_selective("---------\n");   
        
        for(i = begin_score1; i < indel_point; i++)
        {
            if(mapUnsafe(ref1[i]) == mapUnsafe(read[i]))
         print_selective(" %c --- %c\n", ref1[i], read[i]);  
         else
         print_selective(" %c  X  %c\n", ref1[i], read[i]);  
        }
                 print_selective("---------\n");   
        for(i = indel_point; i < length; i++)
        {
              if(mapUnsafe(ref1[i]) == mapUnsafe(read[i]))
        print_selective(" %c --- %c\n", ref1[i], read[i]);   
                  else
         print_selective(" %c  X  %c\n", ref1[i], read[i]);   
        }
        
         print_selective("\n- \t-\n\n");   
    
         for(i = 0; i < indel_point; i++)
         {
                  if(mapUnsafe(ref2[i]) == mapUnsafe(read[i]))
          print_selective(" %c --- %c\n", ref2[i], read[i]);   
          else
          print_selective(" %c  X  %c\n", ref2[i], read[i]);   
              
         }
         
         
         print_selective("---------\n");   
    
         
         for(i = indel_point; i < max_position_score2+1; i++)
        {
                        if(mapUnsafe(ref2[i]) == mapUnsafe(read[i]))
         print_selective(" %c --- %c\n", ref2[i], read[i]);   
         else
                      print_selective(" %c  X  %c\n", ref2[i], read[i]);   
        }
        
        
         print_selective("-     -\n");   
    
         
        for(i = max_position_score2+1; i < length; i++)
        {
         print_selective(" %c  X  %c\n", ref2[i], read[i]);   
        }
    }
        
   return 0; 
}





int alignSingleInsertion(unsigned int length, char* ref1, char* ref2, char* read, int offset)
{
    
    assert(length > 0);
    
    int i;
    double temp, add;
    
    int begin_score1 = -1;
    int begin_score1_before_indel = 0;
    
    int indel_point = 0;
    double max_score2 = -1.0;
    int max_position_score2 = -1;
    
    double score1 = 0, score2 = 0;
    
    double match_score = 2.0;
    double mismatch_penalty = -4.0;
    
    int verbose = 1;
    
        if(mapUnsafe(ref1[0]) == mapUnsafe(read[0]))
            temp = score1 + match_score;
        else
            temp = score1 + mismatch_penalty;
            
        score1 = 0;
        if(temp > 0)
        {
            begin_score1 = 0;
            score1 = temp;    
        }
           
        if(mapUnsafe(ref2[0 + offset]) == mapUnsafe(read[0 + offset]))
            score2 = match_score;
        else    
            score2 = mismatch_penalty;
 
    
        for(i = 1; i < length; i++)
        {
//             Calculate value for tail part of the read 
            if(mapUnsafe(ref2[i + offset]) == mapUnsafe(read[i + offset]))
                add = match_score;
            else
                add = mismatch_penalty;
          
//             either insert indel here (from score1) or extend alignment
           if(score1 > score2)
           {
            score2 = score1 + add;
            indel_point = i;
            begin_score1_before_indel = begin_score1;
           }else{
            score2 = score2 + add;   
           }
           
           if(score2 > max_score2)
               max_position_score2 = i + offset;
           
            
//             Calcualte score for head of squence.
//            The oder is important, as score2 must be calculated from Score1[i-1][i-1]
            if(mapUnsafe(ref1[i]) == mapUnsafe(read[i]))
                temp = score1 + match_score;
            else
                temp = score1 + mismatch_penalty;
            
            score1 = 0;
            if(temp > 0)
            {
                score1 = temp;
                if(begin_score1 < 0)
                    begin_score1 = i;
            }else{
             begin_score1 = -1;   
            }
            
            
        }
        
        
        if(begin_score1 > indel_point || begin_score1 < 0)
            begin_score1 = begin_score1_before_indel;
        
        
    if(verbose > 1)
        {
        print_selective("ref \tread\n");   
        
    
         print_selective("-Junk head- (going to %d)\n", (begin_score1 -1));   
        for(i = 0; i < begin_score1; i++)
        {
             if(mapUnsafe(ref1[i]) == mapUnsafe(read[i]))
         print_selective(" %c --- %c\n", ref1[i], read[i]);
         else
          print_selective(" %c  X  %c\n", ref1[i], read[i]);
           
        }
         print_selective("--------- Alignment head (starting at: %d)\n", begin_score1);   
        
        for(i = begin_score1; i < indel_point; i++)
        {
            if(mapUnsafe(ref1[i]) == mapUnsafe(read[i]))
         print_selective(" %c --- %c\n", ref1[i], read[i]);  
         else
         print_selective(" %c  X  %c\n", ref1[i], read[i]);  
        }
                 print_selective("--------- Unaligned head (starting at: %d)\n", indel_point);   
        for(i = indel_point; i < length; i++)
        {
              if(mapUnsafe(ref1[i]) == mapUnsafe(read[i]))
        print_selective(" %c --- %c\n", ref1[i], read[i]);   
                  else
         print_selective(" %c  X  %c\n", ref1[i], read[i]);   
        }
        
         print_selective("\n- unaligned tail -\n\n");   
    
         for(i = 0; i < indel_point + offset; i++)
         {
                  if(mapUnsafe(ref2[i]) == mapUnsafe(read[i]))
          print_selective(" %c --- %c\n", ref2[i], read[i]);   
          else
          print_selective(" %c  X  %c\n", ref2[i], read[i]);   
              
         }
         
         
         print_selective("--------- Matching region of tail (starting at: %d)\n", indel_point + offset);   
    
         
         for(i = indel_point + offset; i < max_position_score2+1; i++)
        {
                        if(mapUnsafe(ref2[i]) == mapUnsafe(read[i]))
         print_selective(" %c --- %c\n", ref2[i], read[i]);   
         else
                      print_selective(" %c  X  %c\n", ref2[i], read[i]);   
        }
        
        
         print_selective("-   Junk part of tail  -\n");   
    
         
        for(i = max_position_score2+1; i < length; i++)
        {
         print_selective(" %c  X  %c\n", ref2[i], read[i]);   
        }
    }
        
   return 0; 
}




char* alignAnyIndel(setting s, globalVariables* g, char* seq, unsigned int num);

// ----------------------------------------------------------------------------------------------------------------
// cloned from our gotoh submission

void print_matrix(int ** matrix, int n, int m)
{
  int i,j;

  for (i = 0; i <= n; ++i)
  {
    for (j = 0; j <= m; ++j)
    {
      printf ("%3d   ", matrix[i][j]);
    }
    printf ("\n");
  }
}

void print_dir(int ** matrix, int n, int m)
{
  int i,j,k;
  char bitmask[8];

  for (i = 0; i <= n; ++i)
  {
    for (j = 0; j <= m; ++j)
    {
      for (k = 0; k < 7; ++k) bitmask[k] = '0';
      bitmask[7] = 0;

      if (i > 0 && j > 0)
      {
        if (matrix[i][j] & DIAG) bitmask[0] = '1';
        if (matrix[i][j] & LEFTOPEN) bitmask[1] = '1';
        if (matrix[i][j] & UPOPEN) bitmask[2] = '1';
        if (matrix[i][j] & LEFT) bitmask[3] = '1';
        if (matrix[i][j] & UP) bitmask[4] = '1';
        if (matrix[i][j] & LEFTEXT) bitmask[5] = '1';
        if (matrix[i][j] & UPEXT) bitmask[6] = '1';
      }
      printf (" %s ", bitmask);
    }
    printf ("\n");
  }
  printf("\n BITMASK: DIAG LEFTOPEN UPOPEN LEFT UP LEFTEXT UPEXT\n");
}


int cstack_push(cstack_t ** stack, char c)
{
  cstack_t * new = (cstack_t *)malloc(sizeof(cstack_t));
  if (!new) return 0;

  new->c = c;
  new->next = *stack;
  *stack = new;

  return 1;
}

char cstack_pop(cstack_t ** stack)
{
  char c;
  cstack_t * tmp;

  if (!*stack) return 0; 

  tmp = (*stack);
  c = (*stack)->c;
  *stack = (*stack)->next;
  free(tmp);

  return c;
}

void cstack_clear(cstack_t ** stack)
{
  while (*stack) cstack_pop(stack);
}

int cstack_size(cstack_t ** stack)
{
  int size;
  cstack_t * top = *stack;

  for (size = 0; top; top = top->next, ++size);

  return size;
}

void cstack_print(const char * prefix, cstack_t ** stack)
{
  cstack_t * top = *stack;
  printf("%s: ", prefix);

  for (;top; top = top->next)
  {
    printf("%c", top->c);
  }
  printf("\n\n");
  
}


void backtrackAffine(int ** d, int ** e, const char * a, const char * b, int i, int j, cstack_t ** astack, cstack_t ** bstack)
{
  int k,m;

  if (i > 0 && j > 0)
  {
    if (e[i][j] & LEFT)
    {
       #ifdef DEBUG
       printf("LEFT\n");
       #endif
       k = j;
       /* process the gap extensions */
       while (e[i][k] & LEFTEXT)
       {
         #ifdef DEBUG
         printf("LEFTEXT\n");
         #endif
         cstack_push(astack, a[k-1]);
         cstack_push(bstack, '-');
         --k;
       }
       /* process the gap open */
//       if (k)
//       {
         cstack_push(astack, a[k-1]);
         cstack_push(bstack, '-');
         --k;
//       }
       backtrackAffine(d,e,a,b,i,k,astack,bstack);
       
       /* bring the stack to its original state */
       for (m = 0; m < j-k; ++m)
       {
         cstack_pop(astack);
         cstack_pop(bstack);
       }
    }
    if (e[i][j] & UP)
    {
       #ifdef DEBUG
       printf("UP\n");
       #endif
       k = i;
       /* process the gap extensions */
       while (e[k][j] & UPEXT)
       {
         #ifdef DEBUG
         printf("UPEXT\n");
         #endif
         cstack_push(astack, '-');
         cstack_push(bstack, b[k-1]);
         --k;
       }
       /* process the gap open */
 //      if (k)
 //      {
         cstack_push(astack, '-');
         cstack_push(bstack, b[k-1]);
         --k;
 //      }
       backtrackAffine(d,e,a,b,k,j,astack,bstack);

       /* bring the stack to its original state */
       for (m = 0; m < i-k; ++m)
       {
         cstack_pop(astack);
         cstack_pop(bstack);
       }
    }
    if (e[i][j] & LEFTOPEN)
    {
      #ifdef DEBUG
      printf("LEFTOPEN\n");
      #endif
      /* single indel gap */
      cstack_push(astack, a[j-1]);
      cstack_push(bstack, '-');
      backtrackAffine(d,e,a,b,i,j-1,astack,bstack);
      cstack_pop(astack);
      cstack_pop(bstack);
    }
    if (e[i][j] & UPOPEN)
    {
      #ifdef DEBUG
      printf("UPOPEN\n");
      #endif
      /* single indel gap */
      cstack_push(astack, '-');
      cstack_push(bstack, b[i-1]);
      backtrackAffine(d,e,a,b,i-1,j,astack,bstack);
      cstack_pop(astack);
      cstack_pop(bstack);
    }
    if (e[i][j] & DIAG)
    {
      #ifdef DEBUG
      printf("DIAG\n");
      #endif
      /* match/mismatch */
      cstack_push(astack, a[j-1]);
      cstack_push(bstack, b[i-1]);
      backtrackAffine(d,e,a,b,i-1,j-1,astack,bstack);
      cstack_pop(astack);
      cstack_pop(bstack);
    }
  }
  else
  {
    /* number of push operations we will need to get to cell (0,0) */
    int pop_count = i+j;

    /* get to (0,j) */
    while (i > 0)
    {
      #ifdef DEBUG
      printf("PUSHUP\n");
      #endif
      cstack_push(astack, '-');
      cstack_push(bstack, b[i-1]);
      --i;
    }

    /* get to (0,0) */
    while (j > 0)
    {
      #ifdef DEBUG
      printf("PUSHLEFT\n");
      #endif
      cstack_push(astack, a[j-1]);
      cstack_push(bstack, '-');
      --j;
    }

    /* print sequences */
    cstack_print("A", astack);
    cstack_print("B", bstack);
    printf("\n");

    /* reset stack to its original state */
    for (i = 0; i < pop_count; ++i)
    {
      cstack_pop(astack);
      cstack_pop(bstack);
    }
  }
}


int alignAffine(unsigned int length_a, const char * a, unsigned int length_b, const char * b)
{
  int i,j;
  
//   Originally formulated as minimization problem, so we multiply by -1
  int mismatch = (-1)*(-4);
  int open = (-1)*(-10);
  int ext = (-1)*(-1);
  
  int alen,blen;
  int match_score;
  int ** d;
  int ** p;
  int ** q;
  int ** e;
  
//   alen = strlen(a);
  alen = length_a;
//   blen = strlen(b);
 blen = length_b;

  /* allocate three matrices for the DP algorithm and another matrix
     e to hold the cell pointers */
  d = (int **)malloc((blen+1)*sizeof(int *));
  p = (int **)malloc((blen+1)*sizeof(int *));
  q = (int **)malloc((blen+1)*sizeof(int *));
  e = (int **)malloc((blen+1)*sizeof(int *));

  for (i = 0; i <= blen; ++i)
  {
    d[i] = (int *)calloc(alen+1, sizeof(int));
    p[i] = (int *)calloc(alen+1, sizeof(int));
    q[i] = (int *)calloc(alen+1, sizeof(int));
    e[i] = (int *)calloc(alen+1, sizeof(int));
  }
  
  /* initialize first row and column of D */
  d[0][0] = 0;
  for (i = 1; i <= alen; ++i)
    d[0][i] = open+i*ext;
  for (i = 1; i <= blen; ++i)
    d[i][0] = open+i*ext;

  /* initialize first row of P and first column of Q */
  for (i = 1; i <= alen; ++i)
    p[0][i] = 2*open + i*ext;
  for (i = 1; i <= blen; ++i)
    q[i][0] = 2*open + i*ext;

  /* fill dynamic programming matrix */
  for (i = 1; i <= blen; ++i)
  {
    for (j = 1; j <= alen; ++j)
    {
      p[i][j] = MIN(d[i-1][j] + open + ext, p[i-1][j] + ext);
      q[i][j] = MIN(d[i][j-1] + open + ext, q[i][j-1] + ext);

      match_score = (a[j-1] == b[i-1]) ? 0 : mismatch;
      d[i][j] = MIN(d[i-1][j-1] + match_score, MIN(p[i][j],q[i][j]));

      if (d[i][j] == d[i-1][j-1] + match_score) e[i][j] = e[i][j] | DIAG;
      if (d[i][j] == d[i-1][j] + open + ext) e[i][j] = e[i][j] | UPOPEN;
      if (d[i][j] == d[i][j-1] + open + ext) e[i][j] = e[i][j] | LEFTOPEN;
      if (d[i][j] == p[i-1][j] + ext && i > 1) e[i][j] = e[i][j] | UP;
      if (d[i][j] == q[i][j-1] + ext && j > 1) e[i][j] = e[i][j] | LEFT;

      if (p[i][j] == p[i-1][j]+ext && i > 1) e[i][j] = e[i][j] | UPEXT;
      if (q[i][j] == q[i][j-1]+ext && j > 1) e[i][j] = e[i][j] | LEFTEXT;
    }
  }

  /* use two stacks for storing the aligned sequences when backtracking */
  cstack_t * astack = NULL;
  cstack_t * bstack = NULL;

  /* recursively backtrack and print all optimal alignments */
  backtrackAffine(d,e,a,b,blen,alen,&astack,&bstack);

  printf("D:\n");
  print_matrix(d,blen,alen);
  printf("\n\n");
  printf("P:\n");
  print_matrix(p,blen,alen);
  printf("\n\n");
  printf("Q:\n");
  print_matrix(q,blen,alen);
  printf("\n\n");

  printf("E:\n");
  print_dir(e,blen,alen);
  printf("\n\n");

  /* deallocate memory */
  for (i = 0; i <= blen; ++i)
  {
    free(d[i]);
    free(p[i]);
    free(q[i]);
    free(e[i]);
  }
  free(d); free(p); free(q); free(e);

  return 0;
}
