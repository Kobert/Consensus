//All functions that manipulate or access the hashtable or entries go here.
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "referenceAssembly.h"
#include "ref_slidingWindow.h"
#include "ref_hash.h"
#include "ref_io.h"
#include "ref_math.h"
#include "ref_arg.h"
#include "map_alignment.h"





void cirseqReadsWithKey(setting arg, globalVariables *globalVar, resultsVector *rv, unsigned int hashValue, unsigned int hashTableSize);

void cirseqReads(setting arg, globalVariables *globalVar, resultsVector *rv)
{
   cirseqReadsWithKey( arg, globalVar, rv, _hashValue, _hashTableSize);
}


void duplicateSequence(char *duplicate, char* seq){
       
    int i;
    int len = strlen(seq);
    for(i=0; i<len; i++)
    {
     duplicate[i] = seq[i];   
     duplicate[i+len] = seq[i];   
    }
        
}

void cirseqReadsWithKey(setting arg, globalVariables *globalVar, resultsVector *rv, unsigned int hashValue, unsigned int hashTableSize)
{
  
//     print_selective("csv:\ti\tcount\tlength\tconsens\tname\n");
    
    int i;

  //Open File and check that it worked
      FILE *readsFile;
  readsFile = fopen(arg.readsFile,"r");
  if(!readsFile){
  readsFile = fopen64(arg.readsFile,"r");//For 32 bit systems this may solve the issue
  }
  if(!readsFile){
    fprintf(stderr, "\n[ERROR:] Could not open file: \"%s\"\n", arg.readsFile);
      fprintf(stderr, "         Please make sure the path is correct and you have read permission.\n");
      assert(readsFile);
  }
  

//   FILE *samFile;
//   samFile = fopen("_temp_placeholder.fq","w");
  
  //Allocate space for short reads and complement thereof
 // unsigned int bytesToRead = 10;
  size_t bytesToRead = 10;
  char * seq;
  seq = (char*)malloc(bytesToRead );

  size_t bytesToRead_name = 10;
  char * seq_name;
  seq_name = (char*)malloc(bytesToRead_name );
  
//unsigned int bytesToReadQ = 10;
  size_t bytesToReadQ = 10;
  char * seqQ;
  seqQ = (char*)malloc(bytesToReadQ );
  
  unsigned int 	numReads = 0,
		numMatches = 0,
		numForwardMatches = 0,
		numReverseMatches = 0,
		numShiftedMatches = 0;
  
  
  unsigned int 	hit = 0,
		miss = 0,
		hitPerRound = 0;
  
  // wichPosition stores the position of the reference genome, where the left most base of the short read is aligned.
  int whichPos = -1;
  
  // Variables to keep track of whether a matching pattern has been found, and what configuration was used. 
  unsigned int 	matchForward = 0,
		matchReverse = 0,
		matchShifted = 0;
  unsigned int match = 0;
  
  unsigned int max_fragments = 0;
//   unsigned int mappingQuality = 255;
  
//     unsigned int readDiff = 2;//During the first iteration of the while loop only read the name and then the read
  unsigned int readDiff = 1;//During the first iteration of the while loop only read the name and then the read

//   int cirSeq_reads_found = 0;
  // This loop goes over all short read's dna data
	while(getNthLine(readDiff, readsFile, &seq_name, &bytesToRead_name))
	{
//              printTime();
//              print_selective(" in while: Name %s\n", seq_name);
//              fflush(stderr);
	  getNthLine(readDiff, readsFile, &seq, &bytesToRead);
          
          readDiff = 2;
	  getNthLine(readDiff, readsFile, &seqQ, &bytesToReadQ);
	        readDiff = 1;
                
          char * complement;
	  complement = strdup(seq);
	  
	  char * result;

          char * consensus = strdup(seq);
          char * consensus_phred = (char*)malloc( (strlen(seq)+1)*sizeof(char) );
	  unsigned int shift = 0;

	  
// 	readDiff = 4;	// Reads are expected to fall 4 non-empty lines apart. 
			// (in between are "+", quality profile, and name of next read starting with "@")
// 	readDiff = 3;	// Reads are expected to fall 4 non-empty lines apart, but we want the name next. 
	numReads ++;
	
	matchForward = 0;
	matchReverse = 0;
	matchShifted = 0;
	match = 0;
        max_fragments = 0;
//         mappingQuality = 255;
	
        
        
//         First test whether the read maps to the reference
//          if(arg.referenceFile  && 0<floor( (double)strlen(consensus)/(double)basesPerWindow()) ){
//              
//              	whichPos = placeFragments(arg, globalVar, rv, globalVar->hashTable, globalVar->entryTable, seq, 0, &hit, &miss, &hitPerRound, &max_fragments);
//                 
//                 
//          }
        
        
        
        
        
        
        
        int repeat_number = alignReplicates(arg, rv, consensus, consensus_phred, seq, seqQ, strlen(seq));
        char* original_consensus = strdup(consensus);
        
        if(repeat_number > 1){
            rv->cir_num++;
            rv->cir_total_length += strlen(consensus);
            rv->cir_foldings += repeat_number;
        }
        
        
//         If a reference is provided and we belive there to be a cirSeq read, try to align it to the reference
if(0)//Unused code. Needs new implementation TODO
 if(arg.referenceFile && (repeat_number > 1) && 0<floor( (double)strlen(consensus)/(double)basesPerWindow())   )
 {
     
     
     char* duplicate = (char*)calloc(2*strlen(seq)+1, sizeof(char));
     duplicateSequence(duplicate, consensus);
          
     
//      flip around some pointers to avoid renaming all variables....
     char * temp;
     temp = seq;
     seq = duplicate;
     duplicate = temp;
     
     
     
     free(complement);
     complement = strdup(seq);

     //      printf("seq bf: %s\n",seq);
	// Find the left most position where the read aligns to the reference.
	// The first try is performed on the regular sequence with no offset.
	whichPos = placeFragments(arg, globalVar, rv, globalVar->hashTable, globalVar->entryTable, seq, 0, &hit, &miss, &hitPerRound, &max_fragments);
	
	  if( whichPos >= 0 && whichPos <= (globalVar->referenceSequenceLength - strlen(seq)) )
	  {
	    result = seq;
	      matchForward = 1;
	      match = 1;
	  }
		
	  // Next, test whether the reverse complemented sequence matched the reference.
	  if(!match)
	  {
	   reverseComplementSequence(complement, seq);
    
	   if(strlen(complement) != strlen(seq))
	   {
	   fprintf(stderr, "\nstrlen(complement) = %lu, strlen(seq) =%lu\n", (long unsigned)strlen(complement), strlen(seq));
	   assert(strlen(complement) == strlen(seq));
	   }
	   
	   whichPos = placeFragments(arg, globalVar, rv, globalVar->hashTable, globalVar->entryTable, complement, 0, &hit, &miss, &hitPerRound, &max_fragments);
	  
	
	  if( whichPos >= 0 && whichPos <= (globalVar->referenceSequenceLength - strlen(seq)) )
	  {
	    result = complement;
	      matchReverse = 1;
	      match = 1;
	  }		
	  }

	  // Now, test whether shifting it such that the fragments are right aligned gives a valid matching.
	  // If it was already right aligned during the first iteration, use a centered method instead.
	  if(!match  && arg.doExtendedTesting)
	  {

	    shift = strlen(seq) % basesPerWindow();
	    if(shift == 0)
	    {
	    //  shift = floor(basesPerWindow()/2.0);
	    }
	    
	   whichPos = placeFragments(arg, globalVar, rv, globalVar->hashTable, globalVar->entryTable, seq, shift, &hit, &miss, &hitPerRound, &max_fragments);
	   whichPos = whichPos - shift;
	   
	  if( whichPos >= 0 && whichPos <= (globalVar->referenceSequenceLength - strlen(seq)) )
	  {
	    result = seq;	    
	      matchForward = 1;
	      match = 1;
	  }		
	  }	

	  // Lastly, do the shifting for the reverse complemented sequence
	  if(!match && arg.doExtendedTesting)
	  {

	    
	   whichPos = placeFragments(arg, globalVar, rv, globalVar->hashTable, globalVar->entryTable, complement, shift, &hit, &miss, &hitPerRound, &max_fragments);
	   whichPos = whichPos - shift;
	   
	  if( whichPos >= 0 && whichPos <= (globalVar->referenceSequenceLength - strlen(seq)) )
	  {
	    result = complement;	    
	      matchReverse = 1;
	      match = 1;
	  }		
	  }	
	  
	  
	  // Finally, when all relevant configurations have been testes, or one good configuration found, evaluate the alignment
	  if(match)
	  {
              
//             mappingQuality = mapping_quality(arg, max_fragments, globalVar->referenceSequenceLength, strlen(seq));

	    numMatches++;
	    rv->on_ref_num++; //For Cirseq output

            
            
	    if(matchForward)
	    {
	      numForwardMatches++;
	    }
	    if(matchReverse)
	    {
	      numReverseMatches++; 
	    }
	    if(matchShifted)
	    {
	      numShiftedMatches++;
	    }
	    
	    
	    if(arg.doTrim)
		  {
//                       printf("trimming\n");
		  unsigned int an,bn;
		  trimRegions(arg, globalVar, result, whichPos, &an, &bn);
// 		  readToResult(rv, &(result[an]), whichPos+an, bn -an); 
                    

                  for(i=0; i< (bn-an) ; i++)
                  {
                      consensus[i]=result[an+i];
                  }
                  consensus[i]='\0';
                  
		    if(matchReverse)
		    {
		    reverseSequence(seqQ, seqQ);  
		    }
  


		  
		  }else{
// 		readToResult(rv, result, whichPos, strlen(seq)); 
	    
// 	        readDiff = 2;
// 		getNthLine(readDiff, readsFile, &seqQ, &bytesToReadQ);
		    if(matchReverse)
		    {
		    reverseSequence(seqQ, seqQ);  
		    }

		  
		  }
		   //For Cirseq output
	    rv->on_ref_total_length += strlen(consensus);
		  
		  
	  }else{ 

              

	    
	  }
// 	       printf("seq af: %s\n",seq);
	  
               temp = duplicate;
               duplicate = seq;
               seq = temp;
               
               
               
               free(duplicate);
}

// TODO do actual quality scores
char* temp_qual = strdup(consensus_phred);

// for(i=0; i < strlen(consensus); i++)
// {
//  temp_qual[i] = 'I';   
// }
// temp_qual[i] = '\0';

// printf("%s\n", seq_name);
// printf("%s\n", consensus);
// printf("+\n");
// printf("%s\n", temp_qual);


if(repeat_number > 1){
    
/*    print_selective("csv:\t%d\t%d\t%d\t%s\t%s", cirSeq_reads_found, repeat_number, strlen(consensus), consensus, seq_name);
    if(repeat_number*(strlen(consensus)*0.8)+strlen(consensus)*0.2 > strlen(seq))
    {
     print_selective("\tBAD_Length");   
    }else{
    print_selective("\tGOOD_Length");   
    }
    print_selective("\n");  */ 
    
//     cirSeq_reads_found++;
//     exit(0);
printf("%s\n", seq_name);
printf("%s\n", consensus);
printf("+\n");
printf("%s\n", temp_qual);

assert(strlen(consensus) == strlen(temp_qual));

// printf("\n%s\n", seq);
// printf("%s\n", original_consensus);
// printf("num matches: %u\n\n", numMatches);
}

	  readDiff = 1; // We should be at quality and want to advance to the name.
	  free(temp_qual);
	  free(complement);
          free(consensus);
          free(consensus_phred);
          free(original_consensus);
	}
  

       
  print_selective("\n\tNumber of reads: %u, number of reads mapped to the reference: %u (%.2f%%)\n", numReads, numMatches, 100.0*(double)numMatches/numReads);
  print_selective("\t(Forward %.2f%%, reverse %.2f%% (Shifted: %.2f%%))\n", 100.0*(double)numForwardMatches/numReads,
								 100.0*(double)numReverseMatches/numReads, 
								 100.0*(double)numShiftedMatches/numReads);
  
  
//   print_selective("Indels: %u\n", rv->indels);
  //Free allocated data structures
//   fclose(samFile);
  
  fclose(readsFile);
  
  free(seq_name);
  free(seqQ);
  free(seq);
  
}





