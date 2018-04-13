
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>


#include "referenceAssembly.h" 
#include "ref_arg.h"			// Functions to interpret user input and initialize datastructures go here.
#include "ref_primes.h"			// Functions to test combinations of primenumbers for efficient hashing techniques
#include "ref_slidingWindow.h"		// Contains all functions that Manipulate or access the slidingWindow, or actual bases themselves.
#include "ref_hash.h"			// All functions that manipulate or access the tashtable or entries go here.
#include "ref_io.h"			// Functions that hanlde input, output, as well as result formatting (should) go here.
#include "ref_math.h"			// Anything more complicated than simple multiplication goes here

#include "map_arg.h"
#include "map_alignment.h"

#include "cir_iterate.h"

int main(int argc, char **argv)
{
    
    
//     double temp_error;
// char c = combine_exact(&temp_error, 'A', 0.5, 'A', 0.001);

// print_selective("c: %c, phred: %d (e: %f)\n", c, P2Q(temp_error), temp_error);

// exit(0);
//     Bad Case Zika
//     char string[] = "GGGGGTGGGGTGTTCAGGGCAGAACAGCAAGGCGAAAAGAAGCTTTATGAGAAACTGTTTTGAGTTCTGTTAATAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAACAAAACCTTCGTCTTCTCGCATACCTGAGCTTTTATATGTACGTGGGTGGGNNNN\0";

//     Good Case Zika
//     char string[] = "AGCGAAGCCACGCTGCTGAACATGCTCAACATCTCCCCCTTCTCCTTTGGGCTGGTCAAGACTGGAGACAAAGTGGGAGCCAGCGAAGCCACGCTGCTGAACATGCTCAACATCCCCCCCTTCTCCTTTGGGCTGGTCAAGACTGGAGACAAAGTGGGCGCCAGCGAAGCCACGCTGCTGAACCTGCTCAAACTCTCCCCCTTCTCCTTNTGNCNNGNNNNGNCNNGNNANNNNNNNNNNNNNNNNNNNNNNNNGCCACTCGAANNGCTCNACNTNCNACCACACTCCCTTCGATCTGTCAACACTTGAGCTACACTCCGCAC\0";
// same in reverse
//     char string[] = "CACGCCTCACATCGAGTTCACAACTGTCTAGCTTCCCTCACACCANCNTNCANCTCGNNAAGCTCACCGNNNNNNNNNNNNNNNNNNNNNNNNANNGNNCNGNNNNGNNCNGTNTTCCTCTTCCCCCTCTCAAACTCGTCCAAGTCGTCGCACCGAAGCGACCGCGGGTGAAACAGAGGTCAGAACTGGTCGGGTTTCCTCTTCCCCCCCTACAACTCGTACAAGTCGTCGCACCGAAGCGACCGAGGGTGAAACAGAGGTCAGAACTGGTCGGGTTTCCTCTTCCCCCTCTACAACTCGTACAAGTCGTCGCACCGAAGCGA\0";
    
    //     Good Case Polio
//     char string[] = "NAAGTTTTGTTACAAATTCCAACTTATCTCTAGCTTGTGGGATAATTTTCTCCTTGAGCCAAGATTTGGTTTTCCAGCATTTCTAGTTGTCTAAGTTTTGTTACAAATTCCAACTTATCTCTAGCTTGTGGGATAATTTTCTCCTTGAGCCAAGATTTGGTTTTCCAGCATTTCTAGTTGTCTATGTTTTGTTACAAATTCCAACTTATNTCNANNTNNNNNNNNNANNTNNNNNNNNNNNNNNNNNNNNNNNNTCNNGNANTTNNTNTNNNNNNNNNTTTGTTACAATGTCCAACTTATCTCTAGCTTGTCGGATACTTCTC\0";
//     char string[] = "NTTGAGTTAAAAATTGAAGTGCCTGAGCAGCCAGATGGCATACCGCCCTTGACACAGTATGTTTTATCCTGATAATCAAGTTGTTAATCATTGTGTTAAAAATTGAAGTGCCTGAGCAGCCAGATGGCATACCGCCCTTGACACAGTATGTTTTATCCTGANAATCAAGATGTTAACCATTGAGTTAAAAATTGAAGTGCCTGAGCAGCNNNNNNNNNNNNNNNNNTNNNNNNNNNNNNNNNNNNNNNNNNNNNCANNTNGNNANNNNCNNNNNNNNNATTTGAACTGACTCACCATCTAGTGGACGTCAAGAACTGGACCAG\0";
//       char string[] = "ACGGGTACGGGTACGGGT\0";
//   alignReplicates(string,strlen(string));
//   exit(0);
  

  
 setting arg = cirseq_initArgs(argc, argv);
//  setting arg = map_initArgs(argc, argv);
 
 globalVariables globalVar;
  
 initGlobalVariables(&arg, &globalVar);
  
 resultsVector rv;
 
 initResults(&rv, strlen(globalVar.referenceSequence));
 
  fflush(stdout);
  

  
  
  if(arg.executeReferenceOnly)
  {
      if(arg.num_multi_references == 1)
     testReference(&globalVar, arg);
          
     print_selective("[NOTE:]    \"-x\" set. Exiting now...\n"); 
      
     cleanForExit(&arg, &globalVar,  &rv);
  
  return 0;   
  }
  
  

cirseqReads(arg, &globalVar, &rv);
//   hashMapReads(arg, &globalVar, &rv);


//   print_selective("\nGood hits: %d, multi hits %d bad hits %d (bad + multi %.11f)\n", good_hit, multi_hit, bad_hit, (bad_hit+multi_hit)/(double)(good_hit+bad_hit+multi_hit));
  
  postProcessResults(arg, &rv);

//   printAvgCoverage(arg, rv);
  
 

 
 
 
//  if(arg.doCsvFile)
// {
//     char *csvFileName;
//   if(arg.outFilePrefix)
// {
//   csvFileName = (char*)calloc(strlen(arg.outFilePrefix)+5, sizeof(char));
// sprintf(csvFileName, "%s.csv", arg.outFilePrefix);
// }else{
//   csvFileName = (char*)calloc(10, sizeof(char));
// sprintf(csvFileName, "_temp.csv");
// }
// FILE * csvFile = fopen(csvFileName,"w");
// 
// printCSV(csvFile, rv);
// fclose(csvFile);
// print_selective("\ncsv output written to:           \"%s\"\n", csvFileName);
// free(csvFileName);
// }



// char *htmlFileName;
//     if(arg.outFilePrefix)
//     {
//     htmlFileName = (char*)calloc(strlen(arg.outFilePrefix)+6, sizeof(char));
//     sprintf(htmlFileName, "%s.html", arg.outFilePrefix);
//     }else{
//     htmlFileName = (char*)calloc(11, sizeof(char));
//     sprintf(htmlFileName, "_temp.html");
//     }
// FILE * htmlFile = fopen(htmlFileName,"w");
// 
// printHtml(htmlFile, arg, rv);
// fclose(htmlFile);
// print_selective("\nHtml output written to:           \"%s\"\n", htmlFileName);
// free(htmlFileName);




if(arg.doTrim)
{
print_selective("\n\t Avg trimmed ratio:   %.4f (Average per read)\n",globalVar.avgRatioTrimmed);

print_selective("\t Total trimmed ratio: %.4f (Fraction of sites)\n\n",(double)globalVar.trimmed/(globalVar.trimmed + globalVar.kept));
} 
  
  
  
  print_selective("\n\t Number of cirseq reads identified: %u\n", rv.cir_num );  
  print_selective("\n\t Average length of cirseq consensus sequences: %.2f\n", (rv.cir_total_length/(double)rv.cir_num) );
  print_selective("\t Average number of replicates per cirseq read: %.2f\n", (rv.cir_foldings/(double)rv.cir_num));

  print_selective("\nAverage error rate for cirseq reads:  %f\n", rv.mean_errors);
//     print_selective("Average Phred score for cirseq reads: %f\n", rv.mean_phred);
//   if(arg.referenceFile){
//    print_selective("\n\t Resulting average coverage on reference: %.2f\n", rv.on_ref_total_length/(double)globalVar.referenceSequenceLength);   
//   }
  
  print_selective("\n");
  
  printTime();
  print_selective(" Freeing data-structures...\n");
  cleanForExit(&arg, &globalVar,  &rv);
  
  printTime();
  print_selective(" Exiting.\n");
  
  return 0;

}



