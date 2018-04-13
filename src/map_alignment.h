
#ifndef _MAP_ALIGN   /* Include guard */
#define _MAP_ALIGN

#define DIAG         1
#define UPOPEN       2
#define LEFTOPEN     4
#define UP           8
#define LEFT        16
#define UPEXT       32
#define LEFTEXT     64

typedef struct cstack_s
{
  char c;
  char q; //optional quality as phred
  struct cstack_s * next;
} cstack_t;

char combine_exact(double * result_p_error, char c, double p_error_c, char s, double p_error_s);

int alignReplicates(setting arg, resultsVector *rv, char* result, char * result_phred, char *seq, char* seqQ, unsigned int seq_length);

int alignSingleDeletion(unsigned int length, char* ref1, char* ref2, char* read);

int alignSingleInsertion(unsigned int length, char* ref1, char* ref2, char* read, int offset);

int alignCirceq_call(setting s, globalVariables* g, int* positive_placements, int num_placements, char * seq, char* seqQ);

int alignAffine(unsigned int length_a, const char * a, unsigned int length_b, const char * b);
#endif