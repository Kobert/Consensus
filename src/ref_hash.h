#ifndef _HASH   /* Include guard */
#define _HASH



void freeTable(hashEntry* entryList, unsigned int length);

//void addToHash(unsigned int *slidingWindow, unsigned int sequencePosition, hashEntry **hashTable, hashEntry *entryTable, unsigned int *itemsInTable, unsigned int *hit, unsigned int *chained, unsigned int hashValue, unsigned int hashTableSize);
int testReference(globalVariables* g, setting s);

int populateHashTable(setting * arg, globalVariables* globalVar, char *seq, hashEntry ***hashTable, hashEntry **entryTable, unsigned int *itemsInTable);

int populateHashTableWithKey(setting * arg, globalVariables* globalVar, char *seq,  hashEntry ***hashTable, hashEntry **entryTable, unsigned int *itemsInTable, unsigned int hashValue, unsigned int hashTableSize);

int getNthLine(int n, FILE * file, char ** read, size_t * bytes);

int getNextNonEmptyLine(FILE * file, char ** read, size_t * bytes);

int getNextRead(FILE * file, char ** read, unsigned int * bytes, unsigned int * firstPosition);

int getNextQuality(FILE * file, char ** read, unsigned int * bytes, unsigned int * firstPosition);

void hashMapReads(setting arg, globalVariables *globalVar, resultsVector *rv);

int placeFragments(setting s, globalVariables* g, resultsVector* rv, hashEntry ** hashTable, hashEntry * entryTable, char * seq, unsigned int firstPosition, unsigned int * hit, unsigned int *miss, unsigned int * hitPerRound, unsigned int *max_fragments);
#endif
