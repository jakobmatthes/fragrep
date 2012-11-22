#ifndef ALN_UTIL_H
#define ALN_UTIL_H

extern "C"
{

	extern int read_alignment(FILE *faln, char *AlignedSeqs[], char *names[]);

	extern int read_clustal(FILE *clust, char *AlignedSeqs[], char *names[]);

	extern int read_stockholm(FILE *clust, char *AlignedSeqs[], char *names[]);

	extern /*@only@*/ /*@notnull@*/ char *consensus(const char *AS[]);

	extern /*@only@*/ /*@notnull@*/ char *consens_mis(const char *AS[]);

}

#endif
