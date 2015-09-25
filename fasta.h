#include <stdio.h>

#define FASTA_MAXLEN 512	/* Limite de caracteres em uma linha do FASTA */

typedef struct fastafile_s {
  FILE *fp;
  char  buffer[FASTA_MAXLEN];
} FASTA;

extern FASTA *abrirFasta(char *seqfile);
extern int    lerFasta(FASTA *fp, char **ret_seq, char **ret_name, unsigned *ret_L);
extern void   encerraFasta(FASTA *ffp);
