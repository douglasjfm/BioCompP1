#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "fasta.h"

FASTA* abrirFasta(char *seqfile)
{
    FASTA *ffp;

    ffp = malloc(sizeof(FASTA));
    ffp->fp = fopen(seqfile, "r");
    if (ffp->fp == NULL)
    {
        free(ffp);
        return NULL;
    }
    if ((fgets(ffp->buffer, FASTA_MAXLEN, ffp->fp)) == NULL)
    {
        free(ffp);
        return NULL;
    }
    return ffp;
}

int lerFasta (FASTA *ffp, char **ret_seq, char **ret_name, unsigned *ret_L)
{
    char *s;
    char *name;
    char *seq;
    int   n;
    int   nalloc;


    if (ffp->buffer[0] != '>') return 0;

    s  = strtok(ffp->buffer+1, " \t\n");
    name = malloc(sizeof(char) * (strlen(s)+1));
    strcpy(name, s);

    seq = malloc(sizeof(char) * 128);
    nalloc = 128;
    n = 0;
    while (fgets(ffp->buffer, FASTA_MAXLEN, ffp->fp))
    {
        if (ffp->buffer[0] == '>') break;

        for (s = ffp->buffer; *s != '\0'; s++)
        {
            if (! isalpha(*s)) continue;

            seq[n] = *s;
            n++;
            if (nalloc == n)
            {
                nalloc += 128;
                seq = realloc(seq, sizeof(char) * nalloc);
            }
        }
    }
    seq[n] = '\0';

    if (ret_name) *ret_name = name;
    *ret_seq  = seq;
    *ret_L    = n;
    return 1;
}

void encerraFasta(FASTA *ffp)
{
    fclose(ffp->fp);
    free(ffp);
}
