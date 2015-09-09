#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*!
    @author douglas medeiros <djfm at cin.ufpe.br>

    1) Instrucoes:

        comando:
            ./P1 <cadeia 1> <cadeia 2>
        saída:
            P1 Console <\Nova Linha>
            Cadeias OK! <\Nova Linha>
            <\Nova Linha>
            --cadeia 1 alinhada--<\Nova Linha>
            --cadeia 2 alinhada--<\Nova Linha>
            score: <Inteiro>
            <Fim>

    2) Exemplo de uso:
        /caminho/do/executavel$ ./P1 acgat gat
        P1 Console
        Cadeias OK!

        acgat
        --gat
        score: -24
*/

#define GAP '-'

/*Pontuações*/
#define SCR_M  5
#define SCR_S -1
#define SCR_R -6
#define SCR_I -8

/*Dierções*/
#define DIRSUB 0x01
#define DIRESQ 0x02
#define DIRDIA 0x04

#define goUP(x) (x & DIRSUB)
#define goLF(x) (x & DIRESQ)
#define goDG(x) (x & DIRDIA)

/*max3: o valor maximo de uma tripla (x, y, z)*/
#define max3(x,y,z) ((x > y) && (x > z)) ? x : ((y > z) ? y : z)

int** allocintmtx(int i, int j)
{
    int **ret;
    ret = (int **) calloc (i, sizeof(int*));
    int c;
    for (c = 0; c < i; c++)
        ret[c] = (int *) calloc(j, sizeof(int));

    return ret;
}

void freemtx(int **m, int i)
{
    int c;
    for (c = 0; c < i; c++)
        if (c < i)
            free(m[c]);
    free(m);
}

/*!
    Função para verificar se a cadeia c contem apenas os caracteres representando os Nucleotideos.
*/
int p1_valida_cadeia (char *c)
{
    if (c)
    {
        int l = strlen(c);
        int i;
        for (i = 0; i < l; i++)
        {
            char nuc = c[i];
            if (nuc != 'a' && nuc != 'c' && nuc != 'g' && nuc != 't' && nuc != 'u')
                return 0;
        }
        return 1;
    }
    return 0;
}

/*!
    Função para encontrar o melhor alinhamento entre as cadeias s e t,
    dado que s e t sejam cadeias validas de nucleotideos (a,c,g,t).

    @param s ponteiro para a cadeia s
    @param t ponteiro para a cadeia t, a ser alinhada com s
    @param als referencia para retornar a cadeia s alinhada
    @param alt referencia para retornar a cadeia t alinhada
*/
void p1_alinhar_s_t (char *s, char *t, char **als, char **alt)
{
    int slen, tlen, lignlen;
    int i, j, k, maxij, diagv, esqv, cimav, globalscr;
    char dir;
    int **m, **mdir;
    char *slign, *tlign;
    int mrows, mcols;

    /*Alocação e  inicialização da matriz de scores*/
    slen = strlen(s);
    tlen = strlen(t);

    mrows = tlen + 1;
    mcols = slen + 1;

    m = allocintmtx(mrows,mcols);
    mdir = allocintmtx(mrows,mcols);

    for (i = 1; i < mcols; i++)
    {
        m[0][i] = i * SCR_R;
        mdir[0][i] = DIRESQ;
    }

    for (i = 1; i < mrows; i++)
    {
        m[i][0] = i * SCR_I;
        mdir[i][0] = DIRSUB;
    }

    for (i = 1; i < mrows; i++)
        for (j = 1; j < mcols; j++)
        {
            esqv = m[i][j-1] + SCR_R;
            diagv = (s[j-1] == t[i-1]) ? (m[i-1][j-1] + SCR_M) : (m[i-1][j-1] + SCR_S);
            cimav = m[i-1][j] + SCR_I;

            maxij = max3(esqv,diagv,cimav);

            dir = 0;

            if (maxij == esqv)
                dir = DIRESQ;
            if (maxij == cimav)
                dir = dir | DIRSUB;
            if (maxij == diagv)
                dir = dir | DIRDIA;

            mdir[i][j] = dir;
            m[i][j] = maxij;
        }

    /*Alocação das strings s e t alinhadas*/

    lignlen = mrows + mcols;

    slign = (char *) calloc(lignlen + 1, sizeof(char));
    tlign = (char *) calloc(lignlen + 1, sizeof(char));

    i = mrows - 1;
    j = mcols - 1;

    slign[lignlen] = '\0';
    tlign[lignlen] = '\0';

    k = lignlen - 1;

    /*Traceback*/

    globalscr = 0;

    do
    {
        dir = mdir[i][j];
        globalscr += m[i][j];

        /*Caminha no traceback segundo a prioridade das operacoes*/
        if (goDG(dir))
        {
            slign[k] = s[j-1];
            tlign[k] = t[i-1];
            i--;
            j--;
            k--;
        }
        else if (goLF(dir))
        {
            slign[k] = s[j-1];
            tlign[k] = GAP;
            j--;
            k--;
        }
        else if (goUP(dir))
        {
            slign[k] = GAP;
            tlign[k] = t[i-1];
            i--;
            k--;
        }
    }
    while(dir);

    if (k+1)
    {
        strcpy(slign,slign + k + 1);
        strcpy(tlign,tlign + k + 1);
    }

    printf("\n%s\n%s\nscore: %d\n",slign,tlign,globalscr);

    if (als && alt)
    {
        *als = slign;
        *alt = tlign;
    }
    else
    {
        free(slign);
        free(tlign);
    }

    /*Codigo para imprimir as matrizes m e mdir (mdir: mapa de traceback).*/

//    for (i = 0; i < mrows; i++)
//    {
//        for (j = 0; j < mcols; j++)
//        {
//            printf("%03d ", m[i][j]);
//        }
//        printf("\n");
//    }
//    printf("\n");
//    for (i = 0; i < mrows; i++)
//    {
//        for (j = 0; j < mcols; j++)
//        {
//            printf("%03d ", mdir[i][j]);
//        }
//        printf("\n");
//    }

    freemtx(mdir,mrows);
    freemtx(m,mrows);
}

int main (int argc, char *argv[])
{
    if (argc > 2)
    {
        printf("P1 Console\n");
        if (p1_valida_cadeia(argv[1]) && p1_valida_cadeia(argv[2]))
        {
            printf("Cadeias OK!\n");
            p1_alinhar_s_t(argv[1],argv[2],NULL,NULL);
        }
        else
        {
            printf("Cadeias Invalidas!\n");
        }
    }
    else
    {//Executa exemplo para depuracao
        p1_alinhar_s_t("atggcatatcccatacaactaggattccaagatgcaacatcaccaatcatagaaga","cacttcatcaagaagatactaaccactacaacgtagaaccttagga",NULL,NULL);
    }
    return 0;
}
