#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "fasta.h"
#include "matrices.h"

/*!
    @author douglas medeiros <djfm at cin.ufpe.br>

    1) Instrucoes:

        comando:
            ./P1 <cadeia 1> <cadeia 2> <tabela>
            Onde <tabela> pode ser pam10 (default), pam120, pam250, blosum30, blosum62, blosum90;
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

/*Direções*/
#define DIRSUB 0x01
#define DIRESQ 0x02
#define DIRDIA 0x04

#define DIR_A 0x20
#define DIR_B 0x40
#define DIR_C 0x80

#define goUP(x) (x & DIRSUB)
#define goLF(x) (x & DIRESQ)
#define goDG(x) (x & DIRDIA)

#define goA(x) (x & DIR_A)
#define goB(x) (x & DIR_B)
#define goC(x) (x & DIR_C)

#define intmax(a,b) ((a > b) ? a : b)

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

int escolherTabela(char *opt)
{
    if (strcpy(opt,"pam120") == 0)
    {
        printf("Tabela: pam120\n");
        return 1;
    }
    if (strcpy(opt,"pam250") == 0)
    {
        printf("Tabela: pam250\n");
        return 2;
    }
    if (strcpy(opt,"blosum30") == 0)
    {
        printf("Tabela: blosum30\n");
        return 3;
    }
    if (strcpy(opt,"blosum62") == 0)
    {
        printf("Tabela: blosum62\n");
        return 4;
    }
    if (strcpy(opt,"blosum90") == 0)
    {
        printf("Tabela: blosum90\n");
        return 5;
    }

    printf("Tabela default: pam10\n");
    return 0;
}

int gappen (int k, _int8 tab[24][24])
{
    int g = tab[0][23];
    if (k > 1)
        return (2 * g + g * (k - 1));
    if (k == 1)
        return 2 * g;
    return 0;
}

/*!
    Função para mapear os caracteres que representam os aminoacidos nos valores do tipo enumeracao Aminoacidos_t
    @param a Caracter representando um aminoacido, espera-se que seja minusculo
    @return Aminoacido_t valor da enumeração correspondente
*/
Aminoacidos_t getAminCode (char a)
{
    Aminoacidos_t aminoCode;
    switch(a)
    {
    case 'a' :
        aminoCode = A;
        break;
    case 'r' :
        aminoCode = R;
        break;
    case 'n' :
        aminoCode = N;
        break;
    case 'd' :
        aminoCode = D;
        break;
    case 'c' :
        aminoCode = C;
        break;
    case 'q' :
        aminoCode = Q;
        break;
    case 'e' :
        aminoCode = E;
        break;
    case 'g' :
        aminoCode = G;
        break;
    case 'h' :
        aminoCode = H;
        break;
    case 'i' :
        aminoCode = I;
        break;
    case 'l' :
        aminoCode = L;
        break;
    case 'k' :
        aminoCode = K;
        break;
    case 'm' :
        aminoCode = M;
        break;
    case 'f' :
        aminoCode = F;
        break;
    case 'p' :
        aminoCode = P;
        break;
    case 's' :
        aminoCode = S;
        break;
    case 't' :
        aminoCode = T;
        break;
    case 'w' :
        aminoCode = W;
        break;
    case 'y' :
        aminoCode = Y;
        break;
    case 'v' :
        aminoCode = V;
        break;
    case 'b' :
        aminoCode = B;
        break;
    case 'z' :
        aminoCode = Z;
        break;
    case 'x' :
        aminoCode = X;
        break;
    default  :
        aminoCode = O;
        break;
    }
    return aminoCode;
}

/*!
    Função para validar se a cadeia c contem apenas os caracteres representando os aminoacidos em minusculo
    @param c Ponteiro para a cadeia a ser validada
    @return um se a cadeia é válida, zero se não
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
            if (!((nuc >= 'a' && nuc <= 'z') || (nuc >= 'A' && nuc <= 'Z')))
                return 0;
            if (nuc >= 'A' && nuc <= 'Z')
                c[i] = nuc + 32;
        }
        return 1;
    }
    return 0;
}

int maxmtx (int **arr, int m, int n, int* i, int* j)
{
    int k,l,bf = 0;
    *i = 0;
    *j = 0;
    for (k = 0; k < m; k++)
        for (l = 0; l < n; l++)
            if (arr[k][l] > bf)
            {
                bf = arr[k][l];
                *i = k;
                *j = l;
            }
    return bf;
}

/*!
    Função para encontrar o melhor alinhamento local entre as cadeias s e t,
    dado que s e t sejam cadeias validas de aminoacidos.

    @param s ponteiro para a cadeia s
    @param t ponteiro para a cadeia t, a ser alinhada com s
    @param als referencia para retornar a cadeia s alinhada
    @param alt referencia para retornar a cadeia t alinhada
    @param tab Matriz 24x24 representando a funcao de substituicao (pam10, blosum30...)
*/
void p1_alinhar_s_t (char *s, char *t, char **als, char **alt, _int8 tab[24][24])
{
    int slen, tlen, lignlen;
    int i, j, k, x, y, maxij, diagv, esqv, cimav, globalscr;
    char dir;
    int **a, **b, **c, **mdir;
    char *slign, *tlign;
    int mrows, mcols;

    /* Alocação e  inicialização das matrizes a,b e c */
    /* Obs.: calloc() ja inicializa com zero, conveniente para o caso de alinhamento local apenas */

    slen = strlen(s);
    tlen = strlen(t);

    mrows = tlen + 1;
    mcols = slen + 1;

    a = allocintmtx(mrows,mcols);
    b = allocintmtx(mrows,mcols);
    c = allocintmtx(mrows,mcols);
    mdir = allocintmtx(mrows,mcols);

    /* iteracao */
    i = 0; //Indice da cadeia t
    j = 0; //Indice da cadeia s
    k = 0; //Indice para blocos maiores que 1
    x = y = 1; //Indices para as matrizes

    while (x < mrows && y < mcols)
    {
        a[x][y] = tab[getAminCode(s[i]) + 0][getAminCode(t[j]) + 0] + intmax(a[x-1][y-1],intmax(b[x-1][y-1],c[x-1][y-1]));

        k = 0;
        b[x][y] = intmax(a[x][y - k] + gappen(k,tab),c[x][y - k] + gappen(k,tab));
        for (k = 1; k <= y; k++)
            if (((a[x][y - k] + gappen(k,tab)) > b[x][y]) || ((c[x][y - k] + gappen(k,tab)) > b[x][y]))
                b[x][y] = intmax(a[x][y - k],c[x][y - k]);

        k = 0;
        c[x][y] = intmax(a[x - k][y] + gappen(k,tab),b[x - k][y] + gappen(k,tab));
        for (k = 1; k <= y; k++)
            if (((a[x - k][y] + gappen(k,tab)) > c[x][y]) || ((b[x - k][y] + gappen(k,tab)) > c[x][y]))
                c[x][y] = intmax(a[x - k][y],b[x - k][y]);

        maxij = intmax(a[x][y],intmax(b[x][y],c[x][y]));

        if (maxij == a[x][y]) mdir[x][y] = DIR_A;
        else if (maxij == b[x][y]) mdir[x][y] = DIR_B;
        else mdir[x][y] = DIR_C;

        x++;y++;
        i++;j++;

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
        globalscr += a[i][j];

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
    freemtx(a,mrows);
    freemtx(b,mrows);
    freemtx(c,mrows);
}

int main (int argc, char *argv[])
{
    if (argc > 2)
    {
        FASTA *fasta1, *fasta2;
        char *cadeia1, *cadeia2;
        unsigned len1, len2;

        fasta1 = abrirFasta(argv[1]);
        fasta2 = abrirFasta(argv[2]);

        lerFasta(fasta1,&cadeia1,NULL,&len1);
        lerFasta(fasta2,&cadeia2,NULL,&len2);

        printf("P1 Console\n");
        if (p1_valida_cadeia(cadeia1) && p1_valida_cadeia(cadeia2))
        {
            int subst = 0;
            printf("Cadeias OK!\n");
            if (argc > 3) subst = escolherTabela(argv[3]);
            p1_alinhar_s_t(cadeia1,cadeia2,NULL,NULL,tab[subst]);
        }
        else
        {
            printf("Cadeias Invalidas!\n");
        }

        encerraFasta(fasta1);
        encerraFasta(fasta2);
        free(cadeia1);
        free(cadeia2);
    }
    else
    {
        /* Executa exemplo para depuracao */
        printf("%d",pam10[23][0]);
        p1_alinhar_s_t("atggcatatcccatacaactaggattccaagatgcaacatcaccaatcatagaaga","cacttcatcaagaagatactaaccactacaacgtagaaccttagga",NULL,NULL,tab[0]);
    }
    return 0;
}
