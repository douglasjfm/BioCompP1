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
#define DIRDIA 0x02
#define DIRESQ 0x04

#define DIR_A 0xA0
#define DIR_B 0xB0
#define DIR_C 0xC0

#define dirmat(a,b,c) ((a > b) ? (a > c ? DIR_A : DIR_C) : (b > c ? DIR_B : DIR_C))

#define goUP(x) (x & DIRSUB)
#define goLF(x) (x & DIRESQ)
#define goDG(x) (x & DIRDIA)

#define goA(x) (x & DIR_A)
#define goB(x) (x & DIR_B)
#define goC(x) (x & DIR_C)

#define intmax(a,b) ((a > b) ? a : b)

#define piso -99//0x80000001

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
            if (arr[k][l] >= bf)
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
    int i, j, x, y, maxij, h, g;
    int score, maxa, ia, ja, maxb, ib, jb, maxc, ic, jc;
    int **a, **b, **c, **tracem;
    unsigned **mdir;
    char *slign, *tlign;
    int mrows, mcols;

    /* Alocação e  inicialização das matrizes a,b e c */

    slen = strlen(s);
    tlen = strlen(t);

    mrows = tlen + 1;
    mcols = slen + 1;

    a = allocintmtx(mrows,mcols);
    b = allocintmtx(mrows,mcols);
    c = allocintmtx(mrows,mcols);
    mdir = (unsigned **)allocintmtx(mrows,mcols);

    a[0][0] = 0;

    b[0][0] = piso;
    for(i = 1; i < mrows; i++)
    {
        a[i][0] = b[i][0] = piso;
        c[i][0] = 0;
    }

    c[0][0] = piso;
    for(j = 1; j < mcols; j++)
    {
        b[0][j] = 0;
        a[0][j] = c[0][j] = piso;
    }

    /* iteracao */

    i = 0; //Indice da cadeia t
    j = 0; //Indice da cadeia s
    x = y = 1; //Indices para as matrizes
    h = tab[0][23]; //Preco de abertura do gap
    g = tab[0][23]; //Preco da extensao de gap

    for(x = 1; x < mrows; x++)
        for(y = 1; y < mcols; y++)
        {
            a[x][y] = tab[getAminCode(s[i]) + 0][getAminCode(t[j]) + 0] + intmax(a[x-1][y-1],intmax(b[x-1][y-1],c[x-1][y-1]));
            b[x][y] = intmax(intmax(h + g + a[x][y-1], g + b[x][y-1]),h + g + c[x][y-1]);
            c[x][y] = intmax(intmax(h + g + a[x-1][y], g + c[x-1][y]),h + g + b[x-1][y]);

            maxij = intmax(a[x][y],intmax(b[x][y],c[x][y]));

            if (maxij == a[x][y])
            {
                _int8 mat = dirmat(a[x-1][y-1],b[x-1][y-1],c[x-1][y-1]);
                mdir[x][y] = (DIRDIA | mat) & 0xff;
            }
            else if (maxij == b[x][y])
            {
                _int8 mat = dirmat(h + a[x][y-1],b[x][y-1],h + c[x][y-1]);
                mdir[x][y] = (DIRESQ | mat) & 0xff;
            }
            else
            {
                _int8 mat = dirmat(h + a[x-1][y],h + b[x-1][y],c[x-1][y]);
                mdir[x][y] = (DIRSUB | mat) & 0xff;
            }

            i = x - 1;
            j = y - 1;
        }

    /*Alocação das strings s e t alinhadas*/

    lignlen = mrows + mcols;

    slign = (char *) calloc(lignlen + 1, sizeof(char));
    tlign = (char *) calloc(lignlen + 1, sizeof(char));

    i = mrows - 1;
    j = mcols - 1;

    slign[lignlen] = '\0';
    tlign[lignlen] = '\0';

    /*Traceback*/

    maxa = maxmtx(a,mrows,mcols,&ia,&ja);
    maxb = maxmtx(b,mrows,mcols,&ib,&jb);
    maxc = maxmtx(c,mrows,mcols,&ic,&jc);

    tracem = (maxa > maxb) ? (maxa > maxc ? a : c) : (maxb > maxc ? b : c);
    i = (tracem == a) ? ia : (tracem == b ? ib : ic);
    j = (tracem == a) ? ja : (tracem == b ? jb : jc);

    x = y = lignlen - 1;
    score = 0;
    while(tracem[i][j] != 0 && i > 0 && j > 0)
    {
        score += tracem[i][j];

        if (goDG(mdir[i][j]))
        {
            slign[x] = s[j-1];
            tlign[y] = t[i-1];
            i--;
            j--;
            if (goA(mdir[i][j])) tracem = a;
            else if (goB(mdir[i][j])) tracem = b;
            else tracem = c;
        }
        else if (goUP(mdir[i][j]))
        {
            slign[x] = GAP;
            tlign[y] = t[i-1];
            i--;
            if (goA(mdir[i][j])) tracem = a;
            else if (goB(mdir[i][j])) tracem = b;
            else tracem = c;
        }
        else
        {
            slign[x] = s[j-1];
            tlign[y] = GAP;
            j--;
            if (goA(mdir[i][j])) tracem = a;
            else if (goB(mdir[i][j])) tracem = b;
            else tracem = c;
        }
        x--;
        y--;
    }

    strcpy(slign,slign + x + 1);
    strcpy(tlign,tlign + y + 1);

    printf("%s\n%s\nscore: %d\n",slign,tlign,score);

    /*Codigo para imprimir as matrizes m e mdir (mdir: mapa de traceback).*/

    for (i = 0; i < mrows; i++)
    {
        for (j = 0; j < mcols; j++)
        {
            printf("%03d ", a[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    for (i = 0; i < mrows; i++)
    {
        for (j = 0; j < mcols; j++)
        {
            printf("%03d ", b[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    for (i = 0; i < mrows; i++)
    {
        for (j = 0; j < mcols; j++)
        {
            printf("%03d ", c[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    for (i = 0; i < mrows; i++)
    {
        for (j = 0; j < mcols; j++)
        {
            printf("%02x ", mdir[i][j]);
        }
        printf("\n");
    }

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
            p1_alinhar_s_t(cadeia1,cadeia2,NULL,NULL,tabela[subst]);
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
        //p1_alinhar_s_t("atggcatatcccatacaactaggattccaagatgcaacatcaccaatcatagaaga","cacttcatcaagaagatactaaccactacaacgtagaaccttagga",NULL,NULL,tab[0]);
        p1_alinhar_s_t("send","end",NULL,NULL,pam10);
    }
    return 0;
}
