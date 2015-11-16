#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef struct Arvore
{
	unsigned o;
	unsigned o2;
	unsigned p;
	struct Arvore *e;
	struct Arvore *d;
}Arvore;

typedef struct Conjunto
{
	unsigned o;
	unsigned rnk;
	struct Conjunto *r;
}Conjunto;

typedef struct Altura
{
	unsigned o1;
	unsigned o2;
	float alt;
}Altura;

typedef struct AU
{
	unsigned o1;
	unsigned o2;
	float alt;
	float pe;
	float pd;
	struct AU *e;
	struct AU *d;
}AU;

/*
	Lista encadeada para a arvore geradora de peso minimo
*/
typedef struct MST
{
	unsigned o1; //Obj 1 da aresta
	unsigned o2; //Obj 2 da aresta
	unsigned p; //peso da aresta
	unsigned cw; //cut-weight
	struct MST *prx; //ref da proxima aresta
}MST;

typedef struct heap
{
	int d;
	int v[2];
}heap;

void troca(heap* a, heap* b)
{
	int c = a->d, v0 = a->v[0], v1 = a->v[1];
	a->d = b->d;
	a->v[0] = b->v[0];
	a->v[1] = b->v[1];
	b->d = c;
	b->v[0] = v0;
	b->v[1] = v1;
}
void heapify(heap* q, int i, int n)
{
	int l = 2 * i + 1, r = 2 * (i + 1), menor = i;
	if (l<n && q[l].d > q[i].d) menor = l;
	if (r<n && q[r].d > q[menor].d) menor = r;
	if (menor != i)
	{
		troca(q + i, q + menor);
		heapify(q, menor, n);
	}
}
void buildheap(heap* vt, int n)
{
	int i;
	for (i = (n - 1) / 2; i >= 0; i--) heapify(vt, i, n);
}

void heapsort(heap *vt, int n, int *h)
{
	int i;
	buildheap(vt, n);
	for (i = n - 1; i >= 0; i--) {
		heap aux;
		aux.d = vt[0].d;
		aux.v[0] = vt[0].v[0];
		aux.v[1] = vt[0].v[1];
		vt[0].d = vt[i].d;
		vt[0].v[0] = vt[i].v[0];
		vt[0].v[1] = vt[i].v[1];
		vt[i].d = aux.d;
		vt[i].v[0] = aux.v[0];
		vt[i].v[1] = aux.v[1];
		heapify(vt, 0, i);
	}
}

int **intquadmtx(unsigned n)
{
	int **m;
	unsigned i;

	m = (int **) malloc(n * sizeof(int*));

	for (i = 0; i < n; i++)
		m[i] = (int *) calloc(n,sizeof(int));

	return m;
}

MST* insmst(MST *h, unsigned i, unsigned j, unsigned p)
{
	MST *novo = NULL, *aux = h;
	if (!h)
	{
		h = (MST*) malloc(sizeof(MST));
		h->o1 = i;
		h->o2 = j;
		h->p = p;
		h->cw = 0;
		h->prx = NULL;
		return h;
	}

	novo = (MST*)malloc(sizeof(MST));
	novo->o1 = i;
	novo->o2 = j;
	novo->p = p;
	novo->cw = 0;
	novo->prx = NULL;

	while (aux->prx) aux = aux->prx;
	aux->prx = novo;

	return h;
}

int conectada(int v[2], MST *T)
{
	int con1, con2;
	
	if (!T)
		return 0;

	con1 = con2 = 0;
	while (T && (!con1 || !con2))
	{
		if (T->o1 == v[0] || T->o2 == v[0])
			con1 = 1;
		if (T->o1 == v[1] || T->o2 == v[1])
			con2 = 1;
		T = T->prx;
	}

	if (con1 && con2)
		return 1;
	return 0;
}

MST* criamst(int** g, unsigned n)
{
	MST *T = NULL;
	heap *ordem;
	int tamheap;
	unsigned i, j, k = 0;

	ordem = (heap *)calloc((n*n - n)/2,sizeof(heap));

	for (i = 0; i < n; i++)
		for (j = i; j < n; j++)
			if (i != j)
			{
				ordem[k].d = g[i][j];
				ordem[k].v[0] = i;
				ordem[k].v[1] = j;
				k++;
			}
	tamheap = 0;
	heapsort(ordem, (n*n - n) / 2, &tamheap);

	k = (n*n - n) / 2;
	for (i = 0; i < k; i++)
		if (!conectada(ordem[i].v,T))
			T = insmst(T, ordem[i].v[0], ordem[i].v[1], ordem[i].d);

	free(ordem);
	return T;
}

Conjunto* makeSet(unsigned x)
{
	Conjunto *c;
	c = (Conjunto*)malloc(sizeof(Conjunto));
	c->o = x;
	c->rnk = 1;
	c->r = c;
	return c;
}

Conjunto* findSet(Conjunto *c)
{
	if (c != c->r)
		c->r = findSet(c->r);
	return c->r;
}

void unionSet(Conjunto *a, Conjunto *b)
{
	Conjunto *la, *lb;

	la = findSet(a);
	lb = findSet(b);
	
	if (lb->rnk > la->rnk)
	{
		la->r = lb;
		lb->rnk += la->rnk;
	}
	else
	{
		lb->r = la;
		la->rnk = lb->rnk;
	}
}

void iniciaFloresta(Conjunto **f, unsigned n)
{
	unsigned i;

	for (i = 0; i < n; i++)
		f[i] = makeSet(i);
}

void atualizaArvore(Arvore *t, Arvore **folhas, Arvore *root)
{
	if (t)
	{
		if (!(t->e || t->d))
		{
			folhas[t->o] = root;
			return;
		}
		atualizaArvore(t->e, folhas, root);
		atualizaArvore(t->d, folhas, root);
	}
}

Arvore* criaR(MST *T, Conjunto **conjuntos, unsigned n)
{
	Arvore **folhas, *R;
	MST *aresta = T;
	unsigned i;

	folhas = (Arvore**)calloc(n, sizeof(Arvore*));
	for (i = 0; i < n; i++)
	{
		folhas[i] = (Arvore*)calloc(1, sizeof(Arvore));
		folhas[i]->o = i;
		folhas[i]->o2 = i;
		folhas[i]->p = 0;
		folhas[i]->e = folhas[i]->d = NULL;
	}

	while (aresta)
	{
		Conjunto *A, *B;
		A = findSet(conjuntos[aresta->o1]);
		B = findSet(conjuntos[aresta->o2]);
		if (A != B)
		{
			Arvore *r, *ra, *rb;
			ra = folhas[A->o];
			rb = folhas[B->o];
			r = (Arvore*)calloc(1, sizeof(Arvore));
			r->o = aresta->o1;
			r->o2 = aresta->o2;
			r->p = aresta->p;
			r->e = ra;
			r->d = rb;
			unionSet(A, B);
			atualizaArvore(r, folhas, r);
			R = r;
		}
		aresta = aresta->prx;
	}
	return R;
}

Arvore* lca(Arvore ***caminhos, unsigned o1, unsigned o2)
{
	Arvore **c1 = caminhos[o1], **c2 = caminhos[o2];
	Arvore *lca = c1[0];
	unsigned i = 0;

	while (c1[i] == c2[i])
	{
		lca = c1[i];
		i++;
	}
	return lca;
}

void preOrd(Arvore *r, Arvore ***caminhos, Arvore **cam, int l)
{
	cam[l] = r;
	if (!r->e && !r->d)
	{
		unsigned i;
		for (i = 0; i <= l; i++)
			caminhos[r->o][i] = cam[i];
		return;
	}
	preOrd(r->e, caminhos, cam, l + 1);
	preOrd(r->d, caminhos, cam, l + 1);
}

unsigned** calcCW(int **L, MST *T, Arvore *R, unsigned n)
{
	Arvore ***caminhos, **cam;
	unsigned i,h = log2(n) + 1;
	MST *aux = T;
	unsigned **CW;

	caminhos = (Arvore***)calloc(n,sizeof(Arvore**));
	CW = (unsigned**)calloc(n, sizeof(unsigned*));
	for (i = 0; i < n; i++)
	{
		caminhos[i] = (Arvore**)calloc(h, sizeof(Arvore*));
		CW[i] = (unsigned*)calloc(n, sizeof(unsigned));
	}

	/*Captura os caminhos ate cada folha em caminhos*/
	cam = (Arvore**)calloc(h, sizeof(Arvore*));
	preOrd(R, caminhos, cam, 0);

	/*Inicializa CW*/
	for (i = 0; i < n; i++)
		for (h = i; h < n; h++)
			CW[i][h] = 0;

	for (i = 0; i < n; i++)
		for (h = i; h < n; h++)
		{
			if (i != h)
			{
				Arvore *e = lca(caminhos, i, h);
				if (L[i][h] > CW[e->o][e->o2])
					CW[e->o][e->o2] = L[i][h];
			}
		}

	return CW;
}

void atualizaArvoreU(AU *t, AU **folhas, AU *root)
{
	if (t)
	{
		if (!(t->e || t->d))
		{
			folhas[t->o1] = root;
			return;
		}
		atualizaArvoreU(t->e, folhas, root);
		atualizaArvoreU(t->d, folhas, root);
	}
}

Arvore* criaU(Arvore *R, heap *T, unsigned **CW, unsigned n)
{
	Conjunto **floresta;
	AU **folhas;
	unsigned i, na = n - 1;

	floresta = (Conjunto**)calloc(n,sizeof(Conjunto*));
	folhas = (AU**)calloc(n, sizeof(AU*));
	iniciaFloresta(floresta, n);
	for (i = 0; i < n; i++)
	{
		folhas[i] = (Arvore*)calloc(1, sizeof(Arvore));
		folhas[i]->o1 = i;
		folhas[i]->o2 = i;
		folhas[i]->alt = 0;
		folhas[i]->pd = 0;
		folhas[i]->pe = 0;
		folhas[i]->e = folhas[i]->d = NULL;
	}
	for (i = 0; i < na; i++)
	{
		Conjunto *A = findSet(floresta[T[i].v[0]]);
		Conjunto *B = findSet(floresta[T[i].v[1]]);

		if (A != B)
		{
			AU *ua = folhas[A->o], *ub = folhas[B->o];
			AU *u = (Arvore*)malloc(sizeof(Arvore));
			u->o1 = T[i].v[0];
			u->o2 = T[i].v[1];
			u->e = ua;
			u->d = ub;
			u->alt = (float) (CW[T[i].v[0]][T[i].v[1]] / 2);
			u->pe = u->alt - ua->alt;
			u->pd = u->alt - ub->alt;
			atualizaArvoreU(u, folhas, u);
			unionSet(A, B);
		}
	}
	return folhas[0];
}

void printmst(MST *t)
{
	printf("vertices: v1, v2,...vN\n");
	printf("MST:\n");
	while (t)
	{
		printf("v%d-v%d: %u\n", t->o1,t->o2,t->p);
		t = t->prx;
	}
	printf("\n");
}

void printr(Arvore *r)
{
	if (!r) return;
	printr(r->e);
	printf("v%d-v%d ",r->o,r->o2);
	printr(r->d);
}

void printcw(unsigned **cw, MST *T)
{
	printf("CW\n");
	while (T)
	{
		printf("v%u-v%u: %u\n",T->o1,T->o2,cw[T->o1][T->o2]);
		T = T->prx;
	}
}

void printu(AU *u)
{
	if (!u)return;
	printu(u->e);
	printf("v%u-v%u:%.02f:%.02f ",u->o1,u->o2,u->pe,u->pd);
	printu(u->d);
}

int main(int argc, char *argv[])
{
	char fname[100];
	int **H, **L;
	unsigned n, i, j, **CW;
	MST *T;
	Conjunto **floresta;
	Arvore *R;
	AU *U;
	heap *T2;

	if (argc > 1)
		strcpy(fname, argv[1]);
	else
		strcpy(fname, "in.txt");

	if (!freopen(fname, "r", stdin))
		return 0x1;
	freopen("djfm_UT_Outputs.TXT","w",stdout);

	while (scanf("%u", &n) == 1) {

		H = intquadmtx(n);
		L = intquadmtx(n);

		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				scanf("%d", (L[i]) + j);
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
			{
				scanf("%d", (H[i]) + j);
				if (H[i][j] < L[i][j])
				{
					printf("N�o � poss�vel calcular uma �rvore ultram�trica obedecendo aos limites dados.\n");
					return 0x2;
				}
			}

		printf("N = %u\n\n", n);

		/* Contrucao MST T */
		T = criamst(H, n);

		printmst(T);

		/*Contrucao de R*/
		floresta = (Conjunto**)calloc(n, sizeof(Conjunto*));
		iniciaFloresta(floresta, n);
		R = criaR(T, floresta, n);

		printf("Arvore R em pre ordem:\n");
		printr(R);
		printf("\n");

		/*Cut-weight*/
		CW = calcCW(L, T, R, n);

		printcw(CW, T);

		/*Construcao U*/
		T2 = (heap*)calloc(n - 1, sizeof(heap));
		i = 0;
		while (T)
		{
			MST *x = T->prx;
			T2[i].v[0] = T->o1;
			T2[i].v[1] = T->o2;
			T2[i].d = CW[T->o1][T->o2];
			free(T);
			T = x;
			i++;
		}
		heapsort(T2, n - 1, &i);

		/*Criacao de U*/
		U = criaU(R, T2, CW, n);

		printf("\nArvore Ultrametrica U, em pre ordem (aresta:pesoEsq:pesoDir):\n");
		printu(U);

		/*Fim*/
		for (i = 0; i < n; i++)
		{
			free(L[i]);
			free(H[i]);
		}

		free(L);
		free(H);

		printf("\n\n");
	}

	fclose(stdin);
	fclose(stdout);
    return 0;
}
