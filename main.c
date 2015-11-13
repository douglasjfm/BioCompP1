#include <stdlib.h>
#include <stdio.h>
#include <string.h>

char **charquadmtx(unsigned n)
{
	char **m;
	unsigned i;

	m = (char **) malloc(n * sizeof(char*));

	for (i = 0; i < n; i++)
		m[i] = (char *) calloc(n,sizeof(char));

	return m;
}

int main (int argc, char *argv[])
{
	char fname[100];
	char **H, **L;
	unsigned n,i,j;

	system("dir > test.txt");

	if (argc > 1)
		strcpy(fname, argv[1]);
	else
		strcpy(fname,"in.txt");

	if (!freopen(fname,"r", stdin))
		exit(0x1);

	scanf("%u",&n);

	H = charquadmtx(n);
	L = charquadmtx(n);

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			scanf("%d",(L[i])+j);
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			scanf("%d", (H[i]) + j);

	printf("ok");
	fclose(stdin);
    return 0;
}
