#include <stdio.h>
#include <stdlib.h>

int rest_aln_usage(void)
{
    fprintf(stderr, "rest aln [option] ref.fa read.fq/fa\n");
    return 0;
}

int rest_aln(int argc, char *argv[])
{
    return rest_aln_usage();
    return 0;
}

