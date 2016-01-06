#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "rest_aln.h"

int rest_aln_usage(void)
{
    fprintf(stderr, "rest aln [option] ref.fa read.fq/fa\n");
    return 0;
}

rest_aln_para *rest_init_ap(void)
{
    rest_aln_para *ap = (rest_aln_para*)calloc(1, sizeof(rest_aln_para));
    ap->seed_len = REST_SEED_LEN;
    return ap;
}

int rest_aln_core(const char *ref_fn, const char *read_fn, rest_aln_para *rest_ap)
{
    return 0;
}

int rest_aln(int argc, char *argv[])
{
    int c;
    rest_aln_para *rest_ap = rest_init_ap();

    while ((c = getopt(argc, argv, "l:")) >= 0) {
        switch (c)
        {
            case 'l': rest_ap->seed_len = atoi(optarg); break;
            default: return rest_aln_usage();
        }
    }
    if (argc - optind != 3) return rest_aln_usage();

    rest_aln_core(argv[optind+1], argv[optind+2], rest_ap);
    return 0;
}
