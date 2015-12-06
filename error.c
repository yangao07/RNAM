#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <errno.h>

int err_exit(const char *header, const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    fprintf(stderr, "[%s] ", header);
    vfprintf(stderr, fmt, args);
    fprintf(stderr, "\n");
    va_end(args);
    exit(EXIT_FAILURE);
}

gzFile err_gzopen(const char *func, char *file, char *m)
{
    gzFile fp;
    if ((fp = gzopen(file, m)) == 0) {
        err_exit(func, "fail to open file '%s'\n", file);
    } else return fp;
}

void err_gzclose(const char *func, gzFile fp)
{
    int ret = gzclose(fp);
    if (ret != Z_OK) {
        err_exit(func, "gzclose fail %s", Z_ERRNO == ret ? strerror(errno) : zError(ret));
    }
}
