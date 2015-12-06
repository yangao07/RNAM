#ifndef _ERROR_H_
#define _ERROR_H_

gzFile err_gzopen(const char *func, char *file, char *m);
void err_gzclose(const char *func, gzFile fp);


#endif
