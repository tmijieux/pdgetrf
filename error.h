#ifndef PDGETRF_ERROR_H
#define PDGETRF_ERROR_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>


#define PDGETRF_ERROR_LIST(ERROR)                  \
    ERROR(SUCCESS)                              \


#define PDGETRF_ERROR_TO_ENUM(error) PDGETRF_ERROR_##error,

enum matprod_error_code {
    PDGETRF_ERROR_LIST(PDGETRF_ERROR_TO_ENUM)
};

char const *matprod_errmsg(int errcode);

#ifndef __GNUC__
#define __PRETTY_FUNCTION__    __FUNCDNAME__
#endif


#define __FILENAME__ (strrchr(__FILE__, '/') ?                  \
                      strrchr(__FILE__, '/') + 1 : __FILE__)

#define matprod_error(format_, ...)                                \
    do {                                                        \
        fprintf(stderr, "ERROR: %s:%d|%s: ",  __FILENAME__ ,    \
                __LINE__, __PRETTY_FUNCTION__);                 \
        fprintf(stderr, (format_), ##__VA_ARGS__);              \
    } while(0)

#define matprod_fatal(format_, ...)                                        \
    do {                                                                \
        fprintf(stderr, "FATAL ERROR: %s:%d|%s: ",  __FILENAME__ ,      \
                __LINE__, __PRETTY_FUNCTION__);                         \
        fprintf(stderr, (format_), ##__VA_ARGS__);                      \
        exit(EXIT_FAILURE);                                             \
    } while(0)

#define matprod_warning(format_, ...)                              \
    do {                                                        \
        fprintf(stderr, "WARNING: %s:%d|%s: ",  __FILENAME__ ,  \
                __LINE__, __PRETTY_FUNCTION__);                 \
        fprintf(stderr, (format_), ##__VA_ARGS__);              \
    } while(0)

#ifdef DEBUG
#define matprod_debug(format_, ...)                                        \
    do {                                                                \
        fprintf(stderr, "DEBUG: %s:%d|%s: " format_, __FILENAME__ ,     \
                __LINE__, __PRETTY_FUNCTION__, ##__VA_ARGS__);          \
    } while(0)

#else // DEBUG
#define matprod_debug(format_, ...) ((void) (format_))
#endif // DEBUG

#endif // PDGETRF_ERROR_H
