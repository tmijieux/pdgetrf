#include "error.h"

#define PDGETRF_ERROR_TO_STRING(error) [PDGETRF_ERROR_##error] = #error,

static char const *matprod_errmsg_[] = {
    PDGETRF_ERROR_LIST(PDGETRF_ERROR_TO_STRING)
};

char const *matprod_errmsg(int error_code)
{
    return matprod_errmsg_[error_code > 0 ? error_code : -error_code];
}
