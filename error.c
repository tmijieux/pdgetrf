#include "error.h"

#define MATPROD_ERROR_TO_STRING(error) [MATPROD_ERROR_##error] = #error,

static char const *matprod_errmsg_[] = {
    MATPROD_ERROR_LIST(MATPROD_ERROR_TO_STRING)
};

char const *matprod_errmsg(int error_code)
{
    return matprod_errmsg_[error_code > 0 ? error_code : -error_code];
}
