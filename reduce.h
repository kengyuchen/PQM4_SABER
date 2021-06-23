#include <stdio.h>
#include <stdint.h>
#include "SABER_params.h"

#define Q_INV 4253089537

// int64_t barrett_reduce(int64_t a);
int64_t montgomery_reduce(const int64_t a);
// int64_t normal_reduce(const int64_t a);
