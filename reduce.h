#include <stdio.h>
#include <stdint.h>
#include "SABER_params.h"

#define Q 25166081
#define Q_INV 4253089537 // Q_INV = Q^-1 mod R
#define R 16733526 // R = 2^32 mod Q

int64_t montgomery_reduce(const int64_t a);
