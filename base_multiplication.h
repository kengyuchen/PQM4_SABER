#include "NTT.h"
void my_karatsuba(const int32_t *a, const int32_t *b, int64_t *result, const int len);
void base_multiplication_layer(int32_t *result_t, int32_t *at, int32_t *bt, const int layer);