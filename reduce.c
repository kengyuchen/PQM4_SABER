#include <stdio.h>
#include <stdint.h>
#include "reduce.h"

/*
    Function: Montgomery Reduction
    
    Perform Mongomery Reduction on Z_Q.
    Notice that montgomery reduce produce A/R mod Q,
    not A mod Q, so the result should be multiplied
    by R mod Q somewhere.
*/
int64_t montgomery_reduce(const int64_t a)
{
    if (a == 0) return 0;
    int32_t l;
    int64_t t;
    l = a & 0xffffffff;        
    l = (l * Q_INV) & 0xffffffff;
    t = (int64_t)l*Q;                     
    t = (a >> 32) - (t >> 32);                    
    return t;
}