#include <stdio.h>
#include <stdint.h>
#include "SABER_params.h"
#include "NTT.h"
#include "reduce.h" 

/*************************************************
* Name:        montgomery_reduce
*
* Description: Montgomery reduction; given a 32-bit integer a, computes
*              16-bit integer congruent to a * R^-1 mod q,
*              where R=2^16
*
* Arguments:   - int32_t a: input integer to be reduced;
*                           has to be in {-q2^15,...,q2^15-1}
*
* Returns:     integer in {-q+1,...,q-1} congruent to a * R^-1 modulo q.
**************************************************/
int64_t montgomery_reduce(const int64_t a)
{
  int64_t l, t;
  l = a& 0xffffffff;
  l = l * (Q_INV);               // l = a * M' 
  l = l& 0xffffffff;             // l = l mod R
  t = (int64_t)l*Q;              // t = a * M'*M 
  t = t& 0xffffffff;             // t = t mod R
  t = a - t;                     // a - lm
  t >>= 32;                      // (a - lm) / R
  if (t > Q) t = t-Q;
  if (t < -Q) t = t+Q;
  return t;
}

// int64_t normal_reduce(const int64_t a){
//     return (int64_t)(a%Q);
// }