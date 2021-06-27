#include <stdint.h>
#include "base_multiplication.h"
// Define whether to use Karatsuba algorithm
#define Karatsuba

/*
	Karatsuba algorithm: Improved Form
	
	Consider (a0+a1x)(b0+b1x) = (a0+a1)(b0+b1)x + (a0b0-a1b1x)(1-x)
	h1 = (a0+a1)(b0+b1)
	h2 = a0b0
	h3 = a1b1

*/
void my_karatsuba(const int32_t *a, const int32_t *b, int64_t *result, const int len){
	if (len == 2){
		result[0] = MUL64(a[0],b[0]);
		result[2] = MUL64(a[1],b[1]);
		result[1] = MUL64((a[0]+a[1]), (b[0]+b[1])) - result[0] - result[2];
		return;
	}
	 
	int64_t h1[len-1], h2[len-1], h3[len-1];
	int32_t f[len/2], g[len/2];

	for (int i = 0; i < len/2; i++){
		f[i] = a[i] + a[i + len/2]; // f = a0+a1
		g[i] = b[i] + b[i + len/2]; // g = b0+b1
	}

	my_karatsuba(f, g, h1, len/2); // h = (a0+a1)(b0+b1)
	my_karatsuba(a, b, h2, len/2);
	my_karatsuba(a + len/2, b + len/2, h3, len/2);

	// result = h2 + (h1-h2-h3)x^4 + h3x^8
	for (int i = 0; i < (len*2 - 1); i++){
		if (i < len/2){
			result[i] = h2[i];
		} else if (i < len - 1){
			result[i] = h2[i] + h1[i - len/2] - h2[i - len/2] - h3[i - len/2];
		} else if (i < len){
			result[i] = h1[i - len/2] - h2[i - len/2] - h3[i - len/2];
		} else if (i < len/2 + len - 1){
			result[i] = h1[i - len/2] - h2[i - len/2] - h3[i - len/2] + h3[i - len];
		} else{
			result[i] = h3[i - len];
		}
	}
}

void base_multiplication_layer(int32_t *result_t, int32_t *at, int32_t *bt, const int layer){
	if (layer == 7){
		for (int i = 0; i < SABER_N; i += 4){
			int32_t xi = root_table[bit_reverse[64 + i/4]];
#ifdef Karatsuba
			result_t[i] = (MUL64(at[i], bt[i]) + MUL64(at[i+1], MODQ_MUL(bt[i+1], xi))) % Q;
		    result_t[i+1] = (MUL64((at[i]+at[i+1]), (bt[i]+bt[i+1])) - MUL64(at[i], bt[i]) - MUL64(at[i+1], bt[i+1])) % Q;
			result_t[i+2] = (MUL64(at[i+2], bt[i+2]) - MUL64(at[i+3], MODQ_MUL(bt[i+3], xi))) % Q;
		    result_t[i+3] = (MUL64((at[i+2]+at[i+3]), (bt[i+2]+bt[i+3])) - MUL64(at[i+2], bt[i+2]) - MUL64(at[i+3], bt[i+3])) % Q;
#else
			result_t[i] = (result_t[i] + at[i]*bt[i] + at[i+1]* MODQ_MUL(bt[i+1], xi)) % Q;
	        result_t[i+1] = (result_t[i+1] + (at[i]*bt[i+1] + at[i+1]*bt[i]) ) % Q;
	        result_t[i+2] = (result_t[i+2] + at[i+2]*bt[i+2] - at[i+3]* MODQ_MUL(bt[i+3], xi)) % Q;
	        result_t[i+3] = (result_t[i+3] + (at[i+2]*bt[i+3] + at[i+3]*bt[i+2]) ) % Q;
#endif
	    }
	} else if (layer == 6){
    	for (int i = 0; i < SABER_N; i += 8){
	    	int32_t xi = root_table[bit_reverse[32 + i/8]];
	    	int64_t result_low[8] = {0}, result_high[8] = {0};
#ifdef Karatsuba
			my_karatsuba(at+i, bt+i, result_low, 4);
			for (int j = 4; j < 8; j++){
				result_t[i+j-4] = (result_low[j-4] + MODQ_MUL((result_low[j] % Q), xi)) % Q;
			}
			my_karatsuba(at+i+4, bt+i+4, result_high, 4);
			for (int j = 4; j < 8; j++){
				result_t[i+j] = (result_high[j-4] - MODQ_MUL((result_high[j] % Q) , xi)) % Q;
			}
#else
			for (int j = 0; j < 4; j++){
				for (int k = 0; k < 4; k++){
					result_low[j+k] += MODQ_MUL(at[i+j], bt[i+k]);
				}
			}
			for (int j = 4; j < 8; j++){
				result_t[i+j-4] = (result_low[j-4] + MODQ_MUL(result_low[j], xi)) % Q;
			}
			for (int j = 0; j < 4; j++){
				for (int k = 0; k < 4; k++){
					result_high[j+k] += MODQ_MUL(at[i+j+4], bt[i+k+4]);
				}
			}
			for (int j = 4; j < 8; j++){
				result_t[i+j] = (result_high[j-4] - MODQ_MUL(result_high[j], xi)) % Q;
			}
#endif
	    }
	} else if (layer == 5){
		for (int i = 0; i < SABER_N; i += 16){
			int32_t xi = root_table[bit_reverse[16 + i/16]];
			int64_t result_low[16] = {0}, result_high[16] = {0};
#ifdef Karatsuba
		    my_karatsuba(at+i, bt+i, result_low, 8);
		    for (int j = 8; j < 16; j++){
		    	result_t[i+j-8] += (result_low[j-8] + MODQ_MUL((result_low[j] % Q) , xi)) % Q;
		    }
		    my_karatsuba(at+i+8, bt+i+8, result_high, 8);
		    for (int j = 8; j < 16; j++){
		    	result_t[i+j] += (result_high[j-8] - MODQ_MUL((result_high[j] % Q), xi)) % Q;
		    }
#else
			for (int j = 0; j < 8; j++){
		        for (int k = 0; k < 8; k++){
		            result_low[j+k] += at[i+j]*bt[i+k];
		        }
		    }
		    for (int j = 8; j < 16; j++){
		    	result_t[i+j-8] += (result_low[j-8] + MODQ_MUL((result_low[j] % Q) , xi)) % Q;
		    }
		    
		    for (int j = 0; j < 8; j++){
		        for (int k = 0; k < 8; k++){
		            result_high[j+k] += at[i+j+8]*bt[i+k+8];
		        }
		    }
		    for (int j = 8; j < 16; j++){
		    	result_t[i+j] += (result_high[j-8] - MODQ_MUL((result_high[j] % Q), xi)) % Q;
		    }
#endif
		}
	}
}