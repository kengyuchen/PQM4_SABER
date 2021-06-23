#include <stdio.h>
#include <stdint.h>
#include "SABER_params.h"
#include "NTT.h"

/*
	Function: Cooley–Tukey Butterfly
	
	Perform in-place Cooley–Tukey Butterfly transform on a
	polynomial with degree len.

	In-place means *in and *out can be the same.
*/
void CT_Butterfly(int64_t* out, const int64_t* in, const int8_t c_index, const int len){
	if (len <= 1){
		return;
	}
	int half = len >> 1;
	for (int i = 0; i < half; i++){
		int64_t a = in[i];
		int64_t b = in[i + half];
		out[i] = a + MODQ_MUL(b, root_table[c_index]);
		out[i + half] = a - MODQ_MUL(b, root_table[c_index]);
	}
	return;
}

/*
	Function: Gentleman–Sande Butterfly
	
	Perform in-place Gentleman–Sande Butterfly transform on a
	polynomial with degree len.

	In-place means *in and *out can be the same.
*/
void GS_Butterfly(int64_t* out, const int64_t* in, const int8_t c_index, const int len){
	if (len <= 1){
		return;
	}
	int half = len >> 1;
	for (int i = 0; i < half; i++){
		int64_t a = in[i];
		int64_t b = in[i + half];
		out[i] = a + b;
		out[i + half] = MODQ_MUL(a - b, inv_root_table[c_index]);
	}
	return;
}

void NTT_forward(int64_t *out, const int16_t *in){
	for (int i = 0; i < SABER_N; i++){
		out[i] = (int64_t)in[i];
	}

	for (int i = 1; i < 128; i++){
		if (i < 2){
			CT_Butterfly(out, out, bit_reverse[i], SABER_N);
		} else if (i < 4){
			CT_Butterfly(out + (SABER_N >> 1) * (i-2), out + (SABER_N >> 1) * (i-2), bit_reverse[i], SABER_N >> 1);
		} else if (i < 8){
			CT_Butterfly(out + (SABER_N >> 2) * (i-4), out + (SABER_N >> 2) * (i-4), bit_reverse[i], SABER_N >> 2);
		} else if (i < 16){
			CT_Butterfly(out + (SABER_N >> 3) * (i-8), out + (SABER_N >> 3) * (i-8), bit_reverse[i], SABER_N >> 3);
		} else if (i < 32){
			CT_Butterfly(out + (SABER_N >> 4) * (i-16), out + (SABER_N >> 4) * (i-16), bit_reverse[i], SABER_N >> 4);
		} else if (i < 64){
			CT_Butterfly(out + (SABER_N >> 5) * (i-32), out + (SABER_N >> 5) * (i-32), bit_reverse[i], SABER_N >> 5);
		} else if (i < 128) {
			CT_Butterfly(out + (SABER_N >> 6) * (i-64), out + (SABER_N >> 6) * (i-64), bit_reverse[i], SABER_N >> 6);
		}
	}
	
	return;
}

void NTT_inv(int16_t *out, const int64_t *in){
	int64_t out_64[SABER_N];
	for (int i = 0; i < SABER_N; i++){
		out_64[i] = in[i];
	}

	for (int i = 127; i > 0; i--){
		if (i < 2){
			GS_Butterfly(out_64, out_64, bit_reverse[i], SABER_N);
		} else if (i < 4){
			GS_Butterfly(out_64 + (SABER_N >> 1) * (i-2), out_64 + (SABER_N >> 1) * (i-2), bit_reverse[i], SABER_N >> 1);
		} else if (i < 8){
			GS_Butterfly(out_64 + (SABER_N >> 2) * (i-4), out_64 + (SABER_N >> 2) * (i-4), bit_reverse[i], SABER_N >> 2);
		} else if (i < 16){
			GS_Butterfly(out_64 + (SABER_N >> 3) * (i-8), out_64 + (SABER_N >> 3) * (i-8), bit_reverse[i], SABER_N >> 3);
		} else if (i < 32){
			GS_Butterfly(out_64 + (SABER_N >> 4) * (i-16), out_64 + (SABER_N >> 4) * (i-16), bit_reverse[i], SABER_N >> 4);
		} else if (i < 64){
			GS_Butterfly(out_64 + (SABER_N >> 5) * (i-32), out_64 + (SABER_N >> 5) * (i-32), bit_reverse[i], SABER_N >> 5);
		} else if (i < 128){
			GS_Butterfly(out_64 + (SABER_N >> 6) * (i-64), out_64 + (SABER_N >> 6) * (i-64), bit_reverse[i], SABER_N >> 6);
		} 
	}
	for (int i = 0; i < SABER_N; i++){
		out_64[i] = MODQ_MUL(out_64[i], inv_128);
	}

	for (int i = 0; i < SABER_N; i++){
		if (out_64[i] < 0){
			out_64[i] = out_64[i] + Q;
		}
		if (out_64[i] > (Q / 2)){
			out_64[i] = out_64[i] - Q;
		}
		out[i] = out_64[i] % 8192;
	}

	return;
}

/*
	Function: accumulated polynomial multiplication by NTT

	Perform polynomial multiplicaiton a*b by NTT trick and add
	the result to res.

	a and b are unsigned int, so it first turns a, b signed and
	change the value range from 0 ~ 2^SABER_EQ to 
	-2^(SABER_EQ-1) ~ 2^(SABER_EQ-1).
*/

void poly_mul_acc_NTT(const uint16_t a[SABER_N], const uint16_t b[SABER_N], uint16_t res[SABER_N])
{
    int16_t a_mod[SABER_N], b_mod[SABER_N], result[SABER_N];;
    int64_t at[SABER_N] = {0};
    int64_t bt[SABER_N] = {0};
    int64_t result_t[SABER_N] = {0};

    for (int i = 0; i < SABER_N; i++){
    	a_mod[i] = a[i] % 8192;
    	if (a_mod[i] > 4095){
    		a_mod[i] = a_mod[i] - 8192;
    	}
  		b_mod[i] = b[i] % 8192;
    	if (b_mod[i] > 4095){
    		b_mod[i] = b_mod[i] - 8192;
    	}
    }

    NTT_forward(at, a_mod);
    NTT_forward(bt, b_mod);

    for (int i = 0; i < SABER_N; i += 4){
        result_t[i] = (at[i]*bt[i] + MODQ_MUL(at[i+1], bt[i+1])*root_table[bit_reverse[64 + i/4]]) % Q;
        result_t[i+1] = (at[i]*bt[i+1] + at[i+1]*bt[i]) % Q;
        result_t[i+2] = (at[i+2]*bt[i+2] - MODQ_MUL(at[i+3], bt[i+3])*root_table[bit_reverse[64 + i/4]]) % Q;
        result_t[i+3] = (at[i+2]*bt[i+3] + at[i+3]*bt[i+2]) % Q;
    }

    NTT_inv(result, result_t);

    for (int i = 0; i < SABER_N; i++){
    	if (result[i] < 0){
    		result[i] = result[i] + 8192;
    	}
    }

    for (int i = 0; i < SABER_N; i++){
        res[i] += (uint16_t)result[i];
    }
}


/*
	Function: accumulated polynomial multiplication by NTT without inverse

	res += NTT(a) * NTT(b)
	Note that we perform the inverse of NTT in poly_mul_acc_NTT_inv
*/
void poly_mul_acc_NTT_no_inv(const uint16_t a[SABER_N], const uint16_t b[SABER_N], int64_t* result_t)
{
    int16_t a_mod[SABER_N], b_mod[SABER_N];
    int64_t at[SABER_N] = {0};
    int64_t bt[SABER_N] = {0};
    // int64_t result_t[SABER_N] = {0};

    for (int i = 0; i < SABER_N; i++){
    	a_mod[i] = a[i] % 8192;
    	if (a_mod[i] > 4095){
    		a_mod[i] = a_mod[i] - 8192;
    	}
  		b_mod[i] = b[i] % 8192;
    	if (b_mod[i] > 4095){
    		b_mod[i] = b_mod[i] - 8192;
    	}
    }

    NTT_forward(at, a_mod);
    NTT_forward(bt, b_mod);

    for (int i = 0; i < SABER_N; i += 4){
        result_t[i] = (result_t[i] + at[i]*bt[i] + at[i+1]* MODQ_MUL(bt[i+1], root_table[bit_reverse[64 + i/4]])) % Q;
        result_t[i+1] = (result_t[i+1] + (at[i]*bt[i+1] + at[i+1]*bt[i]) ) %Q ;
        result_t[i+2] = (result_t[i+2] + at[i+2]*bt[i+2] - at[i+3]* MODQ_MUL(bt[i+3], root_table[bit_reverse[64 + i/4]])) % Q;
        result_t[i+3] = (result_t[i+3] + (at[i+2]*bt[i+3] + at[i+3]*bt[i+2]) ) % Q;
    }
}

/*
	Function: inverse NTT

	res = NTT^-1(result_t)
*/
void poly_mul_acc_NTT_inv(const int64_t result_t[SABER_N], uint16_t* res){
	int16_t result[SABER_N];

    NTT_inv(result, result_t);

    for (int i = 0; i < SABER_N; i++){
    	if (result[i] < 0){
    		result[i] = result[i] + 8192;
    	}
    }
    for (int i = 0; i < SABER_N; i++){
        res[i] += (uint16_t)result[i];
    }
}

/***************************************
 * Incomplete NTT
 * The following are NTT with different layers.
***************************************/

void base_multiplication_layer(int64_t *result_t, int64_t *at, int64_t *bt, const int layer){
	if (layer == 7){
		for (int i = 0; i < SABER_N; i += 4){
			int32_t xi = root_table[bit_reverse[64 + i/4]];
	        result_t[i] = (result_t[i] + at[i]*bt[i] + at[i+1]* MODQ_MUL(bt[i+1], xi)) % Q;
	        result_t[i+1] = (result_t[i+1] + (at[i]*bt[i+1] + at[i+1]*bt[i]) ) %Q ;
	        result_t[i+2] = (result_t[i+2] + at[i+2]*bt[i+2] - at[i+3]* MODQ_MUL(bt[i+3], xi)) % Q;
	        result_t[i+3] = (result_t[i+3] + (at[i+2]*bt[i+3] + at[i+3]*bt[i+2]) ) % Q;
	    }
	} else if (layer == 6){
    	for (int i = 0; i < SABER_N; i += 8){
	    	int32_t xi = root_table[bit_reverse[32 + i/8]];
	    	int64_t result_low[8] = {0}, result_high[8] = {0};
	    	for (int j = 0; j < 4; j++){
		        for (int k = 0; k < 4; k++){
		            result_low[j+k] += at[i+j]*bt[i+k]%Q;
		        }
		    }
		    for (int j = 4; j < 8; j++){
		    	result_t[i+j-4] = (result_low[j-4] + (result_low[j] % Q) * xi) % Q;
		    }

		    for (int j = 0; j < 4; j++){
		        for (int k = 0; k < 4; k++){
		            result_high[j+k] += at[i+j+4]*bt[i+k+4]%Q;
		        }
		    }
		    for (int j = 4; j < 8; j++){
		    	result_t[i+j] = (result_high[j-4] - (result_high[j] % Q) * xi) % Q;
		    }
	    }
	} else if (layer == 5){
		for (int i = 0; i < SABER_N; i += 16){
			int32_t xi = root_table[bit_reverse[16 + i/16]];
			int64_t result_low[16] = {0}, result_high[16] = {0};
			for (int j = 0; j < 8; j++){
		        for (int k = 0; k < 8; k++){
		            result_low[j+k] += at[i+j]*bt[i+k];
		        }
		    }
		    for (int j = 8; j < 16; j++){
		    	result_t[i+j-8] += (result_low[j-8] + (result_low[j] % Q) * xi) % Q;
		    }

		    for (int j = 0; j < 8; j++){
		        for (int k = 0; k < 8; k++){
		            result_high[j+k] += at[i+j+8]*bt[i+k+8];
		        }
		    }
		    for (int j = 8; j < 16; j++){
		    	result_t[i+j] += (result_high[j-8] - (result_high[j] % Q) * xi) % Q;
		    }
		}
	}
}

void NTT_forward_layer(int64_t *out, const int16_t *in, const int layer){
	for (int i = 0; i < SABER_N; i++){
		out[i] = (int64_t)in[i];
	}

	for (int i = 1; i < (1 << layer); i++){
		if (i < 2){
			CT_Butterfly(out, out, bit_reverse[i], SABER_N);
		} else if (i < 4){
			CT_Butterfly(out + (SABER_N >> 1) * (i-2), out + (SABER_N >> 1) * (i-2), bit_reverse[i], SABER_N >> 1);
		} else if (i < 8){
			CT_Butterfly(out + (SABER_N >> 2) * (i-4), out + (SABER_N >> 2) * (i-4), bit_reverse[i], SABER_N >> 2);
		} else if (i < 16){
			CT_Butterfly(out + (SABER_N >> 3) * (i-8), out + (SABER_N >> 3) * (i-8), bit_reverse[i], SABER_N >> 3);
		} else if (i < 32){
			CT_Butterfly(out + (SABER_N >> 4) * (i-16), out + (SABER_N >> 4) * (i-16), bit_reverse[i], SABER_N >> 4);
		} else if (i < 64){
			CT_Butterfly(out + (SABER_N >> 5) * (i-32), out + (SABER_N >> 5) * (i-32), bit_reverse[i], SABER_N >> 5);
		} else if (i < 128) {
			CT_Butterfly(out + (SABER_N >> 6) * (i-64), out + (SABER_N >> 6) * (i-64), bit_reverse[i], SABER_N >> 6);
		}
	}
	
	return;
}

void NTT_inv_layer(int16_t *out, const int64_t *in, const int layer){
	int64_t out_64[SABER_N];
	for (int i = 0; i < SABER_N; i++){
		out_64[i] = in[i];
	}

	for (int i = (1 << layer) - 1; i > 0; i--){
		if (i < 2){
			GS_Butterfly(out_64, out_64, bit_reverse[i], SABER_N);
		} else if (i < 4){
			GS_Butterfly(out_64 + (SABER_N >> 1) * (i-2), out_64 + (SABER_N >> 1) * (i-2), bit_reverse[i], SABER_N >> 1);
		} else if (i < 8){
			GS_Butterfly(out_64 + (SABER_N >> 2) * (i-4), out_64 + (SABER_N >> 2) * (i-4), bit_reverse[i], SABER_N >> 2);
		} else if (i < 16){
			GS_Butterfly(out_64 + (SABER_N >> 3) * (i-8), out_64 + (SABER_N >> 3) * (i-8), bit_reverse[i], SABER_N >> 3);
		} else if (i < 32){
			GS_Butterfly(out_64 + (SABER_N >> 4) * (i-16), out_64 + (SABER_N >> 4) * (i-16), bit_reverse[i], SABER_N >> 4);
		} else if (i < 64){
			GS_Butterfly(out_64 + (SABER_N >> 5) * (i-32), out_64 + (SABER_N >> 5) * (i-32), bit_reverse[i], SABER_N >> 5);
		} else if (i < 128){
			GS_Butterfly(out_64 + (SABER_N >> 6) * (i-64), out_64 + (SABER_N >> 6) * (i-64), bit_reverse[i], SABER_N >> 6);
		} 
	}
	for (int i = 0; i < SABER_N; i++){
		if (layer == 7){
			out_64[i] = MODQ_MUL(out_64[i], inv_128);
		} else if (layer == 6){
			out_64[i] = MODQ_MUL(out_64[i], inv_64);
		} else if (layer == 5){
			out_64[i] = MODQ_MUL(out_64[i], inv_32);
		}
	}

	for (int i = 0; i < SABER_N; i++){
		if (out_64[i] < 0){
			out_64[i] = out_64[i] + Q;
		}
		if (out_64[i] > (Q / 2)){
			out_64[i] = out_64[i] - Q;
		}
		out[i] = out_64[i] % 8192;
	}

	return;
}

/*
	Function: accumulated polynomial multiplication by NTT

	Perform polynomial multiplicaiton a*b by NTT trick and add
	the result to res.

	a and b are unsigned int, so it first turns a, b signed and
	change the value range from 0 ~ 2^SABER_EQ to 
	-2^(SABER_EQ-1) ~ 2^(SABER_EQ-1).
*/

void poly_mul_acc_NTT_layer(const uint16_t a[SABER_N], const uint16_t b[SABER_N], uint16_t res[SABER_N], const int layer)
{
    int16_t a_mod[SABER_N], b_mod[SABER_N], result[SABER_N];;
    int64_t at[SABER_N] = {0};
    int64_t bt[SABER_N] = {0};
    int64_t result_t[SABER_N] = {0};

    for (int i = 0; i < SABER_N; i++){
    	a_mod[i] = a[i] % 8192;
    	if (a_mod[i] > 4095){
    		a_mod[i] = a_mod[i] - 8192;
    	}
  		b_mod[i] = b[i] % 8192;
    	if (b_mod[i] > 4095){
    		b_mod[i] = b_mod[i] - 8192;
    	}
    }

    NTT_forward_layer(at, a_mod, layer);
    NTT_forward_layer(bt, b_mod, layer);

    base_multiplication_layer(result_t, at, bt, layer);

    NTT_inv_layer(result, result_t, layer);

    for (int i = 0; i < SABER_N; i++){
    	if (result[i] < 0){
    		result[i] = result[i] + 8192;
    	}
    }

    for (int i = 0; i < SABER_N; i++){
        res[i] += (uint16_t)result[i];
    }
}


/*
	Function: accumulated polynomial multiplication by NTT without inverse

	res += NTT(a) * NTT(b)
	Note that we perform the inverse of NTT in poly_mul_acc_NTT_inv
*/
void poly_mul_acc_NTT_layer_no_inv(const uint16_t a[SABER_N], const uint16_t b[SABER_N], int64_t* result_t, const int layer)
{
    int16_t a_mod[SABER_N], b_mod[SABER_N];
    int64_t at[SABER_N] = {0};
    int64_t bt[SABER_N] = {0};
    //int64_t result_t[SABER_N] = {0};

    for (int i = 0; i < SABER_N; i++){
    	a_mod[i] = a[i] % 8192;
    	if (a_mod[i] > 4095){
    		a_mod[i] = a_mod[i] - 8192;
    	}
  		b_mod[i] = b[i] % 8192;
    	if (b_mod[i] > 4095){
    		b_mod[i] = b_mod[i] - 8192;
    	}
    }

    NTT_forward_layer(at, a_mod, layer);
    NTT_forward_layer(bt, b_mod, layer);

    base_multiplication_layer(result_t, at, bt, layer);
}

/*
	Function: inverse NTT

	res = NTT^-1(result_t)
*/
void poly_mul_acc_NTT_layer_inv(const int64_t result_t[SABER_N], uint16_t* res, const int layer){
	int16_t result[SABER_N];

    NTT_inv_layer(result, result_t, layer);

    for (int i = 0; i < SABER_N; i++){
    	if (result[i] < 0){
    		result[i] = result[i] + 8192;
    	}
    }
    for (int i = 0; i < SABER_N; i++){
        res[i] += (uint16_t)result[i];
    }
}


// Just for debug print
void print_NTT(const int64_t a[SABER_N], const char tag[64]){
	printf("%s\n", tag);
	for (int i = 0; i < SABER_N; i++){
		printf("%ld ", a[i]);
	}
	printf("\n");
}
void print_short(const int16_t a[SABER_N], const char tag[64]){
	printf("%s\n", tag);
	for (int i = 0; i < SABER_N; i++){
		printf("%d ", a[i]);
	}
	printf("\n");
}