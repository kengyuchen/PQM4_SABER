/*
	Karatsuba algorithm: Improved Form
	
	Consider (a0+a1x)(b0+b1x) = (a0+a1)(b0+b1)x + (a0b0-a1b1x)(1-x)
	h1 = (a0+a1)(b0+b1)
	h2 = a0b0
	h3 = a1b1

*/
void my_karatsuba(const int64_t *a, const int64_t *b, int64_t *result, const int len){
	if (len == 2){
		result[0] = a[0]*b[0];
		result[2] = a[1]*b[1];
		result[1] = (a[0]+a[1])*(b[0]+b[1]) - result[0] - result[2];
		return;
	}
	 
	int64_t h1[len-1], h2[len-1], h3[len-1];
	int64_t f[len/2], g[len/2];

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