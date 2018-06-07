#define _FFT_H_

#include "complex_simple.h"

complex* DFT_naive_1(complex* x, int N);
complex* DFT_naive_2(complex* x, int N);
complex* FFT_simple(complex* x, int N /* must be a power of 2 */);

#endif /* #ifndef _FFT_H_ */

