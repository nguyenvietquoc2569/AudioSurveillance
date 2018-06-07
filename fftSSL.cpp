
#include "fftSSL.h"
#include <QtGlobal>
#include <qmath.h>
const double PI  =3.141592653589793238462;

FFTSSL::FFTSSL(int size)
{
    mLogSize = FFTSSL::log2(size);
    Q_ASSERT((1 << mLogSize) == size);
    mSize = size;
    int le = size >> 1;
    mLookUpTable = (double*) malloc(sizeof(double)*2 * (le - 1));
    double theta = PI / le;
    double wrecur[2] = { qCos(theta), -qSin(theta) };
    double w[2] = { wrecur[0], wrecur[1] };
    for (int i = 0; i < (le - 1); ++i) {
        mLookUpTable[2 * i] = (double) wrecur[0];
        mLookUpTable[2 * i + 1] = (double) wrecur[1];
        wrecur[0] = wrecur[0] * w[0] - wrecur[1] * w[1];
        wrecur[1] = wrecur[1] * w[0] + wrecur[0] * w[1];
    }
}

FFTSSL::~FFTSSL()
{
    free(mLookUpTable);
}

int FFTSSL::log2(int x)
{
    Q_ASSERT(x>0);
    x--;
    int mask = 1;
    int i=0;
    while (x!=0)
    {
        x=x & (~mask);
        mask = mask << 1;
        i++;
    }
    return i;
}

int FFTSSL::getSize()
{
    return mSize;
}

double *FFTSSL::compute(double *data,int datasize)
{
    Q_ASSERT(datasize==mSize);
    double * result = (double*)malloc(sizeof(double)*mSize*2);
    for(int i=0; i < datasize ; i++)
    {
        result[i<<1] = data[i];
        result[(i<<1) +1 ] = 0.0f;
    }

    for (int l = 0, le = mSize, windex = 1; l < mLogSize; ++l, windex <<= 1) {
                le >>= 1;
                // first iteration with no multiplies
                for (int i = 0; i < mSize; i += (le << 1)) {
                    double temp[2] = { result[i << 1] + result[(i + le) << 1],
                            result[(i << 1) + 1] + result[((i + le) << 1) + 1] };
                    result[(i + le) << 1] = result[i << 1] - result[(i + le) << 1];
                    result[((i + le) << 1) + 1] = result[(i << 1) + 1]
                            - result[((i + le) << 1) + 1];
                    result[i << 1] = temp[0];
                    result[(i << 1) + 1] = temp[1];
                }
                // remaining iterations use stored mLookUpTable
                int uindex = windex - 1;
                for (int j = 1; j < le; ++j, uindex += windex) {
                    double u[2] = { mLookUpTable[uindex << 1],
                            mLookUpTable[(uindex << 1) + 1] };
                    for (int i = j; i < mSize; i += (le << 1)) {
                        double temp[2] = { result[i << 1] + result[(i + le) << 1],
                                result[(i << 1) + 1] + result[((i + le) << 1) + 1] };
                        double tm[2] = { result[i << 1] - result[(i + le) << 1],
                                result[(i << 1) + 1] - result[((i + le) << 1) + 1] };
                        result[(i + le) << 1] = tm[0] * u[0] - tm[1] * u[1];
                        result[((i + le) << 1) + 1] = tm[1] * u[0] + tm[0] * u[1];
                        result[i << 1] = temp[0];
                        result[(i << 1) + 1] = temp[1];
                    }
                }
            }

    for (int i = 1, j = 0; i < (mSize - 1); ++i) {
                int k = mSize >> 1;
                while (k <= j) {
                    j -= k;
                    k >>= 1;
                }
                j += k;
                if (i < j) { // swap result[i] and result[j]
                    double temp[2] = { result[i << 1], result[(i << 1) + 1] };
                    result[i << 1] = result[j << 1];
                    result[(i << 1) + 1] = result[(j << 1) + 1];
                    result[j << 1] = temp[0];
                    result[(j << 1) + 1] = temp[1];
                }
            }
            return result;

}
