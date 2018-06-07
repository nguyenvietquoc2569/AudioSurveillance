#ifndef FFTSSL_H
#define FFTSSL_H

class FFTSSL
{
public:
    FFTSSL(int size);
    ~FFTSSL();
    static int log2(int x);
    int getSize();
    double* compute(double* data,int datasize);

private:
    int mSize;
    int mLogSize;
    double* mLookUpTable;
};

#endif // FFT_H
