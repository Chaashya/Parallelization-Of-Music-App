using System;
using System.Numerics;
using System.Threading.Tasks;

namespace DigitalMusicAnalysis
{
    class FastFourierTransform
    {
        private static int BitReversal(int bits, int n)
        {
            int reversedN = n;
            int count = bits - 1;

            n >>= 1;
            while (n > 0)
            {
                reversedN = (reversedN << 1) | (n & 1);
                count--;
                n >>= 1;
            }

            return ((reversedN << count) & ((1 << bits) - 1));
        }

        public static Complex[] IterativeFFT(Complex[] x, int L, Complex[] twiddles)
        {
            int N = x.Length;
            Complex[] Y = new Complex[N];
            int bits = (int)Math.Log(N, 2);

            Parallel.For(0, N, new ParallelOptions { MaxDegreeOfParallelism = MainWindow.numThreads }, i =>
            {
                int pos = BitReversal(bits, i);
                Y[i] = x[pos];
            });

            for (int ii = 2; ii <= N; ii <<= 1)
            {
                for (int jj = 0; jj < N; jj += ii)
                {
                    for (int kk = 0; kk < ii / 2; kk++)
                    {
                        int e = jj + kk;
                        int o = jj + kk + (ii / 2);
                        Complex even = Y[e];
                        Complex odd = Y[o];

                        Y[e] = even + odd * twiddles[kk * (L / ii)];
                        Y[o] = even + odd * twiddles[(kk + (ii / 2)) * (L / ii)];
                    }
                }
            }

            return Y;
        }
    }
}
