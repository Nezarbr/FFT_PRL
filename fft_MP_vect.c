#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include <getopt.h>
#include <sys/time.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <sys/time.h>
#include <immintrin.h>

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint16_t u16;

double cutoff = 500;
u64 seed = 0;
u64 size = 0;
char *filename = "output.wav";

/******************** pseudo-random function (SPECK-like) ********************/

#define ROR(x, r) ((x >> r) | (x << (64 - r)))
#define ROL(x, r) ((x << r) | (x >> (64 - r)))
#define R(x, y, k) (x = ROR(x, 8), x += y, x ^= k, y = ROL(y, 3), y ^= x)
u64 PRF(u64 seed, u64 IV, u64 i)
{
        u64 y = i;
        u64 x = 0xBaadCafeDeadBeefULL;
        u64 b = IV;
        u64 a = seed;
        R(x, y, b);
        for (int i = 0; i < 32; i++) {
                R(a, b, i);
                R(x, y, b);
        }
        return x + i;
}

/************************** Fast Fourier Transform ***************************/
/* This code assumes that n is a power of two !!!                            */
/*****************************************************************************/

static inline __m256d complex_product(__m256d a, __m256d b) {
    __m256d ar = _mm256_unpacklo_pd(a, a);
    __m256d ai = _mm256_unpackhi_pd(a, a);

    __m256d br = _mm256_unpacklo_pd(b, b);
    __m256d bi = _mm256_unpackhi_pd(b, b);

    __m256d mul_r = _mm256_sub_pd(_mm256_mul_pd(ar, br), _mm256_mul_pd(ai, bi));
    __m256d mul_i = _mm256_add_pd(_mm256_mul_pd(ai, br), _mm256_mul_pd(ar, bi));

    return _mm256_shuffle_pd(mul_r, mul_i, 0x0C);
}



void FFT_rec(u64 n, const double complex *X, double complex *Y, u64 stride) {
    if (n == 1) {
        Y[0] = X[0];
        return;
    }

        #pragma omp task  if (n>1024)
        FFT_rec(n / 2, X, Y, 2 * stride);
        #pragma omp task  if (n>1024)
        FFT_rec(n / 2, X + stride, Y + n / 2, 2 * stride);

        #pragma omp taskwait



    double complex omega_n = cexp(-2.0 * I * M_PI / n);
    double complex omega = 1.0;

    __m256d omega_n_vect = _mm256_set_pd(cimag(omega_n), creal(omega_n), cimag(omega_n), creal(omega_n));
    __m256d omega_n_vect_sq = complex_product(omega_n_vect, omega_n_vect);


    __m256d omega_vect = _mm256_set_pd(cimag(omega_n), creal(omega_n), cimag(omega), creal(omega));

    

    for (u64 i = 0; i < n / 2; i += 2) {
        __m256d vect_odd = _mm256_loadu_pd((double *)&(Y[i]));
        __m256d vect_even = _mm256_loadu_pd((double *)&(Y[i + n / 2]));

        __m256d prod = complex_product(vect_even, omega_vect);
        vect_even = prod;

        __m256d res_add = _mm256_add_pd(vect_odd, vect_even);
        __m256d res_sub = _mm256_sub_pd(vect_odd, vect_even);

        _mm256_storeu_pd((double *)&(Y[i]), res_add);
        _mm256_storeu_pd((double *)&(Y[i + n / 2]), res_sub);

        prod = complex_product(omega_vect, omega_n_vect_sq);
        omega_vect = prod;
    }
}


void FFT(u64 n, const double complex * X, double complex *Y)
{
        /* sanity check */
        if ((n & (n - 1)) != 0) 
                errx(1, "size is not a power of two (this code does not handle other cases)");
        #pragma omp parallel
        #pragma omp single
        FFT_rec(n, X, Y, 1);                        /* stride == 1 initially */
}

/* Computes the inverse Fourier transform, but destroys the input */
void iFFT(u64 n, double complex * X, double complex *Y)
{
        for (u64 i = 0; i < n; i++)
                X[i] = conj(X[i]);

        FFT(n, X, Y);

        for (u64 i = 0; i < n; i++)
                Y[i] = conj(Y[i]) / n;
}

/******************* utility functions ********************/

double wtime()
{
        struct timeval ts;
        gettimeofday(&ts, NULL);
        return (double) ts.tv_sec + ts.tv_usec / 1e6;
}

void process_command_line_options(int argc, char **argv)
{
        struct option longopts[5] = {
                {"size", required_argument, NULL, 'n'},
                {"seed", required_argument, NULL, 's'},
                {"output", required_argument, NULL, 'o'},
                {"cutoff", required_argument, NULL, 'c'},
                {NULL, 0, NULL, 0}
        };
        char ch;
        while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
                switch (ch) {
                case 'n':
                        size = atoll(optarg);
                        break;
                case 's':
                        seed = atoll(optarg);
                        break;
                case 'o':
                        filename = optarg;
                        break;
                case 'c':
                        cutoff = atof(optarg);
                        break;
                default:
                        errx(1, "Unknown option\n");
                }
        }
        /* validation */
        if (size == 0)
                errx(1, "missing --size argument");
}

/* save at most 10s of sound output in .WAV format */
void save_WAV(char *filename, u64 size, double complex *C)
{
        assert(size < 1000000000);
        FILE *f = fopen(filename, "w");
        if (f == NULL)
                err(1, "fopen");
        printf("Writing <= 10s of audio output in %s\n", filename);
        u32 rate = 44100;        // Sample rate
        u32 frame_count = 10*rate;
        if (size < frame_count)
                frame_count = size;
        u16 chan_num = 2;        // Number of channels
        u16 bits = 16;           // Bit depth
        u32 length = frame_count*chan_num*bits / 8;
        u16 byte;
        double multiplier = 32767;

        /* WAVE Header Data */
        fwrite("RIFF", 1, 4, f);
        u32 chunk_size = length + 44 - 8;
        fwrite(&chunk_size, 4, 1, f);
        fwrite("WAVE", 1, 4, f);
        fwrite("fmt ", 1, 4, f);
        u32 subchunk1_size = 16;
        fwrite(&subchunk1_size, 4, 1, f);
        u16 fmt_type = 1;  // 1 = PCM
        fwrite(&fmt_type, 2, 1, f);
        fwrite(&chan_num, 2, 1, f);
        fwrite(&rate, 4, 1, f);
        // (Sample Rate * BitsPerSample * Channels) / 8
        uint32_t byte_rate = rate * bits * chan_num / 8;
        fwrite(&byte_rate, 4, 1, f);
        uint16_t block_align = chan_num * bits / 8;
        fwrite(&block_align, 2, 1, f);
        fwrite(&bits, 2, 1, f);

        /* Marks the start of the data */
        fwrite("data", 1, 4, f);
        fwrite(&length, 4, 1, f);  // Data size
        for (u32 i = 0; i < frame_count; i++)
        {
                byte = creal(C[i]) * multiplier;
                fwrite(&byte, 2, 1, f);
                byte = cimag(C[i]) * multiplier;
                fwrite(&byte, 2, 1, f);
        }
        fclose(f);
}



double wallclock_time()
{
	struct timeval tmp_time;
	gettimeofday(&tmp_time, NULL);
	return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6);
}

/*************************** main function *********************************/

int main(int argc, char **argv)
{


    process_command_line_options(argc, argv);


    /* starting timer */
	double start = wallclock_time();
        
        /* generate white noise */
    
    double complex *B = malloc(size * sizeof(*B));
    double complex *C = malloc(size * sizeof(*C));
    double complex *A = malloc(size * sizeof(*A));

    printf("Generating white noise...\n");

    #pragma omp parallel for schedule(auto)
    for (u64 i = 0; i < size; i++) {
            double real = 2 * (PRF(seed, 0, i) * 5.42101086242752217e-20) - 1;
            double imag = 2 * (PRF(seed, 1, i) * 5.42101086242752217e-20) - 1;
            A[i] = real + imag * I;
        }

    printf("Forward FFT...\n");
    FFT(size, A, B);

        /* damp fourrier coefficients */
    printf("Adjusting Fourier coefficients...\n");

    #pragma omp parallel for schedule(auto)
    for (u64 i = 0; i < size; i++) {
        double tmp = sin(i * 2 * M_PI / 44100);
        B[i] *= tmp * cexp(-i*2*I*M_PI / 4 / 44100);
        B[i] *= (i+1) / exp((i * cutoff) / size);
    }
        
    printf("Inverse FFT...\n");
    iFFT(size, B, C);

    printf("Normalizing output...\n");
    double max = 0;


    #pragma omp parallel for reduction(max:max)
    for (u64 i = 0; i < size; i++)
        max = fmax(max, cabs(C[i]));

    printf("max = %g\n", max);

    #pragma omp parallel for schedule(auto)
    for (u64 i = 0; i < size; i++)
            C[i] /= max;
    

    double end = wallclock_time();
    fprintf(stderr, "Total computing time: %g sec\n", end - start);

    if (filename != NULL)
        save_WAV(filename, size, C);
        
    exit(EXIT_SUCCESS);
}