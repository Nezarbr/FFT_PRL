#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define SAMPLE_RATE 44100
#define DURATION 10 // seconds
#define FILENAME "white_noise.wav"

typedef uint64_t u64;

// Pseudo-random number generator function
u64 PRF(u64 seed, u64 i) {
    // A simple pseudo-random number generator
    // (In practice, use a better PRG for higher quality noise)
    return seed + i * 6364136223846793005ULL + 1442695040888963407ULL;
}

// Function to write WAV file header
void write_wav_header(FILE *f, int sampleRate, int totalSamples) {
    int byteRate = sampleRate * 2 * 2; // 16-bit mono
    int blockSize = 2 * 2; // 16-bit mono
    int dataSize = totalSamples * 2 * 2; // 16-bit mono, total number of samples
    fwrite("RIFF", 4, 1, f); // ChunkID
    int chunkSize = 36 + dataSize;
    fwrite(&chunkSize, 4, 1, f); // ChunkSize
    fwrite("WAVE", 4, 1, f); // Format
    fwrite("fmt ", 4, 1, f); // Subchunk1ID
    int subChunk1Size = 16; // for PCM
    fwrite(&subChunk1Size, 4, 1, f); // Subchunk1Size
    short audioFormat = 1; // PCM = 1
    fwrite(&audioFormat, 2, 1, f); // AudioFormat
    short numChannels = 2; // Mono = 1, Stereo = 2
    fwrite(&numChannels, 2, 1, f); // NumChannels
    fwrite(&sampleRate, 4, 1, f); // SampleRate
    fwrite(&byteRate, 4, 1, f); // ByteRate
    fwrite(&blockSize, 2, 1, f); // BlockAlign
    short bitsPerSample = 16;
    fwrite(&bitsPerSample, 2, 1, f); // BitsPerSample
    fwrite("data", 4, 1, f); // Subchunk2ID
    fwrite(&dataSize, 4, 1, f); // Subchunk2Size
}

int main() {
    FILE *f = fopen(FILENAME, "wb");
    if (!f) {
        perror("Error opening file");
        return EXIT_FAILURE;
    }
    
    int totalSamples = SAMPLE_RATE * DURATION;
    write_wav_header(f, SAMPLE_RATE, totalSamples);

    u64 seed = time(NULL); // Use current time as seed
    for (int i = 0; i < totalSamples; i++) {
        double sample = (double)PRF(seed, i) / (double)UINT64_MAX;
        sample = sample * 2.0 - 1.0; // Normalize to [-1,1]
        short intSample = (short)(sample * 32767.0); // Convert to 16-bit sample
        fwrite(&intSample, sizeof(short), 1, f); // Left channel
        fwrite(&intSample, sizeof(short), 1, f); // Right channel
    }

    fclose(f);
    printf("White noise has been generated and saved to %s\n", FILENAME);

    return EXIT_SUCCESS;
}
