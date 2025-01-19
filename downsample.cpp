#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <stdint.h>

template <typename T>
T clamp(const T& value, const T& minValue, const T& maxValue) {
    return (value < minValue) ? minValue : (value > maxValue) ? maxValue : value;
}

// Simple low-pass filter (FIR filter)
std::vector<double> lowPassFilter(const std::vector<double>& input, int filterSize) {
    std::vector<double> output(input.size(), 0.0);
    double alpha = 1.0 / filterSize; // simple averaging filter

    for (size_t i = 0; i < input.size(); ++i) {
        for (int j = 0; j < filterSize; ++j) {
            if (i >= j) {
                output[i] += alpha * input[i - j];
            }
        }
    }
    return output;
}

// Downsampling function
std::vector<double> downsample(const std::vector<double>& input, int originalRate, int targetRate) {
    int decimationFactor = originalRate / targetRate;
    std::vector<double> filtered = lowPassFilter(input, 5);
    std::vector<double> output;

    for (size_t i = 0; i < filtered.size(); i += decimationFactor) {
        output.push_back(filtered[i]);
    }
    return output;
}

// Function to read multi-channel PCM data from a file
std::vector<double> readFromPCMFile(const std::string& filename, size_t sampleCount, int channels) {
    std::ifstream inFile(filename, std::ios::binary);
    if (!inFile) {
        std::cerr << "Error opening file for reading: " << filename << std::endl;
        return {};
    }

    std::vector<double> samples;
    samples.reserve(sampleCount * channels);

    // Read interleaved 16-bit PCM samples
    int16_t sample[2]; // Adjust for the number of channels
    while (inFile.read(reinterpret_cast<char*>(&sample), sizeof(sample))) {
        for (int i = 0; i < channels; ++i) {
            samples.push_back(static_cast<double>(sample[i]) / 32768.0); // Normalize to [-1.0, 1.0]
        }
    }

    inFile.close();
    return samples;
}

// Function to write multi-channel output to a PCM file
void writeToPCMFile(const std::string& filename, const std::vector<double>& samples, int channels) {
    std::ofstream outFile(filename, std::ios::binary);
    if (!outFile) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    size_t sampleCount = samples.size() / channels; // Total samples divided by number of channels
    for (size_t i = 0; i < sampleCount; ++i) {
        for (int ch = 0; ch < channels; ++ch) {
            // Clamp and convert to 16-bit PCM
            int16_t pcmSample = static_cast<int16_t>(clamp(samples[i * channels + ch] * 32767.0, -32768.0, 32767.0));
            outFile.write(reinterpret_cast<const char*>(&pcmSample), sizeof(pcmSample));
        }
    }

    outFile.close();
}

int main() {
    // Specify the input and output file paths
    std::string inputFilePath = "output_48khz_2.pcm"; // Replace with your input file path
    std::string outputFilePath = "output_8khz_2_v1.pcm"; // Output file path
    int originalRate = 48000; // Original sample rate
    int targetRate = 8000;    // Target sample rate
    int channels = 2;         // Number of audio channels (e.g., 2 for stereo)

    // Read PCM data from the input file
    std::vector<double> input = readFromPCMFile(inputFilePath, 44100, channels); // Adjust sample count as needed

    // Check if input is valid
    if (input.empty()) {
        return 1; // Exit if there was an error reading the file
    }

    // Downsample each channel separately
    std::vector<double> output;
    for (int ch = 0; ch < channels; ++ch) {
        // Extract channel samples
        std::vector<double> channelSamples;
        for (size_t i = ch; i < input.size(); i += channels) {
            channelSamples.push_back(input[i]);
        }
        // Downsample the channel
        std::vector<double> downsampledChannel = downsample(channelSamples, originalRate, targetRate);
        output.insert(output.end(), downsampledChannel.begin(), downsampledChannel.end());
    }

    // Write the downsampled audio data to a PCM file
    writeToPCMFile(outputFilePath, output, channels);

    std::cout << "Downsampling complete. Output written to " << outputFilePath << std::endl;

    return 0;
}