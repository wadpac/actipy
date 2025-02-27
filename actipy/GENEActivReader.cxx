#include <sstream>
#include <string>
#include <limits>
#include <fstream>
#include <vector>
#include <iomanip>  // get_time

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

/**
 * Replicates bin file header, also calculates and returns
 * x/y/z gain/offset values along with number of pages of data in file bin
 * format described in GENEActiv manual ("Decoding .bin files", pg.27)
 * http://www.geneactiv.org/wp-content/uploads/2014/03/
 * geneactiv_instruction_manual_v1.2.pdf
 */
int parseBinFileHeader(std::istream& input_file, int fileHeaderSize, int linesToAxesCalibration,
                        double (&gainVals)[3], int (&offsetVals)[3]) {
    // skip the first i lines in the file
    auto max_streamsize = std::numeric_limits<std::streamsize>::max();
    for (int i = 0; i < linesToAxesCalibration; i++) {
        input_file.ignore(max_streamsize, '\n');
    }
    // read axes calibration lines for gain and offset values
    // data like -> x gain:25548 \n x offset:574 ... Volts:300 \n Lux:800
    std::string line;
    input_file.ignore(max_streamsize, ':');
    input_file >> gainVals[0]; // xGain
    input_file.ignore(max_streamsize, ':');
    input_file >> offsetVals[0]; // xOffset
    input_file.ignore(max_streamsize, ':');
    input_file >> gainVals[1]; // y
    input_file.ignore(max_streamsize, ':');
    input_file >> offsetVals[1]; // y
    input_file.ignore(max_streamsize, ':');
    input_file >> gainVals[2]; // z
    input_file.ignore(max_streamsize, ':');
    input_file >> offsetVals[2]; // z

    int volts, lux, numBlocksTotal;

    input_file.ignore(max_streamsize, ':');
    input_file >> volts; // volts
    input_file.ignore(max_streamsize, ':');
    input_file >> lux; // lux

    input_file.ignore(max_streamsize, '\n'); // 9 blank
    input_file.ignore(max_streamsize, '\n'); // 10 memory status header

    input_file.ignore(max_streamsize, ':');
    input_file >> numBlocksTotal; // 11
    input_file.ignore(max_streamsize, '\n');

    // ignore remaining header lines in bin file
    for (int i = 0; i < fileHeaderSize - linesToAxesCalibration - 11; i++) {
        input_file.ignore(max_streamsize, '\n');
    }
    return numBlocksTotal;
}


int getSignedIntFromHex(const std::string &hex) {
    // input hex base is 16
    int rawVal = std::stoll(hex, nullptr, 16);
    int unsignedLimit = 4096; // 2^[length*4] #i.e. 3 hexBytes (12 bits)
                                // limit = 4096
    int signedLimit = 2048; // 2^[length*(4-1)] #i.e. 3 hexBytes - 1 bit (11
                            // bits) limit = 2048
    if (rawVal > signedLimit) {
        rawVal = rawVal - unsignedLimit;
    }
    return rawVal;
}


std::tuple<py::dict, py::array_t<long>, py::array_t<float>, py::array_t<float>, py::array_t<float>, py::array_t<float>> readFile(std::string accFile, bool verbose, std::size_t start = 0, std::size_t end = 0) {

    py::dict info;

    int fileHeaderSize = 59;
    int linesToAxesCalibration = 47;
    int blockHeaderSize = 9;
    int statusOK = -1;
    double sampleRate = -1;
    int errCounter = 0;

    std::vector<long> time_array;
    std::vector<float> x_array, y_array, z_array, T_array;

    auto max_streamsize = std::numeric_limits<std::streamsize>::max();

    try {
        std::ifstream input_file(accFile);
        // Read header to determine mfrGain and mfrOffset values
        double mfrGain[3];
        int mfrOffset[3];
        int numBlocksTotalint = parseBinFileHeader(input_file, fileHeaderSize, linesToAxesCalibration, mfrGain, mfrOffset);
        if (numBlocksTotalint < 0) {
            fprintf(stderr, "WARNING: numBlocksTotal read in from header is negative, file corrupted?\n");
        }
        std::size_t numBlocksTotal = numBlocksTotalint;

        std::size_t blockCount = 0;
        std::string header;
        long blockTime = 0;  // Unix millis
        double temperature = 0.0;
        double freq = 0.0;
        std::string data;
        std::string timeFmtStr = "Page Time:%Y-%m-%d %H:%M:%S:";

        if (start == 0) {
            start = 1;
        }
        if (end == 0) {
            end = numBlocksTotal;
        }

        std::string line;
        while (std::getline(input_file, line)) {
            ++blockCount;
            if (blockCount >= start && blockCount <= end) {
                // header: "Recorded Data" (0), serialCode (1), seq num (2),
                // blockTime (3), unassigned (4), temp (5), batteryVolt (6),
                // deviceStatus (7), freq (8), data (9)
                for (int i = 1; i < blockHeaderSize; i++) {
                    try {
                        std::getline(input_file, header);
                        if (i == 3) {
                            std::tm tm = {};
                            std::stringstream ss(header);
                            ss >> std::get_time(&tm, timeFmtStr.c_str());
                            int milliseconds;
                            ss >> milliseconds;
                            auto tp = std::chrono::system_clock::from_time_t(std::mktime(&tm));
                            // Note: the timezone variable may not be portable to all OS'es!
                            blockTime = std::chrono::duration_cast<std::chrono::milliseconds>(tp.time_since_epoch()).count() + milliseconds - timezone * 1000;

                            // The above could be replaced by the following OS-portable C++20 when
                            // all compilers support it:
                            // std::chrono::utc_time<std::chrono::seconds> tp;
                            // std::stringstream ss(header);
                            // ss >> std::chrono::parse(timeFmtStr, tp);
                            // int milliseconds;
                            // ss >> milliseconds;
                            // blockTime = std::chrono::duration_cast<std::chrono::milliseconds>(tp.time_since_epoch()).count() + milliseconds;
                        } else if (i == 5) {
                            std::stringstream ss(header);
                            ss.ignore(max_streamsize, ':');
                            ss >> temperature;
                        } else if (i == 8) {
                            std::stringstream ss(header);
                            ss.ignore(max_streamsize, ':');
                            ss >> freq;
                        }
                    } catch (const std::exception &e) {
                        errCounter++;
                        fprintf(stderr, "header error: %s\n", e.what());
                        continue;
                    }
                }
                sampleRate = freq;

                // now process hex data
                std::getline(input_file, data);

                // raw reading values
                std::size_t hexPosition = 0;
                int xRaw = 0;
                int yRaw = 0;
                int zRaw = 0;
                double x = 0.0;
                double y = 0.0;
                double z = 0.0;
                double t = 0.0;

                int i = 0;
                while (hexPosition < data.size() - 1) {
                    try {
                        std::stringstream ss;
                        xRaw = getSignedIntFromHex(data.substr(hexPosition, 3));
                        yRaw = getSignedIntFromHex(data.substr(hexPosition + 3, 3));
                        zRaw = getSignedIntFromHex(data.substr(hexPosition + 6, 3));

                        // Update values to calibrated measure (taken from GENEActiv manual)
                        x = (xRaw * 100. - mfrOffset[0]) / mfrGain[0];
                        y = (yRaw * 100. - mfrOffset[1]) / mfrGain[1];
                        z = (zRaw * 100. - mfrOffset[2]) / mfrGain[2];

                        t = (double)blockTime + (double)i * (1.0 / freq) * 1000;  // Unix millis

                        time_array.push_back(t);
                        x_array.push_back(x);
                        y_array.push_back(y);
                        z_array.push_back(z);
                        T_array.push_back(temperature);

                        hexPosition += 12;
                        i++;
                    } catch (const std::exception &e) {
                        errCounter++;
                        fprintf(stderr, "data error at i = %d: %s\n", i, e.what());
                        break;  // rest of this block could be corrupted
                    }
                }
            } else if (blockCount < start) {
                // skip this block
                for (int i = 1; i < blockHeaderSize; i++) {  // header
                    input_file.ignore(max_streamsize, '\n');
                }
                input_file.ignore(max_streamsize, '\n');     // hexdata
            } else {
                // after end, no need to scan further
                break;
            }

            // Progress bar
            if (verbose) {
                if ((blockCount % 10000 == 0) || (blockCount == numBlocksTotal)) {
                    printf("Reading file... %lu%%\r", (blockCount * 100 / numBlocksTotal));
                }
            }

        }
        statusOK = 1;
        info["numBlocksTotal"] = numBlocksTotal;
    } catch (const std::exception &e) {
        fprintf(stderr, "an error occurred while reading!\n%s\n", e.what());
        statusOK = 0;
    }

    info["ReadOK"] = statusOK;
    info["ReadErrors"] = errCounter;
    info["SampleRate"] = sampleRate;

    return {info, py::cast(time_array), py::cast(x_array),
            py::cast(y_array), py::cast(z_array), py::cast(T_array)};
}

PYBIND11_MODULE(GENEActivReaderCPP, m) {
    // N.B.: don't use 'read' as the C++ function name, it is already used by
    // some include and will not compile (so we use readFile):
    m.def("read", &readFile, "Read a GENEActive file",
          // declaring arguments to be able to declare the default values for start and end
          py::arg("accFile"), py::arg("verbose"), py::arg("start") = 0, py::arg("end") = 0);
}
