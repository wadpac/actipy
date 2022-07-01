#include <sstream>
#include <string>
#include <limits>
#include <fstream>
#include <vector>
#include <iomanip>  // get_time

#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

// geplakt voorbeeld van https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html#vectorizing-functions om even uit te proberen
py::array_t<double> add_arrays(py::array_t<double> input1, py::array_t<double> input2) {
    py::buffer_info buf1 = input1.request(), buf2 = input2.request();

    if (buf1.ndim != 1 || buf2.ndim != 1)
        throw std::runtime_error("Number of dimensions must be one");

    if (buf1.size != buf2.size)
        throw std::runtime_error("Input shapes must match");

    /* No pointer is passed, so NumPy will allocate the buffer */
    auto result = py::array_t<double>(buf1.size);

    py::buffer_info buf3 = result.request();

    double *ptr1 = static_cast<double *>(buf1.ptr);
    double *ptr2 = static_cast<double *>(buf2.ptr);
    double *ptr3 = static_cast<double *>(buf3.ptr);

    for (long idx = 0; idx < buf1.shape[0]; idx++)
        ptr3[idx] = ptr1[idx] + ptr2[idx];

    return result;

    // // ander voorbeeld, hoe je snel (zonder bounds checks) itereert over dimensies:
    // auto r = x.unchecked<3>(); // x must have ndim = 3; can be non-writeable
    // double sum = 0;
    // for (py::ssize_t i = 0; i < r.shape(0); i++)
    //     for (py::ssize_t j = 0; j < r.shape(1); j++)
    //         for (py::ssize_t k = 0; k < r.shape(2); k++)
    //             sum += r(i, j, k);
    // return sum;
}

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


std::tuple<py::dict, py::array_t<long>, py::array_t<double>, py::array_t<double>, py::array_t<double>, py::array_t<double>> readFile(std::string accFile, bool verbose) {

    py::dict info;

    int fileHeaderSize = 59;
    int linesToAxesCalibration = 47;
    int blockHeaderSize = 9;
    int statusOK = -1;
    double sampleRate = -1;
    int errCounter = 0;

    // NpyWriter writer = new NpyWriter(outFile, ITEM_NAMES_AND_TYPES);
    // py::array_t<long> time_array;
    // py::array_t<double> x_array, y_array, z_array, T_array;
    std::vector<long> time_array;
    std::vector<double> x_array, y_array, z_array, T_array;

    auto max_streamsize = std::numeric_limits<std::streamsize>::max();

    try {
        std::ifstream input_file(accFile);
        // Read header to determine mfrGain and mfrOffset values
        double mfrGain[3];
        int mfrOffset[3];
        int numBlocksTotal = parseBinFileHeader(input_file, fileHeaderSize, linesToAxesCalibration, mfrGain, mfrOffset);

        std::cout << "WARNING: Remove iostream dependency before publishing!\n";

        std::cout << "numBlocksTotal: " << numBlocksTotal << std::endl;
        std::cout << "mfrGain[0]: " << mfrGain[0] << std::endl;
        std::cout << "mfrOffset[0]: " << mfrOffset[0] << std::endl;
        std::cout << "mfrGain[1]: " << mfrGain[1] << std::endl;
        std::cout << "mfrOffset[1]: " << mfrOffset[1] << std::endl;
        std::cout << "mfrGain[2]: " << mfrGain[2] << std::endl;
        std::cout << "mfrOffset[2]: " << mfrOffset[2] << std::endl;

        int blockCount = 0;
        std::string header;
        long blockTime = 0;  // Unix millis
        double temperature = 0.0;
        double freq = 0.0;
        std::string data;
        std::string timeFmtStr = "Page Time:%Y-%m-%d %H:%M:%S:";

        std::string line;
        while (std::getline(input_file, line)) {
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
                        
                        // Note: not doing anything with time zone here, not sure if necessary, but
                        // Java version specifies UTC. I'm not sure whether that does anything, because
                        // I'd say it assumes that the data is in UTC, so if that assumption breaks,
                        // any specified timezone here won't fix that, since no timezone information
                        // is read from the file anywhere anyway.
                        blockTime = std::chrono::duration_cast<std::chrono::milliseconds>(tp.time_since_epoch()).count() + milliseconds;
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
                    std::cerr << "header error: ";
                    std::cerr << e.what() << '\n';
                    continue;
                }
            }
            sampleRate = freq;

            // now process hex data
            std::getline(input_file, data);

            // raw reading values
            int hexPosition = 0;
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
                    // todo *** read in light[36:46] (10 bits to signed int) and
                    // button[47] (bool) values...

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
                    std::cerr << "data error at i = " << i << ": ";
                    std::cerr << e.what() << '\n';
                    break;  // rest of this block could be corrupted
                }
            }

            // Progress bar
            blockCount++;
            if (verbose) {
                if ((blockCount % 10000 == 0) || (blockCount == numBlocksTotal)) {
                    printf("Reading file... %d%%\r", (blockCount * 100 / numBlocksTotal));
                }
            }

        }
        statusOK = 1;
    } catch (const std::exception &e) {
        std::cerr << "an error occurred while reading!\n" << e.what();
        statusOK = 0;
    }

    info["ReadOK"] = statusOK;
    info["ReadErrors"] = errCounter;
    info["SampleRate"] = sampleRate;

    // return {std::move(info), std::move(time_array), std::move(x_array),
    //         std::move(y_array), std::move(z_array), std::move(T_array)};
    return {info, py::cast(time_array), py::cast(x_array),
            py::cast(y_array), py::cast(z_array), py::cast(T_array)};
}

PYBIND11_MODULE(GENEActivReaderCPP, m) {
    m.def("add_arrays", &add_arrays, "Add two NumPy arrays");
    // N.B.: don't use 'read' as the C++ function name, it is already used by
    // some include and will not compile (so we use readFile):
    m.def("read", &readFile, "Read a GENEActive file");
}
