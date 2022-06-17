#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

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

    for (size_t idx = 0; idx < buf1.shape[0]; idx++)
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

PYBIND11_MODULE(GENEActivReaderCPP, m) {
    m.def("add_arrays", &add_arrays, "Add two NumPy arrays");
}
