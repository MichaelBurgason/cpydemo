#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>  // Include numpy header
#include <pybind11/stl.h>  // Include stl header
#include "CPYDemo.h"  // Include the header file

namespace py = pybind11;

PYBIND11_MODULE(noise_demo, m) {  // Changed module name to noisedemo
    m.doc() = "Noisedemo Python binding";
    
    py::class_<CPYDemo>(m, "CPYDemo")  // Changed class name to CPYDemo
        .def(py::init<>())
        .def("print_world", &CPYDemo::printWorld, "Prints 'Hello World' to console")
        .def("count", &CPYDemo::count, "Counts up to n and returns the sum",
             py::arg("n"))
        .def("print_hello", &CPYDemo::printHello, "Prints 'hello'")
        .def("test_function", &CPYDemo::testFunction, "returns the passed double",
             py::arg("n"))
        .def("generate_noise_map", &CPYDemo::generate_noise_map, 
             "Generate a 2D noise map",
             py::arg("width"), py::arg("height"), py::arg("scale"),
             py::arg("offsetX"), py::arg("offsetY"));
}