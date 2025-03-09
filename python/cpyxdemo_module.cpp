#include <pybind11/pybind11.h>
#include "CPYDemo.h"  // This should match your actual header filename

namespace py = pybind11;

PYBIND11_MODULE(cpyxdemo, m) {  // Changed module name to cpyxdemo
    m.doc() = "CPYxDemo Python binding";
    
    py::class_<CPYDemo>(m, "CPYxDemo")  // Changed class name to CPYxDemo
        .def(py::init<>())
        .def("print_world", &CPYDemo::printWorld, "Prints 'Hello World' to console")
        .def("count", &CPYDemo::count, "Counts up to n and returns the sum",
             py::arg("n"))
        .def("print_hello", &CPYDemo::printHello, "Prints 'hello'");
}