# Make sure the build directory is in your Python path
import os
import sys
import time

# Add the absolute path to the build directory
build_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../build'))
sys.path.insert(0, build_dir)

# Now import the renamed module
import cpyxdemo

# Create an instance of the renamed class
demo = cpyxdemo.CPYxDemo()

# Call the print_world method
demo.print_world()

# Show available methods
print("Available methods:", [method for method in dir(demo) if not method.startswith('__')])

# Test parameters
count_to = 1_000_000_000  # 1 billion (adjust as needed for your machine)

# Python implementation of the count function
def python_count(n):
    sum = 0
    for i in range(n):
        sum += i
    return sum

# Timing test for C++ count function
print(f"\nCounting to {count_to} using C++...")
start_time = time.time()
cpp_result = demo.count(count_to)
cpp_time = time.time() - start_time
print(f"C++ result: {cpp_result}")
print(f"C++ time: {cpp_time:.4f} seconds")

# Timing test for Python count function
print(f"\nCounting to {count_to} using Python...")
start_time = time.time()
py_result = python_count(count_to)
py_time = time.time() - start_time
print(f"Python result: {py_result}")
print(f"Python time: {py_time:.4f} seconds")

# Compare the performance
print(f"\nPython is {py_time/cpp_time:.2f}x slower than C++")
print("\nPython script completed successfully!")

demo.print_hello()