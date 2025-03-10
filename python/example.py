import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Make sure the build directory is in your Python path
# Add the absolute path to the build directory
build_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../build'))
sys.path.insert(0, build_dir)

# Now import the renamed module
import noise_demo

# Create an instance of the renamed class
demo = noise_demo.CPYDemo()

# Call the print_world method
demo.print_world()

# Show available methods
print("Available methods:", [method for method in dir(demo) if not method.startswith('__')])

# Test parameters
count_to = 1_000 # 1 billion (adjust as needed for your machine)

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
print(demo.test_function(3.14))

# --------------------------------------------------------
# Noise Map Generation and Visualization
# --------------------------------------------------------

def generate_noise_image(width=512, height=512, scale=0.1, offsetX=0.0, offsetY=0.0, output_filename="noise_map.png"):
    """
    Generate a noise map and save it as an image with a black to white gradient.
    
    Parameters:
    -----------
    width : int
        Width of the noise map
    height : int
        Height of the noise map
    scale : float
        Scale of the noise (higher values = more zoomed out noise)
    offsetX : float
        X offset for the noise sampling
    offsetY : float
        Y offset for the noise sampling
    output_filename : str
        Filename for the output image
    """
    # Generate the noise map using the provided function
    print(f"Generating noise map ({width}x{height}) with scale={scale}...")
    start_time = time.time()
    noise_map = demo.generate_noise_map(width, height, scale, offsetX, offsetY)
    generation_time = time.time() - start_time
    print(f"Noise map generation time: {generation_time:.4f} seconds")
    
    # Convert the C++ vector<vector<double>> to numpy array
    noise_map = np.array(noise_map, dtype=np.float64)
    
    # Create a simple black to white colormap
    colors = [(0, 0, 0), (1, 1, 1)]  # Black to white
    cmap = LinearSegmentedColormap.from_list("BlackToWhite", colors)
    
    # Create the plot
    plt.figure(figsize=(10, 10))
    plt.imshow(noise_map, cmap=cmap)
    plt.axis('off')  # Hide the axes
    plt.title(f"Noise Map (scale={scale}, offset=({offsetX}, {offsetY}))")
    
    # Save the image
    plt.savefig(output_filename, bbox_inches='tight', pad_inches=0, dpi=300)
    print(f"Noise map saved to {output_filename}")
    
    # Display the image
    plt.show()

print("\n" + "="*50)
print("NOISE MAP GENERATION")
print("="*50)

# Generate a single noise map
generate_noise_image(
    width=512,
    height=512,
    scale=0.05,
    offsetX=42.0,
    offsetY=17.0,
    output_filename="perlin_noise_map.png"
)

# Generate multiple noise maps with different scales
scales = [0.01, 0.05, 0.1, 0.2]
for i, scale in enumerate(scales):
    generate_noise_image(
        width=256,  # Smaller size for faster generation of multiple maps
        height=256,
        scale=scale,
        offsetX=42.0,
        offsetY=17.0,
        output_filename=f"noise_map_scale_{scale:.2f}.png"
    )

print("\nNoise map generation completed successfully!")