#include "TerrainGen.h"

double TerrainGen::fractal_noise(double x, double y, int octaves, double persistence, double lacunarity) {
    double total = 0;
    double frequency = 1.0;
    double amplitude = 1.0;
    double max_value = 0;  // Used for normalizing result
    
    for(int i = 0; i < octaves; i++) {
        total += noise.eval(x * frequency, y * frequency) * amplitude;
        
        max_value += amplitude;
        amplitude *= persistence;
        frequency *= lacunarity;
    }
    
    // Normalize the result
    return total / max_value;
}

double TerrainGen::warped_noise(double x, double y, double warp_strength) {
    // Use noise to warp the input coordinates
    double warp_x = noise.eval(x * 0.5, y * 0.5) * warp_strength;
    double warp_y = noise.eval(x * 0.5 + 100, y * 0.5 + 100) * warp_strength;
    
    // Generate noise with warped coordinates
    return noise.eval(x + warp_x, y + warp_y);
}

double TerrainGen::ridged_noise(double x, double y, int octaves, double persistence, double lacunarity) {
    double total = 0;
    double frequency = 1.0;
    double amplitude = 1.0;
    double max_value = 0;
    
    for(int i = 0; i < octaves; i++) {
        // Get absolute value of noise and invert it (1 - abs)
        double value = 1.0 - std::abs(noise.eval(x * frequency, y * frequency));
        
        // Square the value to increase the ridge effect
        value *= value;
        
        total += value * amplitude;
        max_value += amplitude;
        
        amplitude *= persistence;
        frequency *= lacunarity;
    }
    
    return total / max_value;
}

double TerrainGen::terraced_noise(double x, double y, int terraces, double terrace_roughness) {
    double value = fractal_noise(x, y, 6, 0.5, 2.0);
    
    // Add terrace effect
    value *= terraces;
    double integer_part = std::floor(value);
    double fractional_part = value - integer_part;
    
    // Smoothstep the fractional part for terracing
    fractional_part = fractional_part * fractional_part * (3.0 - 2.0 * fractional_part);
    
    // Blend between pure terracing and smooth transitions
    fractional_part = terrace_roughness * fractional_part + 
                      (1.0 - terrace_roughness) * (value - integer_part);
    
    return (integer_part + fractional_part) / terraces;
}

std::vector<std::vector<TerrainCell>> TerrainGen::generate_terrain(int width, int height, const TerrainParams& params) {
    std::vector<std::vector<TerrainCell>> terrain(height, std::vector<TerrainCell>(width));
    OpenSimplexNoise::Noise noise(seed);
    
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            TerrainCell& cell = terrain[y][x];
            
            // Base height with fractals
            cell.height = fractal_noise(
                x * params.height_scale, 
                y * params.height_scale, 
                params.height_octaves, 
                params.height_persistence, 
                params.height_lacunarity
            );
            
            // Continentalness (large, smooth features)
            cell.continentalness = fractal_noise(
                x * params.continental_scale, 
                y * params.continental_scale, 
                params.continental_octaves, 
                0.3, 2.0
            );
            
            // Apply continentalness to height (creates landmasses and oceans)
            // This is a key part for controlling large landmass distribution
            double continental_influence = (cell.continentalness * 2.0) - 0.5;
            cell.height = cell.height * 0.3 + continental_influence * 0.7;
            
            // Generate erosion map
            cell.erosion = fractal_noise(
                x * params.erosion_scale + 500, 
                y * params.erosion_scale + 500, 
                params.erosion_octaves, 
                0.6, 2.0
            );
            
            // Temperature decreases with latitude and altitude
            double latitude_factor = std::abs(height/2.0 - y) / (height/2.0);
            cell.temperature = fractal_noise(
                x * params.temperature_scale + 1000, 
                y * params.temperature_scale + 1000, 
                3, 0.5, 2.0
            );
            // Mix with latitude influence
            cell.temperature = cell.temperature * (1.0 - params.temp_latitude_influence) + 
                              (1.0 - latitude_factor) * params.temp_latitude_influence;
            // Altitude reduces temperature
            cell.temperature -= cell.height * 0.2;
            
            // Humidity
            cell.humidity = fractal_noise(
                x * params.humidity_scale + 2000, 
                y * params.humidity_scale + 2000, 
                params.humidity_octaves, 
                0.4, 2.0
            );
            
            // Optional: Add mountain peaks with ridged noise
            double peaks = ridged_noise(
                x * 0.01 + 3000, 
                y * 0.01 + 3000, 
                4, 0.5, 2.0
            ) * params.peaks_influence;
            
            // Apply erosion - reduces height in erosion-heavy areas
            cell.height -= cell.erosion * 0.2;
            
            // Add mountain peaks to height
            cell.height += peaks;
        }
    }
    
    return terrain;
}