#ifndef TERRAIN_GEN_H
#define TERRAIN_GEN_H

#include "openSimplexNoise/OpenSimplexNoise/OpenSimplexNoise.h"
#include <vector>

struct TerrainParams {
    // Base terrain
    int height_octaves = 6;
    double height_persistence = 0.5;
    double height_lacunarity = 2.0;
    double height_scale = 0.005;
    
    // Continentalness (large-scale features)
    double continental_scale = 0.001;
    int continental_octaves = 2;
    
    // Erosion
    double erosion_scale = 0.01;
    int erosion_octaves = 4;
    
    // Temperature
    double temperature_scale = 0.003;
    double temp_latitude_influence = 0.7;
    
    // Humidity
    double humidity_scale = 0.004;
    int humidity_octaves = 3;
    
    // Peaks intensity
    double peaks_influence = 0.3;
};

struct TerrainCell {
    double height;
    double temperature;
    double humidity;
    double erosion;
    double continentalness;
};



class TerrainGen {
public:
    /**
     * @brief Generates fractal noise using multiple octaves
     * @param x X-coordinate in noise space
     * @param y Y-coordinate in noise space
     * @param octaves Number of noise layers to combine
     * @param persistence How much each octave's amplitude decreases
     * @param lacunarity How much each octave's frequency increases
     * @return Normalized fractal noise value between -1 and 1
     */
    double fractal_noise(double x, double y, int octaves, double persistence, double lacunarity);

    /**
     * @brief Generates warped noise by using noise to distort input coordinates
     * @param x X-coordinate in noise space
     * @param y Y-coordinate in noise space
     * @param warp_strength How much the coordinates should be distorted
     * @return Warped noise value between -1 and 1
     */
    double warped_noise(double x, double y, double warp_strength);

    /**
     * @brief Generates ridged noise for mountain-like features
     * @param x X-coordinate in noise space
     * @param y Y-coordinate in noise space
     * @param octaves Number of noise layers to combine
     * @param persistence How much each octave's amplitude decreases
     * @param lacunarity How much each octave's frequency increases
     * @return Normalized ridged noise value between 0 and 1
     */
    double ridged_noise(double x, double y, int octaves, double persistence, double lacunarity);

    /**
     * @brief Generates terraced noise for plateau-like features
     * @param x X-coordinate in noise space
     * @param y Y-coordinate in noise space
     * @param terraces Number of distinct height levels
     * @param terrace_roughness How smooth the transitions between terraces are (0-1)
     * @return Normalized terraced noise value between 0 and 1
     */
    double terraced_noise(double x, double y, int terraces, double terrace_roughness);

    /**
     * @brief Generates a complete terrain with multiple layers of features
     * @param width Width of the terrain grid
     * @param height Height of the terrain grid
     * @param params Terrain generation parameters
     * @return 2D vector of TerrainCell objects containing height and biome data
     */
    std::vector<std::vector<TerrainCell>> generate_terrain(int width, int height, const TerrainParams& params);

private:
    OpenSimplexNoise::Noise noise;
};

#endif // TERRAIN_GEN_H