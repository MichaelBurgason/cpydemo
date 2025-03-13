#pragma once

#include <vector>
#include <string>
#include <map>
#include "openSimplexNoise/OpenSimplexNoise/OpenSimplexNoise.h"

// Enum defining the terrain parameters that can be modified by transfer functions
enum class TerrainParameter {
    BaseHeight,
    Continentalness,
    RidgeFormation,
    Erosion,
    Temperature,
    Humidity
};

// Structure for a transfer function control point (input→output mapping)
struct TransferControlPoint {
    double inputValue;    // Input value in -1.0 to 1.0 range
    double outputValue;   // Output value in -1.0 to 1.0 range
    double weight;        // Weight of this control point (default = 1.0)
};

// Structure for a complete transfer function spline
struct TransferSpline {
    std::vector<TransferControlPoint> controlPoints;  // Control points defining the transfer function
    double tension = 0.5;        // Controls how "tight" the curve is (0-1)
    double weight = 1.0;         // Global weight of this spline when multiple are applied
    std::string name = "Transfer"; // Name for identification
    TerrainParameter parameter;  // The parameter this transfer function affects
};

// Legacy structure for backward compatibility (spatial splines)
struct SplineControlPoint {
    int x;              // X position in the noise map
    int y;              // Y position in the noise map
    double value;       // Target parameter value at this point (-1.0 to 1.0 range)
    double weight;      // Influence weight (0-1)
    double radius;      // Radius of influence
};

// Legacy structure for backward compatibility
struct TerrainSpline {
    std::vector<SplineControlPoint> controlPoints;
    double tension = 0.5;         // Controls how "tight" the curve is (0-1)
    double bias = 0.0;            // Controls the direction of the curve (-1 to 1)
    double continuity = 0.0;      // Controls the smoothness of the curve (-1 to 1)
    double influence = 1.0;       // Global influence of this spline (0-1)
    std::string name = "Spline";  // Name for identification
    TerrainParameter parameter;   // The parameter this spline affects
    bool intensify = true;        // If true, strengthens the parameter values; if false, weakens them
};

// Structure to hold parameters for terrain generation
struct TerrainParams {
    // Base terrain noise parameters
    int height_octaves = 6;
    double height_persistence = 0.5;
    double height_lacunarity = 2.0;
    double height_scale = 0.005;
    
    // Continentalness parameters
    double continental_scale = 0.001;
    int continental_octaves = 2;
    
    // Ridge formation parameters
    double ridge_scale = 0.003;
    int ridge_octaves = 4;
    double ridge_persistence = 0.5;
    double ridge_influence = 0.6;
    
    // Erosion parameters
    double erosion_scale = 0.01;
    int erosion_octaves = 4;
    
    // Temperature parameters
    double temperature_scale = 0.003;
    double temp_latitude_influence = 0.7;
    
    // Humidity parameters
    double humidity_scale = 0.004;
    int humidity_octaves = 3;
    
    // Peaks intensity
    double peaks_influence = 0.3;
    
    // Transfer function parameters
    std::map<TerrainParameter, std::vector<TransferSpline>> transferFunctions;
    
    // Legacy spline parameters (kept for backward compatibility)
    std::map<TerrainParameter, std::vector<TerrainSpline>> parameterSplines;
    std::map<TerrainParameter, double> splineInfluence = {
        {TerrainParameter::BaseHeight, 0.5},
        {TerrainParameter::Continentalness, 0.7}, 
        {TerrainParameter::RidgeFormation, 0.8},
        {TerrainParameter::Erosion, 0.4},
        {TerrainParameter::Temperature, 0.3},
        {TerrainParameter::Humidity, 0.5}
    };
    
    // Auto-generation options
    bool auto_generate_splines = false; // Disable legacy spline auto-generation
    bool auto_generate_transfer_functions = true; // Enable transfer function auto-generation
    std::map<TerrainParameter, int> auto_spline_count = {
        {TerrainParameter::BaseHeight, 3},
        {TerrainParameter::Continentalness, 2}, 
        {TerrainParameter::RidgeFormation, 4},
        {TerrainParameter::Erosion, 3},
        {TerrainParameter::Temperature, 2},
        {TerrainParameter::Humidity, 2}
    };
    
    // Animation steps
    int animation_frames = 8;
};

// Enumeration for biome types
enum class BiomeType {
    DeepOcean, Ocean, WarmOcean, FrozenOcean, 
    Beach, RockyShore, Marsh,
    Desert, Savanna, Shrubland, TropicalShrubland, 
    Plains, TemperateGrassland, BorealPlains, Tundra,
    TropicalRainforest, TropicalSeasonalForest, TemperateRainforest, 
    DeciduousForest, Forest, MixedForest, Taiga, SnowyTaiga,
    Hills, Badlands, Mountain, SnowyMountain, SnowyCap, IceSheet
};

// Structure to hold terrain data for a single cell
struct TerrainCell {
    double baseNoise;       // Raw noise value (-1 to 1) before parameter modifications
    double height;          // Final height value after all modifications
    
    // Raw parameter values from original noise
    double rawTemperature;     // Raw temperature value (-1 to 1)
    double rawHumidity;        // Raw humidity value (-1 to 1)
    double rawErosion;         // Raw erosion value (-1 to 1)
    double rawContinentalness; // Raw continentalness value (-1 to 1)
    double rawRidgeValue;      // Raw ridge formation value (-1 to 1)
    
    // Transformed parameter values after applying transfer functions
    double temperature;     // Temperature value (-1 to 1, modified by transfer function)
    double humidity;        // Humidity value (-1 to 1, modified by transfer function)
    double erosion;         // Erosion value (-1 to 1, modified by transfer function)
    double continentalness; // Continentalness value (-1 to 1, modified by transfer function)
    double ridgeValue;      // Ridge formation value (-1 to 1, modified by transfer function)
    
    BiomeType biome;        // Final biome classification
};

class TerrainGen2DAnimation {
public:
    TerrainGen2DAnimation(int64_t seed = 16546746132);
    ~TerrainGen2DAnimation() = default;
    
    // Generate a 3D array where:
    // - First dimension: height (y)
    // - Second dimension: width (x)
    // - Third dimension: animation frames
    // Each frame represents a stage in the terrain generation process
    std::vector<std::vector<std::vector<TerrainCell>>> generateAnimatedTerrain(
        int width, 
        int height, 
        const TerrainParams& params = TerrainParams()
    );
    
    // Get the name of each animation frame for display
    std::vector<std::string> getAnimationFrameNames(const TerrainParams& params = TerrainParams());
    
    // Utility functions for terrain generation
    static BiomeType classifyBiome(const TerrainCell& cell);
    
    // Generate default transfer functions for various terrain parameters
    static std::vector<TransferSpline> generateDefaultTransferFunction(TerrainParameter parameter);
    
private:
    OpenSimplexNoise::Noise noise;
    
    // Noise generation functions
    double fractalNoise(double x, double y, int octaves, double persistence, double lacunarity);
    double ridgedNoise(double x, double y, int octaves, double persistence, double lacunarity);
    
    // Terrain generation stages
    void addBaseNoiseLayer(
        std::vector<std::vector<std::vector<TerrainCell>>>& terrain,
        int width, int height, const TerrainParams& params, int frameIndex
    );
    
    void applyTransferFunction(
        std::vector<std::vector<std::vector<TerrainCell>>>& terrain,
        int width, int height, const TerrainParams& params, 
        TerrainParameter parameter, int frameIndex
    );
    
    void calculateFinalHeight(
        std::vector<std::vector<std::vector<TerrainCell>>>& terrain,
        int width, int height, const TerrainParams& params, int frameIndex
    );
    
    void classifyBiomes(
        std::vector<std::vector<std::vector<TerrainCell>>>& terrain,
        int width, int height, int frameIndex
    );
    
    // Transfer function utilities
    double evaluateTransferSpline(
        const TransferSpline& spline, double inputValue
    );
    
    double cardinalSpline(
        double p0, double p1, double p2, double p3, double t, double tension
    );
    
    void updateHeightFromParameter(
        TerrainCell& cell, TerrainParameter parameter
    );
};