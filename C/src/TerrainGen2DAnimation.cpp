#include "TerrainGen2DAnimation.h"
#include <cmath>
#include <algorithm>
#include <random>

TerrainGen2DAnimation::TerrainGen2DAnimation(int64_t seed) : noise(seed) {
}

std::vector<std::vector<std::vector<TerrainCell>>> TerrainGen2DAnimation::generateAnimatedTerrain(
    int width, 
    int height, 
    const TerrainParams& params) {
    
    // Create 3D vector to store animation frames
    int frames = params.animation_frames;
    std::vector<std::vector<std::vector<TerrainCell>>> animatedTerrain(
        height, 
        std::vector<std::vector<TerrainCell>>(
            width, 
            std::vector<TerrainCell>(frames)
        )
    );
    
    // Generate each layer for the animation - each step builds upon the previous
    addBaseNoiseLayer(animatedTerrain, width, height, params, 0);  // Frame 0: Base noise values
    
    if (frames > 1)
        applyTransferFunction(animatedTerrain, width, height, params, TerrainParameter::Continentalness, 1);
        
    if (frames > 2)
        applyTransferFunction(animatedTerrain, width, height, params, TerrainParameter::RidgeFormation, 2);
        
    if (frames > 3)
        applyTransferFunction(animatedTerrain, width, height, params, TerrainParameter::Erosion, 3);
        
    if (frames > 4)
        applyTransferFunction(animatedTerrain, width, height, params, TerrainParameter::Temperature, 4);
        
    if (frames > 5)
        applyTransferFunction(animatedTerrain, width, height, params, TerrainParameter::Humidity, 5);
        
    if (frames > 6)
        calculateFinalHeight(animatedTerrain, width, height, params, 6);
    
    // Final frame with all layers combined and biomes assigned
    if (frames > 7)
        classifyBiomes(animatedTerrain, width, height, frames - 1);
    
    return animatedTerrain;
}

std::vector<std::string> TerrainGen2DAnimation::getAnimationFrameNames(const TerrainParams& params) {
    std::vector<std::string> frameNames;
    
    frameNames.push_back("Base Noise Values");
    frameNames.push_back("Continentalness Transfer Applied");
    frameNames.push_back("Ridge Formation Transfer Applied");
    frameNames.push_back("Erosion Transfer Applied");
    frameNames.push_back("Temperature Distribution");
    frameNames.push_back("Humidity Distribution");
    frameNames.push_back("Final Height Calculation");
    frameNames.push_back("Biome Classification");
    
    // Ensure we have the right number of names
    while (frameNames.size() > params.animation_frames) {
        frameNames.pop_back();
    }
    
    return frameNames;
}

// Fractal Brownian Motion noise implementation (output range: -1 to 1)
double TerrainGen2DAnimation::fractalNoise(double x, double y, int octaves, double persistence, double lacunarity) {
    double total = 0;
    double frequency = 1.0;
    double amplitude = 1.0;
    double maxValue = 0;  // Used for normalizing result
    
    for(int i = 0; i < octaves; i++) {
        total += noise.eval(x * frequency, y * frequency) * amplitude;
        
        maxValue += amplitude;
        amplitude *= persistence;
        frequency *= lacunarity;
    }
    
    // Normalize the result to -1 to 1 range
    return total / maxValue;
}

// Ridged noise for mountain ranges (output range: 0 to 1)
double TerrainGen2DAnimation::ridgedNoise(double x, double y, int octaves, double persistence, double lacunarity) {
    double total = 0;
    double frequency = 1.0;
    double amplitude = 1.0;
    double maxValue = 0;
    
    for(int i = 0; i < octaves; i++) {
        // Get absolute value of noise and invert it (1 - abs)
        double value = 1.0 - std::abs(noise.eval(x * frequency, y * frequency));
        
        // Square the value to increase the ridge effect
        value *= value;
        
        total += value * amplitude;
        maxValue += amplitude;
        
        amplitude *= persistence;
        frequency *= lacunarity;
    }
    
    // Scale to -1 to 1 range (by scaling 0-1 to -1 to 1)
    return (total / maxValue) * 2.0 - 1.0;
}

// Generate the base noise layers for all parameters
void TerrainGen2DAnimation::addBaseNoiseLayer(
    std::vector<std::vector<std::vector<TerrainCell>>>& terrain,
    int width, int height, const TerrainParams& params, int frameIndex) {
    
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            TerrainCell& cell = terrain[y][x][frameIndex];
            
            // Generate raw noise values for each parameter in the -1 to 1 range
            // Base height noise
            cell.baseNoise = fractalNoise(
                x * params.height_scale, 
                y * params.height_scale, 
                params.height_octaves, 
                params.height_persistence, 
                params.height_lacunarity
            );
            
            // Continentalness noise
            cell.rawContinentalness = fractalNoise(
                x * params.continental_scale, 
                y * params.continental_scale, 
                params.continental_octaves, 
                0.3, 2.0
            );
            // Initialize transformed value with raw
            cell.continentalness = cell.rawContinentalness;
            
            // Ridge formation noise
            cell.rawRidgeValue = ridgedNoise(
                x * params.ridge_scale + 500, 
                y * params.ridge_scale + 500, 
                params.ridge_octaves, 
                params.ridge_persistence, 
                2.0
            );
            // Initialize transformed value with raw
            cell.ridgeValue = cell.rawRidgeValue;
            
            // Erosion noise
            cell.rawErosion = fractalNoise(
                x * params.erosion_scale + 1000, 
                y * params.erosion_scale + 1000, 
                params.erosion_octaves, 
                0.6, 2.0
            );
            // Initialize transformed value with raw
            cell.erosion = cell.rawErosion;
            
            // Temperature - base noise plus latitude influence
            double latitudeFactor = (static_cast<double>(y) / height) * 2.0 - 1.0; // -1 at top, +1 at bottom
            cell.rawTemperature = fractalNoise(
                x * params.temperature_scale + 1500, 
                y * params.temperature_scale + 1500, 
                3, 0.5, 2.0
            );
            
            // Mix with latitude influence
            cell.rawTemperature = cell.rawTemperature * (1.0 - params.temp_latitude_influence) + 
                             latitudeFactor * params.temp_latitude_influence;
            // Initialize transformed value with raw
            cell.temperature = cell.rawTemperature;
            
            // Humidity
            cell.rawHumidity = fractalNoise(
                x * params.humidity_scale + 2000, 
                y * params.humidity_scale + 2000, 
                params.humidity_octaves, 
                0.4, 2.0
            );
            // Initialize transformed value with raw
            cell.humidity = cell.rawHumidity;
            
            // Set initial height based only on base noise
            cell.height = cell.baseNoise;
            
            // Initialize biome as default
            cell.biome = BiomeType::Plains;
            
            // Copy to next frame as starting point if not the last frame
            if (frameIndex < params.animation_frames - 1) {
                terrain[y][x][frameIndex + 1] = cell;
            }
        }
    }
}

// Apply transfer function to modify a specific terrain parameter
void TerrainGen2DAnimation::applyTransferFunction(
    std::vector<std::vector<std::vector<TerrainCell>>>& terrain,
    int width, int height, const TerrainParams& params, 
    TerrainParameter parameter, int frameIndex) {
    
    // Get transfer splines for this parameter
    std::vector<TransferSpline> transferSplines;
    if (params.transferFunctions.find(parameter) != params.transferFunctions.end()) {
        transferSplines = params.transferFunctions.at(parameter);
    }
    
    // Auto-generate transfer splines if needed
    if (params.auto_generate_transfer_functions && transferSplines.empty()) {
        transferSplines = generateDefaultTransferFunction(parameter);
    }
    
    // Apply the transfer functions to transform the raw parameter values
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            TerrainCell& prevCell = terrain[y][x][frameIndex - 1];
            TerrainCell& cell = terrain[y][x][frameIndex];
            
            // Copy all data from previous frame
            cell = prevCell;
            
            // Get the raw value for the specified parameter
            double rawValue = 0.0;
            switch (parameter) {
                case TerrainParameter::BaseHeight:
                    rawValue = cell.baseNoise;
                    break;
                case TerrainParameter::Continentalness:
                    rawValue = cell.rawContinentalness;
                    break;
                case TerrainParameter::RidgeFormation:
                    rawValue = cell.rawRidgeValue;
                    break;
                case TerrainParameter::Erosion:
                    rawValue = cell.rawErosion;
                    break;
                case TerrainParameter::Temperature:
                    rawValue = cell.rawTemperature;
                    break;
                case TerrainParameter::Humidity:
                    rawValue = cell.rawHumidity;
                    break;
            }
            
            // Apply transfer functions to transform the raw value
            double transformedValue = rawValue; // Default to raw value
            
            if (!transferSplines.empty()) {
                // Apply each transfer function and blend the results
                double totalWeight = 0.0;
                double weightedSum = 0.0;
                
                for (const auto& spline : transferSplines) {
                    double splineOutput = evaluateTransferSpline(spline, rawValue);
                    
                    // Apply the weight of this transfer function
                    double weight = spline.weight;
                    weightedSum += splineOutput * weight;
                    totalWeight += weight;
                }
                
                // Calculate the final transformed value as the weighted average
                if (totalWeight > 0.0) {
                    transformedValue = weightedSum / totalWeight;
                }
            }
            
            // Clamp the transformed value to valid range
            transformedValue = std::max(-1.0, std::min(1.0, transformedValue));
            
            // Update the specific parameter with the transformed value
            switch (parameter) {
                case TerrainParameter::BaseHeight:
                    cell.baseNoise = transformedValue;
                    break;
                case TerrainParameter::Continentalness:
                    cell.continentalness = transformedValue;
                    break;
                case TerrainParameter::RidgeFormation:
                    cell.ridgeValue = transformedValue;
                    break;
                case TerrainParameter::Erosion:
                    cell.erosion = transformedValue;
                    break;
                case TerrainParameter::Temperature:
                    cell.temperature = transformedValue;
                    break;
                case TerrainParameter::Humidity:
                    cell.humidity = transformedValue;
                    break;
            }
            
            // Update height visualization based on the current parameter
            updateHeightFromParameter(cell, parameter);
        }
    }
    
    // Copy to next frame if not the last frame
    if (frameIndex < params.animation_frames - 1) {
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                terrain[y][x][frameIndex + 1] = terrain[y][x][frameIndex];
            }
        }
    }
}

// Evaluate a transfer spline for a given input value
double TerrainGen2DAnimation::evaluateTransferSpline(const TransferSpline& spline, double inputValue) {
    const auto& points = spline.controlPoints;
    
    // Need a minimum number of points for interpolation
    if (points.size() < 2) {
        return inputValue; // Just return the input value if insufficient control points
    }
    
    // Find the control points to interpolate between
    size_t lowerIndex = 0;
    size_t upperIndex = 1;
    
    // Find the segment that contains the input value
    for (size_t i = 0; i < points.size() - 1; i++) {
        if (inputValue >= points[i].inputValue && inputValue <= points[i + 1].inputValue) {
            lowerIndex = i;
            upperIndex = i + 1;
            break;
        }
    }
    
    // Special case: input value is outside the range of control points
    if (inputValue < points.front().inputValue) {
        return points.front().outputValue;
    }
    if (inputValue > points.back().inputValue) {
        return points.back().outputValue;
    }
    
    // Calculate interpolation parameter t (0 to 1)
    double t = (inputValue - points[lowerIndex].inputValue) / 
               (points[upperIndex].inputValue - points[lowerIndex].inputValue);
    
    // Find the indices for cardinal spline interpolation
    // (use boundary points if needed)
    size_t p0Index = (lowerIndex > 0) ? lowerIndex - 1 : lowerIndex;
    size_t p1Index = lowerIndex;
    size_t p2Index = upperIndex;
    size_t p3Index = (upperIndex < points.size() - 1) ? upperIndex + 1 : upperIndex;
    
    // Get the output values for interpolation
    double p0 = points[p0Index].outputValue;
    double p1 = points[p1Index].outputValue;
    double p2 = points[p2Index].outputValue;
    double p3 = points[p3Index].outputValue;
    
    // Perform cardinal spline interpolation with tension
    return cardinalSpline(p0, p1, p2, p3, t, spline.tension);
}

// Helper function to update the height value based on a parameter change
// This provides a visual representation of how each parameter affects the terrain
void TerrainGen2DAnimation::updateHeightFromParameter(TerrainCell& cell, TerrainParameter parameter) {
    switch (parameter) {
        case TerrainParameter::BaseHeight:
            cell.height = cell.baseNoise;
            break;
        case TerrainParameter::Continentalness:
            // Higher continentalness raises land, lower lowers it
            cell.height = cell.baseNoise * 0.3 + cell.continentalness * 0.7;
            break;
        case TerrainParameter::RidgeFormation:
            // Ridges add height where ridgeValue is positive
            cell.height = cell.height + cell.ridgeValue * 0.3;
            break;
        case TerrainParameter::Erosion:
            // Erosion lowers height in high-erosion areas
            cell.height = cell.height - cell.erosion * 0.15;
            break;
        case TerrainParameter::Temperature:
            // Temperature visualization doesn't change height
            break;
        case TerrainParameter::Humidity:
            // Humidity visualization doesn't change height
            break;
    }
    
    // Clamp height to valid range
    cell.height = std::max(-1.0, std::min(1.0, cell.height));
}

// Calculate final height using all parameters
void TerrainGen2DAnimation::calculateFinalHeight(
    std::vector<std::vector<std::vector<TerrainCell>>>& terrain,
    int width, int height, const TerrainParams& params, int frameIndex) {
    
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            TerrainCell& prevCell = terrain[y][x][frameIndex - 1];
            TerrainCell& cell = terrain[y][x][frameIndex];
            
            // Copy all data from previous frame
            cell = prevCell;
            
            // Calculate final height using all transformed parameters
            // Start with base height from continentalness
            double finalHeight = cell.baseNoise //* 0.4 + cell.continentalness * 0.6;
            
            // Apply ridge formations to add mountains
            finalHeight += cell.ridgeValue * params.ridge_influence;
            
            // Apply erosion to create valleys
            finalHeight -= cell.erosion * 0.2;
            
            // Small humidity influence - higher humidity slightly lowers peaks
            finalHeight -= std::max(0.0, cell.humidity) * 0.05 * finalHeight;
            
            // Clamp final height
            cell.height = std::max(-1.0, std::min(1.0, finalHeight));
        }
    }
    
    // Copy to next frame if not the last frame
    if (frameIndex < params.animation_frames - 1) {
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                terrain[y][x][frameIndex + 1] = terrain[y][x][frameIndex];
            }
        }
    }
}

// Generate default transfer functions for various terrain parameters
std::vector<TransferSpline> TerrainGen2DAnimation::generateDefaultTransferFunction(TerrainParameter parameter) {
    std::vector<TransferSpline> transferSplines;
    TransferSpline spline;
    
    // Set name and parameter
    spline.parameter = parameter;
    spline.weight = 1.0;
    spline.tension = 0.5;
    
    // Configure control points based on the parameter type
    switch (parameter) {
        case TerrainParameter::BaseHeight:
            spline.name = "Default Base Height Transfer";
            // Simple linear mapping
            spline.controlPoints.push_back({-1.0, -1.0, 1.0});
            spline.controlPoints.push_back({-0.5, -0.5, 1.0});
            spline.controlPoints.push_back({0.0, 0.0, 1.0});
            spline.controlPoints.push_back({0.5, 0.5, 1.0});
            spline.controlPoints.push_back({1.0, 1.0, 1.0});
            break;
            
        case TerrainParameter::Continentalness:
            spline.name = "Deep Ocean Trench Transfer";
            
            // Deep ocean trench
            spline.controlPoints.push_back({-1.0, -0.8, 1.0});    // Absolute deepest point
            spline.controlPoints.push_back({-0.9, -1, 2.0});      // Very deep trench
            spline.controlPoints.push_back({-0.8, -0.8, 1.0});    // Start of trench
            
            // Underwater ridge or shelf feature
            spline.controlPoints.push_back({-0.6, -0.6, 1.0});    // Underwater shelf
            spline.controlPoints.push_back({-0.55, -0.3, 1.0});   // Underwater ridge (rises a bit)
            
            // Continental shelf and deeper water
            spline.controlPoints.push_back({-0.5, -0.8, 1.0});    // Deep water
            
            // Quick drop off near shore
            spline.controlPoints.push_back({-0.3, -0.5, 1.0});    // Getting deep quickly  
            spline.controlPoints.push_back({-0.2, -0.2, 1.0});    // Shallow water
            
            // Land and shoreline
            spline.controlPoints.push_back({-0.1, 0.0, 1.0});     // Shoreline
            spline.controlPoints.push_back({0.0, 0.1, 1.0});      // Lowland
            spline.controlPoints.push_back({0.4, 0.2, 1.0});      // Lowland
            spline.controlPoints.push_back({0.55, 0.25, 1.0});      // Highland
            spline.controlPoints.push_back({0.6, 0.4, 1.0});      // Highland
            spline.controlPoints.push_back({0.65, 0.6, 1.0});      // Highland
            spline.controlPoints.push_back({0.7, 0.5, 1.0});      // Highland
            spline.controlPoints.push_back({0.75, 0.8, 1.0});      // Highland
            spline.controlPoints.push_back({1.0, 1.0, 3.0});      // Max height (reduced from 1.0)
            break;
            
        case TerrainParameter::RidgeFormation:
            spline.name = "Ravine and Ridge Formation";
            // Negative ridge values become ravines and river valleys
            spline.controlPoints.push_back({-1.0, 1, 1.0}); // Deep ravines now become elevated
            spline.controlPoints.push_back({-0.25, 0.1, 2.0}); // 
            
            // Near-zero values remain relatively flat (keeping your pattern around 0)
            spline.controlPoints.push_back({-0.2, 0.1, 1.0}); // Minor depressions
            spline.controlPoints.push_back({0.0, 0.0, 1.0}); // Neutral areas
            spline.controlPoints.push_back({0.2, -0.1, 1.0}); // Slight elevations become slight depressions
            
            // Positive ridge values become mountains and hills, but inverted to create ridges
            spline.controlPoints.push_back({0.4, -0.2, 1.0}); // Hills become minor valleys
            spline.controlPoints.push_back({0.5, -0.4, 1.0}); // Hills become moderate valleys
            spline.controlPoints.push_back({0.7, -0.7, 1.0}); // Mountains become deep valleys
            spline.controlPoints.push_back({0.9, -0.8, 1.0}); // High mountains become deep ravines
            spline.controlPoints.push_back({1.0, -1.0, 5.0}); // Peaks become deepest ravines
            break;
            
        case TerrainParameter::Erosion:
            spline.name = "Erosion Pattern";
            // S-curve for erosion effect
            spline.controlPoints.push_back({-1.0, -0.8, 1.0});
            spline.controlPoints.push_back({-0.7, -0.6, 1.0});
            spline.controlPoints.push_back({-0.3, -0.2, 1.0});
            spline.controlPoints.push_back({0.0, 0.0, 1.0});
            spline.controlPoints.push_back({0.3, 0.2, 1.0});
            spline.controlPoints.push_back({0.7, 0.6, 1.0});
            spline.controlPoints.push_back({1.0, 0.8, 1.0});
            break;
            
        case TerrainParameter::Temperature:
            spline.name = "Temperature Bands";
            // Create distinct temperature bands
            spline.controlPoints.push_back({-1.0, -1.0, 1.0});     // Coldest
            spline.controlPoints.push_back({-0.7, -0.8, 1.0});     // Very cold
            spline.controlPoints.push_back({-0.4, -0.4, 1.0});     // Cold
            spline.controlPoints.push_back({-0.1, 0.0, 1.0});      // Cool
            spline.controlPoints.push_back({0.1, 0.2, 1.0});       // Mild
            spline.controlPoints.push_back({0.4, 0.5, 1.0});       // Warm
            spline.controlPoints.push_back({0.7, 0.8, 1.0});       // Hot
            spline.controlPoints.push_back({1.0, 1.0, 1.0});       // Very hot
            break;
            
        case TerrainParameter::Humidity:
            spline.name = "Humidity Distribution";
            // Create more extreme humidity distribution
            spline.controlPoints.push_back({-1.0, -1.0, 1.0});     // Very dry
            spline.controlPoints.push_back({-0.6, -0.8, 1.0});     // Dry 
            spline.controlPoints.push_back({-0.2, -0.3, 1.0});     // Slightly dry
            spline.controlPoints.push_back({0.0, 0.0, 1.0});       // Neutral
            spline.controlPoints.push_back({0.2, 0.3, 1.0});       // Slightly humid
            spline.controlPoints.push_back({0.6, 0.7, 1.0});       // Humid
            spline.controlPoints.push_back({1.0, 1.0, 1.0});       // Very humid
            break;
    }
    
    transferSplines.push_back(spline);
    return transferSplines;
}

// Cardinal spline interpolation between points
double TerrainGen2DAnimation::cardinalSpline(double p0, double p1, double p2, double p3, double t, double tension) {
    // Cardinal spline formula
    double t2 = t * t;
    double t3 = t2 * t;
    
    double s = (1.0 - tension) / 2.0;
    
    double h1 = 2*t3 - 3*t2 + 1;
    double h2 = -2*t3 + 3*t2;
    double h3 = t3 - 2*t2 + t;
    double h4 = t3 - t2;
    
    return h1*p1 + h2*p2 + s * (h3*(p2-p0) + h4*(p3-p1));
}

// Classify biomes based on terrain properties
void TerrainGen2DAnimation::classifyBiomes(
    std::vector<std::vector<std::vector<TerrainCell>>>& terrain,
    int width, int height, int frameIndex) {
    
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            TerrainCell& prevCell = terrain[y][x][frameIndex - 1];
            TerrainCell& cell = terrain[y][x][frameIndex];
            
            // Copy all data from previous frame
            cell = prevCell;
            
            // Assign biome based on terrain properties
            cell.biome = classifyBiome(cell);
        }
    }
}

// Biome classification based on terrain parameters
BiomeType TerrainGen2DAnimation::classifyBiome(const TerrainCell& cell) {
    // Extract all relevant parameters
    double height = cell.height;
    double temp = cell.temperature;
    double humidity = cell.humidity;
    double erosion = cell.erosion;
    double ridge = cell.ridgeValue;
    
    // Water biomes based on depth
    if (height < -0.6) {
        return BiomeType::DeepOcean;  // Add this to BiomeType enum
    }
    
    if (height < -0.1) {
        // Different ocean variants based on temperature
        if (temp > 0.3) {
            return BiomeType::WarmOcean;  // Add this to BiomeType enum
        } 
        if (temp < -0.3) {
            return BiomeType::FrozenOcean;  // Add this to BiomeType enum
        }
        return BiomeType::Ocean;
    }
    
    // Coastal zones (beaches, marshes, etc.)
    if (height < -0.05) {
        // Marsh/swamp in warm humid areas
        if (temp > 0.2 && humidity > 0.3) {
            return BiomeType::Marsh;  // Add this to BiomeType enum
        }
        // Cold rocky shores
        if (temp < -0.2) {
            return BiomeType::RockyShore;  // Add this to BiomeType enum
        }
        // Regular beach for other cases
        return BiomeType::Beach;
    }
    
    // High mountains and peaks
    if (height > 0.7) {
        if (height > 0.85) {
            return BiomeType::SnowyCap;
        }
        if (temp < -0.3) {
            // Cold mountains with snow
            return BiomeType::SnowyMountain;  // Add this to BiomeType enum
        }
        // Regular mountains
        return BiomeType::Mountain;
    }
    
    // Hills and elevated terrain (but not mountains)
    if (height > 0.4 || ridge > 0.5) {
        // Highland with high erosion becomes badlands
        if (erosion > 0.4 && humidity < 0.0) {
            return BiomeType::Badlands;  // Add this to BiomeType enum
        }
        // Regular hills otherwise
        return BiomeType::Hills;  // Add this to BiomeType enum
    }
    
    // Temperature-humidity matrix for main biomes
    // This creates a 2D classification system that considers both factors
    
    // Hot climates (temp > 0.4)
    if (temp > 0.4) {
        if (humidity < -0.6) {
            return BiomeType::Desert;  // Very dry
        }
        if (humidity < -0.2) {
            return BiomeType::Savanna;  // Moderately dry
        }
        if (humidity < 0.2) {
            return BiomeType::TropicalShrubland;  // Add this to BiomeType enum
        }
        if (humidity < 0.6) {
            return BiomeType::TropicalSeasonalForest;  // Add this to BiomeType enum
        }
        return BiomeType::TropicalRainforest;  // Very wet
    }
    
    // Warm climates (temp > 0.2)
    if (temp > 0.2) {
        if (humidity < -0.4) {
            return BiomeType::Shrubland;  // Add this to BiomeType enum
        }
        if (humidity < 0.0) {
            return BiomeType::Plains;
        }
        if (humidity < 0.4) {
            return BiomeType::DeciduousForest;  // Add this to BiomeType enum
        }
        return BiomeType::TemperateRainforest;  // Add this to BiomeType enum
    }
    
    // Temperate climates (temp > 0.0)
    if (temp > 0.0) {
        if (humidity < -0.3) {
            return BiomeType::TemperateGrassland;  // Add this to BiomeType enum
        }
        if (humidity < 0.2) {
            return BiomeType::MixedForest;  // Add this to BiomeType enum
        }
        return BiomeType::Forest;
    }
    
    // Cool climates (temp > -0.3)
    if (temp > -0.3) {
        if (humidity < -0.2) {
            return BiomeType::BorealPlains;  // Add this to BiomeType enum
        }
        return BiomeType::Taiga;
    }
    
    // Cold climates (temp > -0.6)
    if (temp > -0.6) {
        if (humidity < 0.0) {
            return BiomeType::Tundra;
        }
        return BiomeType::SnowyTaiga;  // Add this to BiomeType enum
    }
    
    // Extremely cold climates
    return BiomeType::IceSheet;  // Add this to BiomeType enum
}
