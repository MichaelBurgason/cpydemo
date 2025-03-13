#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "CPYDemo.h"
#include "TerrainGen2DAnimation.h" // Include the terrain generator header

namespace py = pybind11;

PYBIND11_MODULE(noise_demo, m) {
    m.doc() = "Noisedemo Python binding with transfer function terrain generation";
    
    // Existing CPYDemo bindings
    py::class_<CPYDemo>(m, "CPYDemo")
        .def(py::init<>())
        .def("print_world", &CPYDemo::printWorld, "Prints 'Hello World' to console")
        .def("count", &CPYDemo::count, "Counts up to n and returns the sum", py::arg("n"))
        .def("print_hello", &CPYDemo::printHello, "Prints 'hello'")
        .def("test_function", &CPYDemo::testFunction, "returns the passed double", py::arg("n"))
        .def("generate_noise_map", &CPYDemo::generate_noise_map, "Generate a 2D noise map",
            py::arg("width"), py::arg("height"), py::arg("scale"), py::arg("offsetX"), py::arg("offsetY"));
    
    // Biome type enum
    py::enum_<BiomeType>(m, "BiomeType")
    // Water biomes
    .value("DeepOcean", BiomeType::DeepOcean)
    .value("Ocean", BiomeType::Ocean)
    .value("WarmOcean", BiomeType::WarmOcean)
    .value("FrozenOcean", BiomeType::FrozenOcean)
    
    // Coastal biomes
    .value("Beach", BiomeType::Beach)
    .value("RockyShore", BiomeType::RockyShore)
    .value("Marsh", BiomeType::Marsh)
    
    // Hot dry biomes
    .value("Desert", BiomeType::Desert)
    .value("Savanna", BiomeType::Savanna)
    .value("Shrubland", BiomeType::Shrubland)
    .value("TropicalShrubland", BiomeType::TropicalShrubland)
    
    // Grassland biomes
    .value("Plains", BiomeType::Plains)
    .value("TemperateGrassland", BiomeType::TemperateGrassland)
    .value("BorealPlains", BiomeType::BorealPlains)
    
    // Cold flat biomes
    .value("Tundra", BiomeType::Tundra)
    
    // Tropical forest biomes
    .value("TropicalRainforest", BiomeType::TropicalRainforest)
    .value("TropicalSeasonalForest", BiomeType::TropicalSeasonalForest)
    
    // Temperate forest biomes
    .value("TemperateRainforest", BiomeType::TemperateRainforest)
    .value("DeciduousForest", BiomeType::DeciduousForest)
    .value("Forest", BiomeType::Forest)
    .value("MixedForest", BiomeType::MixedForest)
    
    // Boreal forest biomes
    .value("Taiga", BiomeType::Taiga)
    .value("SnowyTaiga", BiomeType::SnowyTaiga)
    
    // Elevated terrain
    .value("Hills", BiomeType::Hills)
    .value("Badlands", BiomeType::Badlands)
    .value("Mountain", BiomeType::Mountain)
    .value("SnowyMountain", BiomeType::SnowyMountain)
    .value("SnowyCap", BiomeType::SnowyCap)
    
    // Ice biomes
    .value("IceSheet", BiomeType::IceSheet)
    
    .export_values();
    
    // TerrainParameter enum
    py::enum_<TerrainParameter>(m, "TerrainParameter")
        .value("BaseHeight", TerrainParameter::BaseHeight)
        .value("Continentalness", TerrainParameter::Continentalness)
        .value("RidgeFormation", TerrainParameter::RidgeFormation)
        .value("Erosion", TerrainParameter::Erosion)
        .value("Temperature", TerrainParameter::Temperature)
        .value("Humidity", TerrainParameter::Humidity)
        .export_values();

    // Legacy SplineControlPoint struct
    py::class_<SplineControlPoint>(m, "SplineControlPoint")
        .def(py::init<>())
        .def_readwrite("x", &SplineControlPoint::x)
        .def_readwrite("y", &SplineControlPoint::y)
        .def_readwrite("value", &SplineControlPoint::value)
        .def_readwrite("weight", &SplineControlPoint::weight)
        .def_readwrite("radius", &SplineControlPoint::radius);
        
    // Legacy TerrainSpline struct
    py::class_<TerrainSpline>(m, "TerrainSpline")
        .def(py::init<>())
        .def_readwrite("controlPoints", &TerrainSpline::controlPoints)
        .def_readwrite("tension", &TerrainSpline::tension)
        .def_readwrite("bias", &TerrainSpline::bias)
        .def_readwrite("continuity", &TerrainSpline::continuity)
        .def_readwrite("intensify", &TerrainSpline::intensify)
        .def_readwrite("influence", &TerrainSpline::influence)
        .def_readwrite("name", &TerrainSpline::name)
        .def_readwrite("parameter", &TerrainSpline::parameter);
    
    // New TransferControlPoint struct
    py::class_<TransferControlPoint>(m, "TransferControlPoint")
        .def(py::init<>())
        .def(py::init<double, double, double>(), 
            py::arg("inputValue"), py::arg("outputValue"), py::arg("weight") = 1.0)
        .def_readwrite("inputValue", &TransferControlPoint::inputValue)
        .def_readwrite("outputValue", &TransferControlPoint::outputValue)
        .def_readwrite("weight", &TransferControlPoint::weight);
        
    // New TransferSpline struct
    py::class_<TransferSpline>(m, "TransferSpline")
        .def(py::init<>())
        .def_readwrite("controlPoints", &TransferSpline::controlPoints)
        .def_readwrite("tension", &TransferSpline::tension)
        .def_readwrite("weight", &TransferSpline::weight)
        .def_readwrite("name", &TransferSpline::name)
        .def_readwrite("parameter", &TransferSpline::parameter);
    
    // TerrainParams struct with new fields
    py::class_<TerrainParams>(m, "TerrainParams")
        .def(py::init<>())
        // Base terrain parameters
        .def_readwrite("height_octaves", &TerrainParams::height_octaves)
        .def_readwrite("height_persistence", &TerrainParams::height_persistence)
        .def_readwrite("height_lacunarity", &TerrainParams::height_lacunarity)
        .def_readwrite("height_scale", &TerrainParams::height_scale)
        
        // Continentalness parameters
        .def_readwrite("continental_scale", &TerrainParams::continental_scale)
        .def_readwrite("continental_octaves", &TerrainParams::continental_octaves)
        
        // Ridge formation parameters
        .def_readwrite("ridge_scale", &TerrainParams::ridge_scale)
        .def_readwrite("ridge_octaves", &TerrainParams::ridge_octaves)
        .def_readwrite("ridge_persistence", &TerrainParams::ridge_persistence)
        .def_readwrite("ridge_influence", &TerrainParams::ridge_influence)
        
        // Erosion parameters
        .def_readwrite("erosion_scale", &TerrainParams::erosion_scale)
        .def_readwrite("erosion_octaves", &TerrainParams::erosion_octaves)
        
        // Temperature parameters
        .def_readwrite("temperature_scale", &TerrainParams::temperature_scale)
        .def_readwrite("temp_latitude_influence", &TerrainParams::temp_latitude_influence)
        
        // Humidity parameters
        .def_readwrite("humidity_scale", &TerrainParams::humidity_scale)
        .def_readwrite("humidity_octaves", &TerrainParams::humidity_octaves)
        
        // Peak intensity (legacy support)
        .def_readwrite("peaks_influence", &TerrainParams::peaks_influence)
        
        // Animation frames
        .def_readwrite("animation_frames", &TerrainParams::animation_frames)
        
        // Legacy auto-generation options
        .def_readwrite("auto_generate_splines", &TerrainParams::auto_generate_splines)
        .def_readwrite("parameterSplines", &TerrainParams::parameterSplines)
        .def_readwrite("splineInfluence", &TerrainParams::splineInfluence)
        .def_readwrite("auto_spline_count", &TerrainParams::auto_spline_count)
        
        // New transfer function options
        .def_readwrite("auto_generate_transfer_functions", &TerrainParams::auto_generate_transfer_functions)
        .def_readwrite("transferFunctions", &TerrainParams::transferFunctions);
    
    // TerrainCell struct with new fields
    py::class_<TerrainCell>(m, "TerrainCell")
        .def(py::init<>())
        // Base fields
        .def_readonly("baseNoise", &TerrainCell::baseNoise)
        .def_readonly("height", &TerrainCell::height)
        .def_readonly("biome", &TerrainCell::biome)
        
        // Raw parameter values
        .def_readonly("rawTemperature", &TerrainCell::rawTemperature)
        .def_readonly("rawHumidity", &TerrainCell::rawHumidity)
        .def_readonly("rawErosion", &TerrainCell::rawErosion)
        .def_readonly("rawContinentalness", &TerrainCell::rawContinentalness)
        .def_readonly("rawRidgeValue", &TerrainCell::rawRidgeValue)
        
        // Transformed parameter values
        .def_readonly("temperature", &TerrainCell::temperature)
        .def_readonly("humidity", &TerrainCell::humidity)
        .def_readonly("erosion", &TerrainCell::erosion)
        .def_readonly("continentalness", &TerrainCell::continentalness)
        .def_readonly("ridgeValue", &TerrainCell::ridgeValue);
    
    // TerrainGen2DAnimation class with updated methods
    py::class_<TerrainGen2DAnimation>(m, "TerrainGen2DAnimation")
        .def(py::init<int64_t>(), py::arg("seed") = 321314)
        .def("generateAnimatedTerrain", &TerrainGen2DAnimation::generateAnimatedTerrain,
            "Generate a 3D array of terrain data with animation frames",
            py::arg("width"), py::arg("height"), py::arg("params") = TerrainParams())
        .def("getAnimationFrameNames", &TerrainGen2DAnimation::getAnimationFrameNames,
            "Get names for each animation frame", py::arg("params") = TerrainParams())
        .def_static("classifyBiome", &TerrainGen2DAnimation::classifyBiome,
            "Classify a terrain cell into a biome type")
        .def_static("generateDefaultTransferFunction", &TerrainGen2DAnimation::generateDefaultTransferFunction,
            "Generate a default transfer function for a parameter", py::arg("parameter"));
    
    // Helper functions to extract NumPy arrays from the 3D terrain data
    m.def("extract_height_data", [](const std::vector<std::vector<std::vector<TerrainCell>>>& terrain) {
        if (terrain.empty() || terrain[0].empty() || terrain[0][0].empty()) {
            throw std::runtime_error("Terrain data is empty");
        }
        
        int height = terrain.size();
        int width = terrain[0].size();
        int frames = terrain[0][0].size();
        
        py::array_t<double> result({height, width, frames});
        auto buf = result.mutable_unchecked<3>();
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                for (int f = 0; f < frames; f++) {
                    buf(y, x, f) = terrain[y][x][f].height;
                }
            }
        }
        
        return result;
    }, "Extract height data as NumPy array from terrain data");
    
    m.def("extract_base_noise_data", [](const std::vector<std::vector<std::vector<TerrainCell>>>& terrain) {
        if (terrain.empty() || terrain[0].empty() || terrain[0][0].empty()) {
            throw std::runtime_error("Terrain data is empty");
        }
        
        int height = terrain.size();
        int width = terrain[0].size();
        int frames = terrain[0][0].size();
        
        py::array_t<double> result({height, width, frames});
        auto buf = result.mutable_unchecked<3>();
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                for (int f = 0; f < frames; f++) {
                    buf(y, x, f) = terrain[y][x][f].baseNoise;
                }
            }
        }
        
        return result;
    }, "Extract base noise data as NumPy array from terrain data");
    
    // Add functions to extract raw parameter values from terrain data
    m.def("extract_raw_temperature_data", [](const std::vector<std::vector<std::vector<TerrainCell>>>& terrain) {
        if (terrain.empty() || terrain[0].empty() || terrain[0][0].empty()) {
            throw std::runtime_error("Terrain data is empty");
        }
        
        int height = terrain.size();
        int width = terrain[0].size();
        int frames = terrain[0][0].size();
        
        py::array_t<double> result({height, width, frames});
        auto buf = result.mutable_unchecked<3>();
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                for (int f = 0; f < frames; f++) {
                    buf(y, x, f) = terrain[y][x][f].rawTemperature;
                }
            }
        }
        
        return result;
    }, "Extract raw temperature data as NumPy array from terrain data");
    
    m.def("extract_raw_humidity_data", [](const std::vector<std::vector<std::vector<TerrainCell>>>& terrain) {
        if (terrain.empty() || terrain[0].empty() || terrain[0][0].empty()) {
            throw std::runtime_error("Terrain data is empty");
        }
        
        int height = terrain.size();
        int width = terrain[0].size();
        int frames = terrain[0][0].size();
        
        py::array_t<double> result({height, width, frames});
        auto buf = result.mutable_unchecked<3>();
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                for (int f = 0; f < frames; f++) {
                    buf(y, x, f) = terrain[y][x][f].rawHumidity;
                }
            }
        }
        
        return result;
    }, "Extract raw humidity data as NumPy array from terrain data");
    
    m.def("extract_raw_erosion_data", [](const std::vector<std::vector<std::vector<TerrainCell>>>& terrain) {
        if (terrain.empty() || terrain[0].empty() || terrain[0][0].empty()) {
            throw std::runtime_error("Terrain data is empty");
        }
        
        int height = terrain.size();
        int width = terrain[0].size();
        int frames = terrain[0][0].size();
        
        py::array_t<double> result({height, width, frames});
        auto buf = result.mutable_unchecked<3>();
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                for (int f = 0; f < frames; f++) {
                    buf(y, x, f) = terrain[y][x][f].rawErosion;
                }
            }
        }
        
        return result;
    }, "Extract raw erosion data as NumPy array from terrain data");
    
    m.def("extract_raw_continentalness_data", [](const std::vector<std::vector<std::vector<TerrainCell>>>& terrain) {
        if (terrain.empty() || terrain[0].empty() || terrain[0][0].empty()) {
            throw std::runtime_error("Terrain data is empty");
        }
        
        int height = terrain.size();
        int width = terrain[0].size();
        int frames = terrain[0][0].size();
        
        py::array_t<double> result({height, width, frames});
        auto buf = result.mutable_unchecked<3>();
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                for (int f = 0; f < frames; f++) {
                    buf(y, x, f) = terrain[y][x][f].rawContinentalness;
                }
            }
        }
        
        return result;
    }, "Extract raw continentalness data as NumPy array from terrain data");
    
    m.def("extract_raw_ridge_value_data", [](const std::vector<std::vector<std::vector<TerrainCell>>>& terrain) {
        if (terrain.empty() || terrain[0].empty() || terrain[0][0].empty()) {
            throw std::runtime_error("Terrain data is empty");
        }
        
        int height = terrain.size();
        int width = terrain[0].size();
        int frames = terrain[0][0].size();
        
        py::array_t<double> result({height, width, frames});
        auto buf = result.mutable_unchecked<3>();
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                for (int f = 0; f < frames; f++) {
                    buf(y, x, f) = terrain[y][x][f].rawRidgeValue;
                }
            }
        }
        
        return result;
    }, "Extract raw ridge formation data as NumPy array from terrain data");
    
    // Keep the existing transformed data extraction functions
    m.def("extract_temperature_data", [](const std::vector<std::vector<std::vector<TerrainCell>>>& terrain) {
        if (terrain.empty() || terrain[0].empty() || terrain[0][0].empty()) {
            throw std::runtime_error("Terrain data is empty");
        }
        
        int height = terrain.size();
        int width = terrain[0].size();
        int frames = terrain[0][0].size();
        
        py::array_t<double> result({height, width, frames});
        auto buf = result.mutable_unchecked<3>();
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                for (int f = 0; f < frames; f++) {
                    buf(y, x, f) = terrain[y][x][f].temperature;
                }
            }
        }
        
        return result;
    }, "Extract temperature data as NumPy array from terrain data");
    
    m.def("extract_humidity_data", [](const std::vector<std::vector<std::vector<TerrainCell>>>& terrain) {
        if (terrain.empty() || terrain[0].empty() || terrain[0][0].empty()) {
            throw std::runtime_error("Terrain data is empty");
        }
        
        int height = terrain.size();
        int width = terrain[0].size();
        int frames = terrain[0][0].size();
        
        py::array_t<double> result({height, width, frames});
        auto buf = result.mutable_unchecked<3>();
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                for (int f = 0; f < frames; f++) {
                    buf(y, x, f) = terrain[y][x][f].humidity;
                }
            }
        }
        
        return result;
    }, "Extract humidity data as NumPy array from terrain data");
    
    m.def("extract_erosion_data", [](const std::vector<std::vector<std::vector<TerrainCell>>>& terrain) {
        if (terrain.empty() || terrain[0].empty() || terrain[0][0].empty()) {
            throw std::runtime_error("Terrain data is empty");
        }
        
        int height = terrain.size();
        int width = terrain[0].size();
        int frames = terrain[0][0].size();
        
        py::array_t<double> result({height, width, frames});
        auto buf = result.mutable_unchecked<3>();
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                for (int f = 0; f < frames; f++) {
                    buf(y, x, f) = terrain[y][x][f].erosion;
                }
            }
        }
        
        return result;
    }, "Extract erosion data as NumPy array from terrain data");
    
    m.def("extract_continentalness_data", [](const std::vector<std::vector<std::vector<TerrainCell>>>& terrain) {
        if (terrain.empty() || terrain[0].empty() || terrain[0][0].empty()) {
            throw std::runtime_error("Terrain data is empty");
        }
        
        int height = terrain.size();
        int width = terrain[0].size();
        int frames = terrain[0][0].size();
        
        py::array_t<double> result({height, width, frames});
        auto buf = result.mutable_unchecked<3>();
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                for (int f = 0; f < frames; f++) {
                    buf(y, x, f) = terrain[y][x][f].continentalness;
                }
            }
        }
        
        return result;
    }, "Extract continentalness data as NumPy array from terrain data");
    
    m.def("extract_ridge_value_data", [](const std::vector<std::vector<std::vector<TerrainCell>>>& terrain) {
        if (terrain.empty() || terrain[0].empty() || terrain[0][0].empty()) {
            throw std::runtime_error("Terrain data is empty");
        }
        
        int height = terrain.size();
        int width = terrain[0].size();
        int frames = terrain[0][0].size();
        
        py::array_t<double> result({height, width, frames});
        auto buf = result.mutable_unchecked<3>();
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                for (int f = 0; f < frames; f++) {
                    buf(y, x, f) = terrain[y][x][f].ridgeValue;
                }
            }
        }
        
        return result;
    }, "Extract ridge formation data as NumPy array from terrain data");
    
    m.def("extract_biome_data", [](const std::vector<std::vector<std::vector<TerrainCell>>>& terrain) {
        if (terrain.empty() || terrain[0].empty() || terrain[0][0].empty()) {
            throw std::runtime_error("Terrain data is empty");
        }
        
        int height = terrain.size();
        int width = terrain[0].size();
        int frames = terrain[0][0].size();
        
        py::array_t<int> result({height, width, frames});
        auto buf = result.mutable_unchecked<3>();
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                for (int f = 0; f < frames; f++) {
                    buf(y, x, f) = static_cast<int>(terrain[y][x][f].biome);
                }
            }
        }
        
        return result;
    }, "Extract biome data as NumPy array from terrain data");

    // Legacy helper functions
    m.def("create_terrain_spline", [](const std::string& name, TerrainParameter parameter, bool intensify, double tension, double influence) {
        TerrainSpline spline;
        spline.name = name;
        spline.parameter = parameter;
        spline.intensify = intensify;
        spline.tension = tension;
        spline.influence = influence;
        spline.bias = 0.0;
        spline.continuity = 0.0;
        return spline;
    }, "Create a new parameter spline (legacy)",
    py::arg("name"), py::arg("parameter"), py::arg("intensify") = true, py::arg("tension") = 0.5, py::arg("influence") = 1.0);
    
    m.def("add_control_point", [](TerrainSpline& spline, int x, int y, double value, double weight, double radius) {
        SplineControlPoint point;
        point.x = x;
        point.y = y;
        point.value = value;
        point.weight = weight;
        point.radius = radius;
        spline.controlPoints.push_back(point);
        return point;
    }, "Add a control point to a spline (legacy)",
    py::arg("spline"), py::arg("x"), py::arg("y"), py::arg("value"), 
    py::arg("weight") = 1.0, py::arg("radius") = 20.0);
    
    m.def("add_parameter_spline", [](TerrainParams& params, TerrainParameter parameter, const TerrainSpline& spline) {
        params.parameterSplines[parameter].push_back(spline);
    }, "Add a spline to a specific parameter (legacy)",
    py::arg("params"), py::arg("parameter"), py::arg("spline"));
    
    // New helper functions for transfer functions
    m.def("create_transfer_spline", [](const std::string& name, TerrainParameter parameter, double tension, double weight) {
        TransferSpline spline;
        spline.name = name;
        spline.parameter = parameter;
        spline.tension = tension;
        spline.weight = weight;
        return spline;
    }, "Create a new transfer function spline",
    py::arg("name"), py::arg("parameter"), py::arg("tension") = 0.5, py::arg("weight") = 1.0);
    
    m.def("add_transfer_control_point", [](TransferSpline& spline, double inputValue, double outputValue, double weight) {
        TransferControlPoint point;
        point.inputValue = inputValue;
        point.outputValue = outputValue;
        point.weight = weight;
        spline.controlPoints.push_back(point);
        return point;
    }, "Add a control point to a transfer function",
    py::arg("spline"), py::arg("inputValue"), py::arg("outputValue"), py::arg("weight") = 1.0);
    
    m.def("add_transfer_function", [](TerrainParams& params, TerrainParameter parameter, const TransferSpline& spline) {
        params.transferFunctions[parameter].push_back(spline);
    }, "Add a transfer function to a specific parameter",
    py::arg("params"), py::arg("parameter"), py::arg("spline"));
    
    // Helper to visualize a transfer function
    m.def("generate_transfer_curve", [](const TransferSpline& spline, int numPoints) {
        std::vector<std::pair<double, double>> curve;
        
        // Implement our own version of the transfer function evaluation
        // (similar to the implementation in TerrainGen2DAnimation)
        auto evaluateTransferCurve = [](const TransferSpline& spline, double inputValue) -> double {
            const auto& points = spline.controlPoints;
            
            // Need a minimum number of points for interpolation
            if (points.size() < 2) {
                return inputValue; // Just return the input value if insufficient control points
            }
            
            // Special case: input value is outside the range of control points
            if (inputValue < points.front().inputValue) {
                return points.front().outputValue;
            }
            if (inputValue > points.back().inputValue) {
                return points.back().outputValue;
            }
            
            // Find the segment that contains the input value
            size_t lowerIndex = 0;
            size_t upperIndex = 1;
            for (size_t i = 0; i < points.size() - 1; i++) {
                if (inputValue >= points[i].inputValue && inputValue <= points[i + 1].inputValue) {
                    lowerIndex = i;
                    upperIndex = i + 1;
                    break;
                }
            }
            
            // Calculate interpolation parameter t (0 to 1)
            double t = (inputValue - points[lowerIndex].inputValue) / 
                      (points[upperIndex].inputValue - points[lowerIndex].inputValue);
            
            // Simple linear interpolation as fallback
            if (spline.controlPoints.size() <= 2) {
                return points[lowerIndex].outputValue * (1.0 - t) + 
                       points[upperIndex].outputValue * t;
            }
            
            // For cubic interpolation, find surrounding points
            size_t p0Index = (lowerIndex > 0) ? lowerIndex - 1 : lowerIndex;
            size_t p1Index = lowerIndex;
            size_t p2Index = upperIndex;
            size_t p3Index = (upperIndex < points.size() - 1) ? upperIndex + 1 : upperIndex;
            
            // Get the output values for interpolation
            double p0 = points[p0Index].outputValue;
            double p1 = points[p1Index].outputValue;
            double p2 = points[p2Index].outputValue;
            double p3 = points[p3Index].outputValue;
            
            // Cardinal spline formula
            double tension = spline.tension;
            double t2 = t * t;
            double t3 = t2 * t;
            
            double s = (1.0 - tension) / 2.0;
            
            double h1 = 2*t3 - 3*t2 + 1;
            double h2 = -2*t3 + 3*t2;
            double h3 = t3 - 2*t2 + t;
            double h4 = t3 - t2;
            
            return h1*p1 + h2*p2 + s * (h3*(p2-p0) + h4*(p3-p1));
        };
        
        // Generate points along the curve
        for (int i = 0; i < numPoints; i++) {
            double inputValue = -1.0 + 2.0 * i / (numPoints - 1);
            double outputValue = evaluateTransferCurve(spline, inputValue);
            curve.push_back(std::make_pair(inputValue, outputValue));
        }
        
        return curve;
    }, "Generate points along a transfer function curve for visualization",
    py::arg("spline"), py::arg("numPoints") = 100);
}