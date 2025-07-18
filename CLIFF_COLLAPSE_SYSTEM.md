# CoastalME Cliff Collapse System Documentation

## Overview

The CoastalME cliff collapse system provides a pluggable architecture for simulating cliff erosion and collapse processes. This system allows different cliff collapse algorithms to be easily integrated and switched between during simulation runs.

## Architecture

The cliff collapse system is built on a modular, extensible architecture consisting of several key components:

### Core Components

#### 1. Abstract Base Class (`CCliffAlgorithm`)
- **Location**: `src/cliff_algorithm.h`, `src/cliff_algorithm.cpp`
- **Purpose**: Defines the interface that all cliff collapse algorithms must implement
- **Key Methods**:
  - `Initialize()`: Initialize the algorithm with configuration parameters
  - `ProcessTimestep()`: Process cliff collapse for a single timestep
  - `Reset()`: Reset algorithm state
  - `GetAlgorithmName()`: Get the algorithm's name
  - `GetAlgorithmDescription()`: Get the algorithm's description

#### 2. Data Structures

##### CCliffData
Encapsulates input data passed to cliff algorithms:
- `nCoast`: Coast index
- `nPointOnCoast`: Point index on the coast
- `dTimeStep`: Current timestep duration
- `CliffToe`: Cliff toe position (CGeom2DPoint)
- `CliffTop`: Cliff top position (CGeom2DPoint)
- `dStillWaterLevel`: Current still water level
- `dWaveHeight`: Wave height at cliff location
- `dWavePeriod`: Wave period
- `dWaveAngle`: Wave approach angle

##### CCliffResults
Encapsulates output data from cliff algorithms:
- `dErosionVolume`: Volume of cliff material eroded
- `bCollapseOccurred`: Whether collapse occurred this timestep
- `dCollapseVolume`: Volume of collapsed material
- `NewCliffToe`: Updated cliff toe position
- `NewCliffTop`: Updated cliff top position
- `mapSedimentDeposition`: Sediment deposition locations and volumes

#### 3. Algorithm Factory (`CCliffAlgorithmFactory`)
- **Location**: `src/cliff_algorithm_factory.h`, `src/cliff_algorithm_factory.cpp`
- **Purpose**: Creates and manages cliff algorithm instances
- **Key Methods**:
  - `CreateAlgorithm()`: Create algorithm instance by name
  - `GetAvailableAlgorithms()`: Get list of available algorithms
  - `IsValidAlgorithm()`: Check if algorithm name is valid

#### 4. Cliff Collapse Manager (`CCliffCollapseManager`)
- **Location**: `src/cliff_collapse_manager.h`, `src/cliff_collapse_manager.cpp`
- **Purpose**: Coordinates between CoastalME simulation and cliff algorithms
- **Key Methods**:
  - `Initialize()`: Initialize with selected algorithm
  - `ProcessAllCliffCollapse()`: Process all coasts and points
  - `ProcessCoastCliffCollapse()`: Process single coast
  - `ProcessPointCliffCollapse()`: Process single point
  - `ExtractCliffData()`: Extract simulation data for algorithms
  - `ApplyCliffResults()`: Apply algorithm results to simulation

## Available Algorithms

### 1. Simple Notch Algorithm (`simple_notch`)
- **Location**: `src/cliff_algorithm_simple_notch.h`, `src/cliff_algorithm_simple_notch.cpp`
- **Description**: A basic notch-based cliff erosion model
- **Features**:
  - Maintains notch depth state for each cliff point
  - Wave energy-based erosion calculation
  - Collapse occurs when notch depth exceeds threshold
  - Configurable erosion rate and collapse threshold

### 2. Legacy Algorithm (`legacy_original`)
- **Location**: `src/cliff_algorithm_legacy.h`, `src/cliff_algorithm_legacy.cpp`
- **Description**: Wrapper for the original CoastalME cliff collapse algorithm
- **Features**:
  - Maintains compatibility with existing simulations
  - Uses original cliff collapse logic from `do_cliff_collapse.cpp`

## Configuration

### Algorithm Selection

The cliff collapse algorithm is selected via the configuration file (`.dat` file) using case number 60:

```
# In your .dat configuration file
60 simple_notch
```

**Available Options**:
- `simple_notch`: Simple notch-based erosion model
- `legacy_original`: Original CoastalME cliff collapse algorithm

### Configuration File Integration

The algorithm selection is integrated into the standard CoastalME configuration system:

- **Case Number**: 60
- **Parameter**: Cliff algorithm name (string)
- **Validation**: Automatically validates against available algorithms
- **Error Handling**: Provides clear error messages with available options

## Performance Optimization

### OpenMP Parallelization

The cliff collapse system includes OpenMP optimization for high-performance processing:

#### Coast-Level Parallelization
```cpp
#ifdef _OPENMP
#pragma omp parallel for
#endif
for (int nCoast = 0; nCoast < nCoasts; nCoast++)
{
    ProcessCoastCliffCollapse(nCoast);
}
```

#### Point-Level Parallelization
```cpp
#ifdef _OPENMP
#pragma omp parallel for
#endif
for (int nPoint = 0; nPoint < nCoastPoints; nPoint++)
{
    ProcessPointCliffCollapse(nCoast, nPoint);
}
```

### Performance Benefits
- **Nested Parallelization**: Utilizes multiple CPU cores at both coast and point levels
- **Scalability**: Performance scales with number of CPU cores
- **Thread Safety**: Each algorithm instance maintains independent state
- **Memory Efficiency**: Minimal memory overhead per thread

## Integration with CoastalME

### Simulation Integration

The cliff collapse system integrates seamlessly with the CoastalME simulation:

1. **Initialization**: `CCliffCollapseManager` is created during simulation startup
2. **Algorithm Selection**: Algorithm is loaded based on configuration file
3. **Timestep Processing**: Called during each simulation timestep
4. **Result Application**: Results are applied to simulation state

### Friend Class Access

The system uses friend class declarations to access private simulation data:
```cpp
class CSimulation
{
    friend class CCliffCollapseManager;
    // ...
};
```

### Memory Management

- **Smart Pointers**: Uses `std::unique_ptr` for automatic memory management
- **RAII**: Resource acquisition is initialization for safe cleanup
- **Exception Safety**: Proper exception handling throughout

## Extending the System

### Adding New Algorithms

To add a new cliff collapse algorithm:

1. **Create Algorithm Class**:
   ```cpp
   class CMyCliffAlgorithm : public CCliffAlgorithm
   {
   public:
       bool Initialize(const string& strAlgorithmName) override;
       CCliffResults ProcessTimestep(const CCliffData& cliffData) override;
       void Reset() override;
       string GetAlgorithmName() const override;
       string GetAlgorithmDescription() const override;
   };
   ```

2. **Register in Factory**:
   ```cpp
   // In cliff_algorithm_factory.cpp
   if (strAlgorithmName == "my_algorithm")
       return std::make_unique<CMyCliffAlgorithm>();
   ```

3. **Add to Available Algorithms**:
   ```cpp
   // In GetAvailableAlgorithms()
   vecAlgorithms.push_back("my_algorithm");
   ```

### Algorithm Development Guidelines

- **State Management**: Maintain algorithm state using member variables
- **Thread Safety**: Ensure algorithm instances are thread-safe
- **Error Handling**: Use CoastalME error return codes (RTN_OK, RTN_ERR_*)
- **Documentation**: Provide clear algorithm descriptions
- **Testing**: Test with different cliff configurations

## File Structure

```
src/
├── cliff_algorithm.h/.cpp              # Abstract base class
├── cliff_algorithm_factory.h/.cpp      # Factory pattern implementation
├── cliff_algorithm_simple_notch.h/.cpp # Simple notch algorithm
├── cliff_algorithm_legacy.h/.cpp       # Legacy algorithm wrapper
├── cliff_collapse_manager.h/.cpp       # Main coordination class
├── cliff_wrapper.h/.cpp               # Legacy wrapper compatibility
└── simulation.h/.cpp                  # Integration with main simulation
```

## Error Handling

The system provides comprehensive error handling:

### Configuration Errors
- **Invalid Algorithm**: Clear error message with available options
- **Missing Configuration**: Default to simple_notch algorithm
- **File Access**: Standard CoastalME file error handling

### Runtime Errors
- **Algorithm Failure**: Graceful degradation to no-collapse mode
- **Data Validation**: Input data validation before processing
- **Memory Allocation**: Safe memory management with smart pointers

### Error Messages
```cpp
// Example error message for invalid algorithm
cerr << ERR << "unknown cliff algorithm '" << strAlgorithmName 
     << "'. Available algorithms: simple_notch, legacy_original" << endl;
```

## Future Enhancements

### Planned Features
1. **External Algorithm Loading**: Dynamic loading of algorithm plugins
2. **Algorithm Parameters**: Configuration-based algorithm parameters
3. **Multi-Algorithm Support**: Different algorithms for different coast segments
4. **Visualization**: Real-time visualization of cliff collapse processes
5. **Calibration Tools**: Automated parameter calibration utilities

### Research Integration
- **Coupling with External Models**: Interface for external cliff models
- **Machine Learning**: Integration with ML-based cliff prediction models
- **Uncertainty Quantification**: Monte Carlo simulation support
- **Validation Framework**: Automated validation against field data

## Technical Requirements

### Build Dependencies
- **C++17 Standard**: Modern C++ features
- **OpenMP**: Parallel processing support
- **CMake**: Build system integration
- **CoastalME Core**: Integration with main simulation

### Runtime Requirements
- **Memory**: Minimal additional memory overhead
- **CPU**: Multi-core processors recommended for performance
- **Storage**: No additional storage requirements

## Testing

### Test Coverage
- **Unit Tests**: Individual algorithm testing
- **Integration Tests**: Full simulation integration
- **Performance Tests**: OpenMP scaling verification
- **Validation Tests**: Comparison with reference results

### Test Data
- **Synthetic Cases**: Simple test cases for algorithm verification
- **Real-World Cases**: Happisburgh and other field sites
- **Performance Benchmarks**: Timing and scaling tests

## Support and Maintenance

### Version Compatibility
- **Backward Compatibility**: Legacy algorithm maintains old behavior
- **Configuration Migration**: Automatic migration tools
- **API Stability**: Stable interface for algorithm development

### Documentation Updates
- **Algorithm Documentation**: Each algorithm includes detailed documentation
- **API Reference**: Complete API documentation with examples
- **User Guide**: Step-by-step configuration and usage guide

For questions or issues with the cliff collapse system, please refer to the main CoastalME documentation or contact the development team.