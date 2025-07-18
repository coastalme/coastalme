# CoastalME Cliff Collapse System - API Reference

## Overview

This document provides a comprehensive API reference for the CoastalME cliff collapse system, intended for developers who want to extend the system with new algorithms or integrate with external models.

## Core Classes

### CCliffAlgorithm (Abstract Base Class)

**Header**: `cliff_algorithm.h`

The abstract base class that defines the interface for all cliff collapse algorithms.

#### Pure Virtual Methods

```cpp
virtual bool Initialize(const string& strAlgorithmName) = 0;
```
- **Purpose**: Initialize the algorithm with the given name
- **Parameters**: `strAlgorithmName` - Name of the algorithm
- **Returns**: `true` if initialization successful, `false` otherwise
- **Notes**: Called once during simulation setup

```cpp
virtual CCliffResults ProcessTimestep(const CCliffData& cliffData) = 0;
```
- **Purpose**: Process cliff collapse for a single timestep
- **Parameters**: `cliffData` - Input data for the cliff point
- **Returns**: `CCliffResults` containing erosion/collapse results
- **Notes**: Called every timestep for each cliff point

```cpp
virtual void Reset() = 0;
```
- **Purpose**: Reset algorithm state to initial conditions
- **Parameters**: None
- **Returns**: void
- **Notes**: Called when simulation is reset

```cpp
virtual string GetAlgorithmName() const = 0;
```
- **Purpose**: Get the algorithm's name
- **Parameters**: None
- **Returns**: Algorithm name as string
- **Notes**: Used for identification and validation

```cpp
virtual string GetAlgorithmDescription() const = 0;
```
- **Purpose**: Get the algorithm's description
- **Parameters**: None
- **Returns**: Algorithm description as string
- **Notes**: Used for documentation and user information

#### Virtual Destructor

```cpp
virtual ~CCliffAlgorithm() = default;
```
- **Purpose**: Virtual destructor for proper cleanup
- **Notes**: Ensures derived classes are properly destroyed

### CCliffData (Input Data Structure)

**Header**: `cliff_algorithm.h`

Encapsulates input data passed to cliff algorithms.

#### Public Members

```cpp
int nCoast;
```
- **Purpose**: Coast index in the simulation
- **Range**: 0 to number of coasts - 1

```cpp
int nPointOnCoast;
```
- **Purpose**: Point index on the coast
- **Range**: 0 to number of points on coast - 1

```cpp
double dTimeStep;
```
- **Purpose**: Current timestep duration
- **Units**: Hours
- **Range**: Positive values

```cpp
CGeom2DPoint CliffToe;
```
- **Purpose**: Current cliff toe position
- **Units**: Meters (x, y coordinates)
- **Notes**: May be (0,0) if no cliff present

```cpp
CGeom2DPoint CliffTop;
```
- **Purpose**: Current cliff top position
- **Units**: Meters (x, y coordinates)
- **Notes**: May be (0,0) if no cliff present

```cpp
double dStillWaterLevel;
```
- **Purpose**: Current still water level
- **Units**: Meters above datum
- **Notes**: Includes tidal elevation

```cpp
double dWaveHeight;
```
- **Purpose**: Significant wave height at cliff location
- **Units**: Meters
- **Range**: Non-negative values

```cpp
double dWavePeriod;
```
- **Purpose**: Wave period
- **Units**: Seconds
- **Range**: Positive values

```cpp
double dWaveAngle;
```
- **Purpose**: Wave approach angle
- **Units**: Degrees
- **Range**: 0-360 degrees

#### Constructor

```cpp
CCliffData();
```
- **Purpose**: Default constructor
- **Initializes**: All numeric values to 0, points to origin

### CCliffResults (Output Data Structure)

**Header**: `cliff_algorithm.h`

Encapsulates output data from cliff algorithms.

#### Public Members

```cpp
double dErosionVolume;
```
- **Purpose**: Volume of cliff material eroded this timestep
- **Units**: Cubic meters
- **Range**: Non-negative values

```cpp
bool bCollapseOccurred;
```
- **Purpose**: Whether cliff collapse occurred this timestep
- **Values**: `true` if collapse occurred, `false` otherwise

```cpp
double dCollapseVolume;
```
- **Purpose**: Volume of material collapsed this timestep
- **Units**: Cubic meters
- **Range**: Non-negative values
- **Notes**: Should be 0 if bCollapseOccurred is false

```cpp
CGeom2DPoint NewCliffToe;
```
- **Purpose**: Updated cliff toe position after processing
- **Units**: Meters (x, y coordinates)
- **Notes**: May be same as input if no change

```cpp
CGeom2DPoint NewCliffTop;
```
- **Purpose**: Updated cliff top position after processing
- **Units**: Meters (x, y coordinates)
- **Notes**: May be same as input if no change

```cpp
map<string, double> mapSedimentDeposition;
```
- **Purpose**: Sediment deposition locations and volumes
- **Key**: Location identifier (e.g., "x123_y456")
- **Value**: Sediment volume in cubic meters
- **Notes**: Used for distributing collapsed material

#### Constructor

```cpp
CCliffResults();
```
- **Purpose**: Default constructor
- **Initializes**: All numeric values to 0, boolean to false, map empty

### CCliffAlgorithmFactory (Factory Class)

**Header**: `cliff_algorithm_factory.h`

Manages creation and registration of cliff algorithms.

#### Static Methods

```cpp
static unique_ptr<CCliffAlgorithm> CreateAlgorithm(const string& strAlgorithmName);
```
- **Purpose**: Create algorithm instance by name
- **Parameters**: `strAlgorithmName` - Name of algorithm to create
- **Returns**: Unique pointer to algorithm instance, nullptr if invalid name
- **Thread Safety**: Thread-safe

```cpp
static vector<string> GetAvailableAlgorithms();
```
- **Purpose**: Get list of available algorithm names
- **Parameters**: None
- **Returns**: Vector of available algorithm names
- **Thread Safety**: Thread-safe

```cpp
static bool IsValidAlgorithm(const string& strAlgorithmName);
```
- **Purpose**: Check if algorithm name is valid
- **Parameters**: `strAlgorithmName` - Name to validate
- **Returns**: `true` if valid, `false` otherwise
- **Thread Safety**: Thread-safe

#### Enum Classes

```cpp
enum class ECliffAlgorithm {
    SIMPLE_NOTCH,
    LEGACY_ORIGINAL
};
```
- **Purpose**: Enumeration of available algorithm types
- **Usage**: Internal factory management

### CCliffCollapseManager (Coordination Class)

**Header**: `cliff_collapse_manager.h`

Coordinates between the CoastalME simulation and cliff algorithms.

#### Public Methods

```cpp
CCliffCollapseManager(CSimulation* pSimulation);
```
- **Purpose**: Constructor
- **Parameters**: `pSimulation` - Pointer to main simulation object
- **Notes**: Stores simulation pointer for data access

```cpp
bool Initialize(const string& strAlgorithmName);
```
- **Purpose**: Initialize with selected algorithm
- **Parameters**: `strAlgorithmName` - Name of algorithm to use
- **Returns**: `true` if successful, `false` otherwise
- **Notes**: Creates and initializes algorithm instance

```cpp
int ProcessAllCliffCollapse();
```
- **Purpose**: Process cliff collapse for all coasts and points
- **Parameters**: None
- **Returns**: CoastalME return code (RTN_OK, RTN_ERR_*)
- **Notes**: OpenMP parallelized for performance

```cpp
int ProcessCoastCliffCollapse(int nCoast);
```
- **Purpose**: Process cliff collapse for a single coast
- **Parameters**: `nCoast` - Coast index
- **Returns**: CoastalME return code (RTN_OK, RTN_ERR_*)
- **Notes**: Processes all points on the coast

```cpp
int ProcessPointCliffCollapse(int nCoast, int nPointOnCoast);
```
- **Purpose**: Process cliff collapse for a single point
- **Parameters**: `nCoast` - Coast index, `nPointOnCoast` - Point index
- **Returns**: CoastalME return code (RTN_OK, RTN_ERR_*)
- **Notes**: Core processing for individual cliff points

```cpp
void Reset();
```
- **Purpose**: Reset algorithm state
- **Parameters**: None
- **Returns**: void
- **Notes**: Delegates to algorithm's Reset() method

```cpp
void SetCliffCollapseEnabled(bool bEnabled);
```
- **Purpose**: Enable or disable cliff collapse processing
- **Parameters**: `bEnabled` - Enable flag
- **Returns**: void

```cpp
bool IsCliffCollapseEnabled() const;
```
- **Purpose**: Check if cliff collapse is enabled
- **Parameters**: None
- **Returns**: `true` if enabled, `false` otherwise

```cpp
string GetAlgorithmName() const;
```
- **Purpose**: Get current algorithm name
- **Parameters**: None
- **Returns**: Algorithm name as string

```cpp
string GetAlgorithmDescription() const;
```
- **Purpose**: Get current algorithm description
- **Parameters**: None
- **Returns**: Algorithm description as string

```cpp
vector<string> GetAvailableAlgorithms() const;
```
- **Purpose**: Get list of available algorithms
- **Parameters**: None
- **Returns**: Vector of algorithm names

```cpp
bool SwitchAlgorithm(const string& strNewAlgorithmName);
```
- **Purpose**: Switch to a different algorithm
- **Parameters**: `strNewAlgorithmName` - Name of new algorithm
- **Returns**: `true` if successful, `false` otherwise
- **Notes**: Replaces current algorithm instance

#### Private Methods

```cpp
CCliffData ExtractCliffData(int nCoast, int nPointOnCoast) const;
```
- **Purpose**: Extract cliff data from simulation
- **Parameters**: `nCoast` - Coast index, `nPointOnCoast` - Point index
- **Returns**: `CCliffData` structure with extracted data
- **Notes**: Accesses simulation data via friend class

```cpp
void ApplyCliffResults(const CCliffResults& results, int nCoast, int nPointOnCoast);
```
- **Purpose**: Apply algorithm results to simulation
- **Parameters**: `results` - Algorithm results, `nCoast` - Coast index, `nPointOnCoast` - Point index
- **Returns**: void
- **Notes**: Updates simulation state with results

```cpp
void ApplyCollapseDeposition(const CCliffResults& results, int nCoast, int nPointOnCoast);
```
- **Purpose**: Apply sediment deposition from collapse
- **Parameters**: `results` - Algorithm results, `nCoast` - Coast index, `nPointOnCoast` - Point index
- **Returns**: void
- **Notes**: Distributes collapsed sediment

```cpp
void UpdateCliffGeometry(const CCliffResults& results, int nCoast, int nPointOnCoast);
```
- **Purpose**: Update cliff geometry after erosion/collapse
- **Parameters**: `results` - Algorithm results, `nCoast` - Coast index, `nPointOnCoast` - Point index
- **Returns**: void
- **Notes**: Updates cliff toe and top positions

## Implemented Algorithms

### CCliffAlgorithmSimpleNotch

**Header**: `cliff_algorithm_simple_notch.h`

A simple notch-based cliff erosion algorithm.

#### Key Features

- Maintains notch depth state for each cliff point
- Wave energy-based erosion calculation
- Collapse occurs when notch depth exceeds threshold

#### State Management

```cpp
map<string, double> m_mapNotchDepth;
```
- **Purpose**: Tracks notch depth for each cliff point
- **Key**: Point identifier (coast_point format)
- **Value**: Current notch depth in meters

#### Algorithm Parameters

```cpp
static const double EROSION_RATE_FACTOR = 0.001;
static const double COLLAPSE_THRESHOLD = 2.0;
```
- **EROSION_RATE_FACTOR**: Controls erosion rate (m/timestep per J/m²)
- **COLLAPSE_THRESHOLD**: Notch depth threshold for collapse (m)

### CCliffAlgorithmLegacy

**Header**: `cliff_algorithm_legacy.h`

Wrapper for the original CoastalME cliff collapse algorithm.

#### Key Features

- Maintains compatibility with existing simulations
- Uses original cliff collapse logic
- Identical results to previous versions

## Development Guidelines

### Creating New Algorithms

#### Step 1: Implement Algorithm Class

```cpp
#include "cliff_algorithm.h"

class CMyCliffAlgorithm : public CCliffAlgorithm
{
private:
    // Algorithm state variables
    map<string, double> m_mapState;
    
public:
    bool Initialize(const string& strAlgorithmName) override;
    CCliffResults ProcessTimestep(const CCliffData& cliffData) override;
    void Reset() override;
    string GetAlgorithmName() const override;
    string GetAlgorithmDescription() const override;
};
```

#### Step 2: Implement Required Methods

```cpp
bool CMyCliffAlgorithm::Initialize(const string& strAlgorithmName)
{
    // Initialize algorithm state
    m_mapState.clear();
    
    // Validate algorithm name
    if (strAlgorithmName != "my_algorithm")
        return false;
    
    return true;
}

CCliffResults CMyCliffAlgorithm::ProcessTimestep(const CCliffData& cliffData)
{
    CCliffResults results;
    
    // Create unique identifier for this point
    string strPointID = to_string(cliffData.nCoast) + "_" + to_string(cliffData.nPointOnCoast);
    
    // Algorithm-specific processing
    // ...
    
    return results;
}

void CMyCliffAlgorithm::Reset()
{
    m_mapState.clear();
}

string CMyCliffAlgorithm::GetAlgorithmName() const
{
    return "my_algorithm";
}

string CMyCliffAlgorithm::GetAlgorithmDescription() const
{
    return "My custom cliff erosion algorithm";
}
```

#### Step 3: Register in Factory

```cpp
// In cliff_algorithm_factory.cpp
unique_ptr<CCliffAlgorithm> CCliffAlgorithmFactory::CreateAlgorithm(const string& strAlgorithmName)
{
    if (strAlgorithmName == "simple_notch")
        return make_unique<CCliffAlgorithmSimpleNotch>();
    else if (strAlgorithmName == "legacy_original")
        return make_unique<CCliffAlgorithmLegacy>();
    else if (strAlgorithmName == "my_algorithm")  // Add this line
        return make_unique<CMyCliffAlgorithm>();
    
    return nullptr;
}

vector<string> CCliffAlgorithmFactory::GetAvailableAlgorithms()
{
    vector<string> vecAlgorithms;
    vecAlgorithms.push_back("simple_notch");
    vecAlgorithms.push_back("legacy_original");
    vecAlgorithms.push_back("my_algorithm");  // Add this line
    
    return vecAlgorithms;
}
```

### Thread Safety Guidelines

- **State Management**: Use point-specific keys for state maps
- **Shared Data**: Avoid shared mutable state between threads
- **Synchronization**: Use thread-local storage when necessary
- **Testing**: Test algorithms with OpenMP enabled

### Performance Guidelines

- **Memory Allocation**: Minimize dynamic memory allocation in ProcessTimestep
- **Computational Complexity**: Keep algorithms efficient for large datasets
- **Vectorization**: Use vectorizable operations where possible
- **Profiling**: Profile algorithms for performance bottlenecks

### Error Handling

- **Input Validation**: Validate input data in ProcessTimestep
- **Return Codes**: Use appropriate CoastalME return codes
- **Exception Safety**: Ensure algorithms don't throw exceptions
- **Graceful Degradation**: Handle invalid data gracefully

## Integration Points

### Simulation Integration

The cliff collapse system integrates with CoastalME at several points:

1. **Initialization**: `CSimulation::CSimulation()` creates manager
2. **Configuration**: `CSimulation::bReadIniFile()` reads algorithm selection
3. **Processing**: `CSimulation::nDoAllWaveEnergyToCoastLandforms()` calls manager
4. **Cleanup**: Automatic cleanup via smart pointers

### Data Access

The manager accesses simulation data through friend class access:

```cpp
class CSimulation
{
    friend class CCliffCollapseManager;
    
    // Manager can access private members:
    vector<CRWCoast> m_VCoast;
    // ... other private members
};
```

### Memory Management

- **Smart Pointers**: Use `std::unique_ptr` for automatic cleanup
- **RAII**: Resource acquisition is initialization
- **Exception Safety**: Proper cleanup on exceptions

## Testing Framework

### Unit Testing

```cpp
// Example unit test structure
void TestMyCliffAlgorithm()
{
    // Create algorithm instance
    auto pAlgorithm = make_unique<CMyCliffAlgorithm>();
    
    // Test initialization
    assert(pAlgorithm->Initialize("my_algorithm"));
    
    // Test processing
    CCliffData testData;
    testData.nCoast = 0;
    testData.nPointOnCoast = 0;
    testData.dTimeStep = 1.0;
    
    CCliffResults results = pAlgorithm->ProcessTimestep(testData);
    
    // Validate results
    assert(results.dErosionVolume >= 0.0);
    
    // Test reset
    pAlgorithm->Reset();
}
```

### Integration Testing

Test algorithms within the full simulation framework:

1. **Small Test Cases**: Simple cliff configurations
2. **Performance Tests**: Large-scale simulations
3. **Validation Tests**: Comparison with reference results
4. **Regression Tests**: Ensure no breaking changes

## Common Patterns

### State Management

```cpp
// Pattern for maintaining per-point state
map<string, StateData> m_mapPointState;

string CreatePointKey(int nCoast, int nPoint)
{
    return to_string(nCoast) + "_" + to_string(nPoint);
}
```

### Wave Energy Calculation

```cpp
double CalculateWaveEnergy(const CCliffData& data)
{
    // Common wave energy calculation
    double dWaveEnergy = 0.125 * WATER_DENSITY * GRAVITY * 
                        data.dWaveHeight * data.dWaveHeight * data.dWavePeriod;
    return dWaveEnergy;
}
```

### Collapse Threshold

```cpp
bool CheckCollapseThreshold(double dCurrentValue, double dThreshold)
{
    return dCurrentValue >= dThreshold;
}
```

## Error Codes

The system uses standard CoastalME error codes:

- `RTN_OK`: Success
- `RTN_ERR_NOCOAST`: Invalid coast index
- `RTN_ERR_BADPARAM`: Invalid parameters
- `RTN_ERR_MEMALLOC`: Memory allocation failure

## Future API Extensions

### Planned Enhancements

1. **Algorithm Parameters**: Configuration-based parameters
2. **External Models**: Interface for external cliff models
3. **Visualization**: Real-time visualization support
4. **Calibration**: Automated parameter calibration

### API Stability

- **Backward Compatibility**: Existing APIs will remain stable
- **Extension Points**: New features will extend existing interfaces
- **Deprecation Policy**: Deprecated features will be marked clearly

---

*This API reference provides comprehensive information for developing with the CoastalME cliff collapse system. For usage information, see the user guide documentation.*