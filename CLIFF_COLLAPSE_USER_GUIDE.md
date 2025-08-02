# CoastalME Cliff Collapse System - User Guide

## Quick Start

### 1. Enable Cliff Collapse in Configuration

To use the new cliff collapse system, add the following line to your `.dat` configuration file:

```
60 simple_notch
```

This line should be placed in the main configuration section of your file.

### 2. Prerequisites

Ensure your simulation has the required components:
- **Consolidated Sediment**: Cliff collapse requires consolidated sediment layers
- **Cliff Collapse Enabled**: The main cliff collapse option must be enabled
- **Coast Definition**: Proper coastline definition with cliff segments

### 3. Basic Configuration Example

```
# Example .dat file section
# ... other configuration parameters ...

# Enable cliff collapse
# (existing CoastalME parameter)
1 1  # Enable cliff collapse processing

# Select cliff collapse algorithm (NEW)
60 simple_notch

# ... rest of configuration ...
```

## Algorithm Selection

### Available Algorithms

#### Simple Notch Algorithm (`simple_notch`)
- **Best for**: General-purpose cliff erosion modeling
- **Features**: 
  - Wave energy-based erosion
  - Notch depth tracking
  - Collapse threshold triggering
- **Computational cost**: Low
- **Use cases**: Long-term simulations, parametric studies

#### Legacy Algorithm (`legacy_original`)
- **Best for**: Maintaining compatibility with existing simulations
- **Features**:
  - Original CoastalME cliff collapse behavior
  - Identical results to previous versions
- **Computational cost**: Low
- **Use cases**: Validation runs, backward compatibility

### Switching Algorithms

To switch between algorithms, simply change the algorithm name in line 60:

```
# For simple notch algorithm
60 simple_notch

# For legacy algorithm
60 legacy_original
```

## Configuration Details

### Configuration File Integration

The cliff collapse algorithm selection is integrated into the standard CoastalME configuration system:

- **Line number**: 60
- **Format**: `60 <algorithm_name>`
- **Validation**: Automatic validation against available algorithms
- **Case sensitivity**: Algorithm names are case-sensitive

### Error Handling

If an invalid algorithm name is provided, CoastalME will:
1. Display an error message listing available algorithms
2. Terminate the simulation
3. Provide guidance on correct algorithm names

Example error message:
```
ERROR: unknown cliff algorithm 'invalid_name'. Available algorithms: simple_notch, legacy_original
```

## Performance Optimization

### OpenMP Support

The cliff collapse system automatically utilizes OpenMP for parallel processing when available:

- **Automatic detection**: System detects OpenMP availability at compile time
- **Multi-core scaling**: Performance scales with number of CPU cores
- **No configuration required**: Optimization is automatic

### Performance Tips

1. **Multi-core Systems**: Use systems with multiple CPU cores for best performance
2. **Memory**: Ensure adequate RAM for large simulations
3. **Compiler Optimization**: Use optimized builds (Release mode)

### Monitoring Performance

To monitor cliff collapse performance:

```bash
# Time the simulation
time ./cme your_config.dat

# Monitor CPU usage
top -p $(pgrep cme)
```

## Algorithm-Specific Configuration

### Simple Notch Algorithm

Currently uses hardcoded parameters, but future versions will support:
- Erosion rate parameters
- Collapse threshold settings
- Wave energy coefficients

### Legacy Algorithm

Uses the same parameters as the original CoastalME cliff collapse system:
- All existing cliff collapse parameters apply
- No additional configuration required

## Common Use Cases

### 1. New Simulations

For new simulations, recommend using the simple notch algorithm:
```
60 simple_notch
```

### 2. Validation Studies

For validating against previous results, use the legacy algorithm:
```
60 legacy_original
```

### 3. Comparative Studies

Run the same simulation with different algorithms to compare results:
```
# Run 1: Simple notch
60 simple_notch

# Run 2: Legacy
60 legacy_original
```

## Troubleshooting

### Common Issues

#### 1. Algorithm Not Found
**Error**: `unknown cliff algorithm 'algorithm_name'`
**Solution**: Check spelling and use one of: `simple_notch`, `legacy_original`

#### 2. Configuration Line Missing
**Error**: Cliff collapse enabled but no algorithm specified
**Solution**: Add line `60 simple_notch` to your configuration file

#### 3. Prerequisites Not Met
**Error**: Cliff collapse algorithm requires consolidated sediment
**Solution**: Ensure consolidated sediment is enabled in your configuration

#### 4. Performance Issues
**Problem**: Slow cliff collapse processing
**Solutions**:
- Verify OpenMP is available and enabled
- Check system resources (CPU, memory)
- Consider using a smaller grid or fewer timesteps for testing

### Debug Mode

For debugging cliff collapse issues:

1. **Enable verbose output**: Check simulation log files
2. **Reduce complexity**: Start with simple test cases
3. **Monitor resources**: Check CPU and memory usage
4. **Validate input data**: Ensure cliff segments are properly defined

## Output and Results

### Standard Output

The cliff collapse system produces the same output files as the original system:
- `cliff_collapse_erosion.csv`
- `cliff_collapse_deposition.csv`
- `cliff_collapse_net.csv`
- Various raster outputs for visualization

### Algorithm-Specific Output

Each algorithm may produce additional information:
- **Simple Notch**: Notch depth information (planned)
- **Legacy**: Standard CoastalME cliff output

### Visualization

Results can be visualized using:
- **QGIS**: Load raster outputs for spatial visualization
- **Excel/CSV**: Analyze time series data
- **Python/R**: Custom analysis scripts

## Best Practices

### 1. Configuration Management

- **Version control**: Keep configuration files in version control
- **Documentation**: Document algorithm choices and reasoning
- **Testing**: Test configurations with small runs first

### 2. Performance

- **Profiling**: Profile performance for large simulations
- **Optimization**: Use Release builds for production runs
- **Monitoring**: Monitor system resources during runs

### 3. Validation

- **Comparison**: Compare results between algorithms
- **Sensitivity**: Test sensitivity to algorithm choice
- **Calibration**: Calibrate against field data when available

## Example Workflows

### Workflow 1: New Study Setup

1. **Start with simple notch**:
   ```
   60 simple_notch
   ```

2. **Run short test simulation**:
   ```bash
   ./cme test_config.dat
   ```

3. **Check results and performance**

4. **Scale up to full simulation**

### Workflow 2: Validation Study

1. **Run with legacy algorithm**:
   ```
   60 legacy_original
   ```

2. **Compare with historical results**

3. **Run with new algorithm**:
   ```
   60 simple_notch
   ```

4. **Compare and analyze differences**

### Workflow 3: Sensitivity Analysis

1. **Prepare base configuration**

2. **Run with different algorithms**:
   ```bash
   # Script to run multiple algorithms
   for alg in simple_notch legacy_original; do
       sed "s/^60 .*/60 $alg/" base_config.dat > config_$alg.dat
       ./cme config_$alg.dat
       mv output/ output_$alg/
   done
   ```

3. **Compare results**

## Migration from Legacy System

### Existing Simulations

To migrate existing simulations to the new system:

1. **Add algorithm selection line**:
   ```
   60 legacy_original
   ```

2. **Verify identical results**

3. **Optionally test new algorithms**:
   ```
   60 simple_notch
   ```

### Configuration Updates

- **No changes required** for existing parameters
- **Add line 60** for algorithm selection
- **All case numbers** remain the same (60 is new)

## Future Features

### Planned Enhancements

1. **Algorithm Parameters**: Configuration-based parameters for each algorithm
2. **External Algorithms**: Support for user-defined algorithms
3. **Real-time Visualization**: Live visualization of cliff collapse
4. **Calibration Tools**: Automated parameter calibration

### Feedback and Development

- **Feature requests**: Contact development team
- **Bug reports**: Use standard CoastalME reporting procedures
- **Algorithm contributions**: Follow development guidelines

## Support

### Documentation

- **Technical documentation**: See `CLIFF_COLLAPSE_SYSTEM.md`
- **API documentation**: See source code comments
- **Examples**: Check test suite configurations

### Getting Help

1. **Check this user guide** for common issues
2. **Review technical documentation** for detailed information
3. **Contact development team** for specific problems
4. **Use CoastalME forums** for community support

### Reporting Issues

When reporting issues, please include:
- Configuration file (`.dat`)
- Error messages
- System information (OS, compiler, OpenMP support)
- Expected vs. actual behavior

---

*This user guide covers the basic usage of the CoastalME cliff collapse system. For technical details and development information, see the technical documentation.*