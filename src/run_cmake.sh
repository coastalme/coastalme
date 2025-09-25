#!/bin/sh

set -e  # Exit on any error

# Default configuration
DEFAULT_BUILD_TYPE="Debug"
DEFAULT_CSHORE_LIBRARY="SHARED"
DEFAULT_CSHORE_INOUT="ARG"
DEFAULT_GENERATOR="Ninja"

# Configuration variables with defaults
BUILD_TYPE="${DEFAULT_BUILD_TYPE}"
CSHORE_LIBRARY="${DEFAULT_CSHORE_LIBRARY}"
CSHORE_INOUT="${DEFAULT_CSHORE_INOUT}"
GENERATOR="${DEFAULT_GENERATOR}"
VERBOSE=false
CLEAN_BUILD=false
BUILD_CSHORE=false
CMAKE_COMPILER_ARGS=""

show_help() {
    cat << EOF
CoastalME Build Script

USAGE:
    $0 [OPTIONS]

OPTIONS:
    -s, --build-cshore        Build CShore library before CMake (default: false)
    -t, --build-type TYPE     Build type: Debug, Release, Prerelease, gcov, Callgrind
                              (default: ${DEFAULT_BUILD_TYPE})
    -l, --library TYPE        CShore library type: STATIC, SHARED
                              (default: ${DEFAULT_CSHORE_LIBRARY})
    -i, --inout METHOD        CShore I/O method: FILE, ARG, BOTH
                              (default: ${DEFAULT_CSHORE_INOUT})
    -g, --generator GENERATOR CMake generator: Ninja, "Unix Makefiles"
                              (default: ${DEFAULT_GENERATOR})
    -c, --clean               Clean build (remove CMakeCache.txt and build files)
    -v, --verbose             Enable verbose output
    -h, --help                Show this help message

EXAMPLES:
    $0                        # Use all defaults (no CShore build)
    $0 -s -t Release -l STATIC # Build CShore with Release and static library
    $0 --build-cshore --clean --verbose # Build CShore with clean and verbose output

NOTES:
    - CShore library build is optional (use -s/--build-cshore to enable)
    - macOS automatically configures gcc-15 for OpenMP support
    - Prerelease builds should not be run with debugger (gdb)
EOF
}

log_info() {
    echo "INFO: $1"
}

log_error() {
    echo "ERROR: $1" >&2
}

validate_build_type() {
    case "$1" in
        Debug|Release|Prerelease|RelWithDebInfo|MinSizeRel|gcov|Callgrind)
            return 0
            ;;
        *)
            log_error "Invalid build type: $1"
            log_error "Valid types: Debug, Release, Prerelease, RelWithDebInfo, MinSizeRel, gcov, Callgrind"
            return 1
            ;;
    esac
}

validate_library_type() {
    case "$1" in
        STATIC|SHARED)
            return 0
            ;;
        *)
            log_error "Invalid library type: $1"
            log_error "Valid types: STATIC, SHARED"
            return 1
            ;;
    esac
}

validate_inout_method() {
    case "$1" in
        FILE|ARG|BOTH)
            return 0
            ;;
        *)
            log_error "Invalid I/O method: $1"
            log_error "Valid methods: FILE, ARG, BOTH"
            return 1
            ;;
    esac
}

parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -t|--build-type)
                BUILD_TYPE="$2"
                validate_build_type "$BUILD_TYPE" || exit 1
                shift 2
                ;;
            -l|--library)
                CSHORE_LIBRARY="$2"
                validate_library_type "$CSHORE_LIBRARY" || exit 1
                shift 2
                ;;
            -i|--inout)
                CSHORE_INOUT="$2"
                validate_inout_method "$CSHORE_INOUT" || exit 1
                shift 2
                ;;
            -g|--generator)
                GENERATOR="$2"
                shift 2
                ;;
            -s|--build-cshore)
                BUILD_CSHORE=true
                shift
                ;;
            -c|--clean)
                CLEAN_BUILD=true
                shift
                ;;
            -v|--verbose)
                VERBOSE=true
                shift
                ;;
            -h|--help)
                show_help
                exit 0
                ;;
            *)
                log_error "Unknown option: $1"
                show_help
                exit 1
                ;;
        esac
    done
}

setup_compiler_args() {
    # macOS compiler setup for OpenMP support
    if [[ $OSTYPE == 'darwin'* ]]; then
        log_info "macOS detected - configuring gcc-15 for OpenMP compatibility"
        export CC=gcc-15
        export CXX=g++-15
        CMAKE_COMPILER_ARGS="-DCMAKE_C_COMPILER=gcc-15 -DCMAKE_CXX_COMPILER=g++-15"
        [ "$VERBOSE" = true ] && log_info "Compiler set to gcc-15/g++-15"
    else
        log_info "Using default compiler, OpenMP may not be available"
        log_info "For OpenMP support, install GCC or install libomp for Clang"
        CMAKE_COMPILER_ARGS=""
    fi
}

clean_build_artifacts() {
    if [ "$CLEAN_BUILD" = true ]; then
        log_info "Cleaning build artifacts..."
        rm -f CMakeCache.txt
        rm -rf CMakeFiles/
        rm -f cmake_install.cmake
        rm -f Makefile
        rm -f ./lib/*
        [ "$VERBOSE" = true ] && log_info "Build artifacts cleaned"
    else
        rm -f CMakeCache.txt
    fi
}

build_cshore_library() {
    log_info "Building CShore library..."
    rm -f ./lib/*
    cd cshore

    if [[ $OSTYPE == 'darwin'* ]]; then
        ./make_cshore_lib.sh "$CSHORE_LIBRARY" "$BUILD_TYPE" "$CSHORE_INOUT"
    else
        ./make_cshore_lib.sh
    fi

    cd ..
    [ "$VERBOSE" = true ] && log_info "CShore library build completed"
}

run_cmake() {
    local cmake_args="-G \"$GENERATOR\" -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCSHORE_LIBRARY=$CSHORE_LIBRARY -DCSHORE_INOUT=$CSHORE_INOUT"

    if [ -n "$CMAKE_COMPILER_ARGS" ]; then
        cmake_args="$cmake_args $CMAKE_COMPILER_ARGS"
    fi

    log_info "Running CMake with configuration:"
    log_info "  Build Type: $BUILD_TYPE"
    log_info "  CShore Library: $CSHORE_LIBRARY"
    log_info "  CShore I/O Method: $CSHORE_INOUT"
    log_info "  Generator: $GENERATOR"
    log_info "  CShore Built: $BUILD_CSHORE"

    if [ "$VERBOSE" = true ]; then
        log_info "CMake command: cmake $cmake_args ."
    fi

    eval "cmake $cmake_args ."
}

run_build() {
    log_info "Building project..."

    case "$GENERATOR" in
        "Ninja")
            ninja
            ;;
        "Unix Makefiles")
            make -j$(nproc 2>/dev/null || echo 4)
            ;;
        *)
            log_error "Unsupported generator for build: $GENERATOR"
            exit 1
            ;;
    esac
}

print_summary() {
    echo ""
    echo "================================================================="
    echo "BUILD COMPLETED SUCCESSFULLY"
    echo "================================================================="
    echo "Configuration:"
    echo "  Build Type: $BUILD_TYPE"
    echo "  CShore Library: $CSHORE_LIBRARY"
    echo "  CShore I/O Method: $CSHORE_INOUT"
    echo "  Generator: $GENERATOR"
    echo "  CShore Built: $BUILD_CSHORE"
    echo "================================================================="
    echo ""
}

main() {
    parse_arguments "$@"

    echo ""
    echo "================================================================="
    echo "CoastalME Build Script"
    echo "================================================================="

    setup_compiler_args
    clean_build_artifacts

    if [ "$BUILD_CSHORE" = true ]; then
        build_cshore_library
        echo ""
        echo "================================================================="
        echo ""
    fi

    run_cmake
    run_build
    print_summary
}

main "$@"

