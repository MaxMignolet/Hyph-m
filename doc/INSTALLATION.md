The build system uses meson. Meson relies on pkg-config (and optionnaly cmake) to find the required dependencies. If Meson is unable to find a dependency, please make sure that it lies within pkg-config search path or cmake search path:
 - For pkg-config:
        `export PKG_CONFIG_PATH="/path/to/dependency/lib/pkgconfig:$PKG_CONFIG_PATH"`
 - For CMake:
        `export CMAKE_LIBRARY_PATH="/path/to/dependency/lib:$CMAKE_LIBRARY_PATH"`
pkg-config requires a `.pc` file. To check that pkg-config is able to find a dependency, run: `pkg-config --path dependency_name` (--path or --cflags or --libs)

Alternatively, a makefile is provided. However the dependencies' paths are hardcoded and needs to be updated to fit your system configuration.
