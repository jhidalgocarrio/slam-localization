rock_find_pkgconfig(Eigen3 eigen3 REQUIRED)

rock_library(rover_localization
    SOURCES sckf.cpp measurement.cpp lesma.cpp
    HEADERS configuration.hpp sckf.hpp measurement.hpp lesma.hpp)

rock_executable(rover_localization_bin main.cpp
    DEPS rover_localization)

