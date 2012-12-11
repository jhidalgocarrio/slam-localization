rock_find_pkgconfig(Eigen3 eigen3 REQUIRED)

rock_library(rover_localization
    SOURCES sckf.cpp lesma.cpp
    HEADERS sckf.hpp lesma.hpp)

rock_executable(rover_localization_bin main.cpp
    DEPS rover_localization)

