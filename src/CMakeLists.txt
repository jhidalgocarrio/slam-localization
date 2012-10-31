rock_find_pkgconfig(Eigen3 eigen3 REQUIRED)

rock_library(asguard_localization
    SOURCES sckf.cpp lesma.cpp
    HEADERS sckf.hpp lesma.hpp)

rock_executable(asguard_localization_bin main.cpp
    DEPS asguard_localization)

