rock_library(asguard_localization
    SOURCES Dummy.cpp
    HEADERS Dummy.hpp)

rock_executable(asguard_localization_bin Main.cpp
    DEPS asguard_localization)

