rock_find_pkgconfig(Eigen3 eigen3 REQUIRED)

set (LOCALIZATION_SRCS
    Sckf.cpp
    Measurement.cpp
    Lesma.cpp
    DataModel.cpp
    DeadReckon.cpp
    )

set (LOCALIZATION_HDRS
    Configuration.hpp
    Sckf.hpp
    Measurement.hpp
    Lesma.hpp
    DataModel.hpp
    DeadReckon.hpp
    DataTypes.hpp
    )

rock_library(rover_localization
    SOURCES ${LOCALIZATION_SRCS}
    DEPS_PKGCONFIG base-types base-lib yaml-cpp envire
    HEADERS ${LOCALIZATION_HDRS})

rock_executable(rover_localization_bin main.cpp
    DEPS rover_localization)

