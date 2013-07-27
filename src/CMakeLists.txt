rock_find_pkgconfig(Eigen3 eigen3 REQUIRED)

set (LOCALIZATION_SRCS
    Util.cpp
    Lesma.cpp
    DataModel.cpp
    DeadReckon.cpp
    )

set (LOCALIZATION_HDRS
    Configuration.hpp
    Util.hpp
    Lesma.hpp
    DataModel.hpp
    DeadReckon.hpp
    DataTypes.hpp
    filters/Usckf.hpp
    filters/MtkWrap.hpp
    filters/State.hpp
    filters/ProcessModels.hpp
    filters/MeasurementModels.hpp
    )

rock_library(rover_localization
    SOURCES ${LOCALIZATION_SRCS}
    DEPS_PKGCONFIG base-types base-lib yaml-cpp odometry envire
    HEADERS ${LOCALIZATION_HDRS})

