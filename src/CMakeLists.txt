set (LOCALIZATION_SRCS
    Lesma.cpp
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
    filters/IIR.hpp
    mtk/SOn.hpp
    )


rock_library(rover_localization
    SOURCES ${LOCALIZATION_SRCS}
    DEPS_PKGCONFIG eigen3 base-types base-lib yaml-cpp odometry envire
    DEPS_CMAKE LAPACK
    HEADERS ${LOCALIZATION_HDRS})

