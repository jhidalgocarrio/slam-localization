set (LOCALIZATION_SRCS
    core/Transform.cpp
    )

set (LOCALIZATION_HDRS
    core/DataModel.hpp
    core/DeadReckon.hpp
    core/Types.hpp
    core/Transform.hpp
    tools/Util.hpp
    tools/Analysis.hpp
    filters/IIR.hpp
    filters/FIR.hpp
    )


rock_library(localization
    SOURCES ${LOCALIZATION_SRCS}
    DEPS_PKGCONFIG eigen3 base-types base-lib #yaml-cpp
    DEPS_CMAKE LAPACK
    HEADERS ${LOCALIZATION_HDRS})

