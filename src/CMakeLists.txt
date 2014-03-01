set (LOCALIZATION_SRCS
    core/Transform.cpp
    )

set (LOCALIZATION_HDRS
    Configuration.hpp
    core/DataModel.hpp
    core/DeadReckon.hpp
    core/Types.hpp
    core/Transform.hpp
    tools/Util.hpp
    tools/Analysis.hpp
    filters/Usckf.hpp
    filters/UsckfError.hpp
    filters/MtkWrap.hpp
    filters/State.hpp
    filters/ProcessModels.hpp
    filters/MeasurementModels.hpp
    filters/IIR.hpp
    )


rock_library(localization
    SOURCES ${LOCALIZATION_SRCS}
    DEPS_PKGCONFIG eigen3 base-types base-lib #yaml-cpp
    DEPS_CMAKE LAPACK
    HEADERS ${LOCALIZATION_HDRS})

