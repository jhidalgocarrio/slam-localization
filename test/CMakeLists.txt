rock_testsuite(UKFoMUnitTest UKFoMUnitTest.cpp
    DEPS rover_localization
    DEPS_CMAKE LAPACK)

rock_testsuite(UsckfUnitTest UsckfUnitTest.cpp
    DEPS rover_localization
    DEPS_CMAKE LAPACK)

rock_testsuite(DataModelUnitTest DataModelUnitTest.cpp
    DEPS rover_localization
    DEPS_CMAKE LAPACK)

