rock_testsuite(UKFoMUnitTest UKFoMUnitTest.cpp
    DEPS rover_localization)

rock_testsuite(UsckfUnitTest UsckfUnitTest.cpp
    DEPS rover_localization)

rock_testsuite(DataModelUnitTest DataModelUnitTest.cpp
    DEPS rover_localization)

