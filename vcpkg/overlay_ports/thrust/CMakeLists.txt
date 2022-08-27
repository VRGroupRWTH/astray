cmake_minimum_required    (VERSION 3.12 FATAL_ERROR)
project                   (thrust LANGUAGES CXX)
include                   (GNUInstallDirs)
add_library               (thrust INTERFACE)
target_include_directories(thrust INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
install(
    TARGETS thrust EXPORT thrust-config
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/
)
install(
    EXPORT thrust-config
    NAMESPACE thrust::
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/thrust
)
install(
    DIRECTORY   ${CMAKE_SOURCE_DIR}/thrust
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp*" PATTERN "*.inl*" PATTERN "*.h*"
)