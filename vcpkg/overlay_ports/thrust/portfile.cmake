vcpkg_from_github(
  OUT_SOURCE_PATH SOURCE_PATH
  REPO thrust/thrust
  REF 1.15.0
  SHA512 0ba2462b8d113c1843bd90adb1fc750fe82f67f9568903479ec5d505d2f82d4d2446e43a7e0f2b81c29871c75d9eed2f87654be21f47dd9f0efb363b24699ef6
  HEAD_REF main
)
file                     (COPY ${CMAKE_CURRENT_LIST_DIR}/CMakeLists.txt DESTINATION ${SOURCE_PATH})
vcpkg_configure_cmake    (SOURCE_PATH ${SOURCE_PATH} PREFER_NINJA)
vcpkg_install_cmake      ()
vcpkg_fixup_cmake_targets()
file                     (REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/include/thrust/cmake)
file                     (REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/debug)
file                     (INSTALL ${SOURCE_PATH}/LICENSE DESTINATION ${CURRENT_PACKAGES_DIR}/share/${PORT} RENAME copyright)
