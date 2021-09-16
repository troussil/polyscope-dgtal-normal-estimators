project(polyscope-dgtal-normal-estimators)
CMAKE_MINIMUM_REQUIRED(VERSION 3.11)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(deps)

##DGtal Specific settings
include("${PROJECT_BINARY_DIR}/deps/DGtal/DGtalConfig.cmake")
include_directories("${DGTAL_INCLUDE_DIRS}")

SET(SRCs
  demo
  )
  
FOREACH(FILE ${SRCs})
  ADD_EXECUTABLE(${FILE} ${FILE})
  TARGET_LINK_LIBRARIES(${FILE} polyscope DGtal ${DGTAL_LIBRARIES})
ENDFOREACH(FILE)
