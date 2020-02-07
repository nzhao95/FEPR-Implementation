cmake_minimum_required(VERSION 3.0)

project(FEPR)

add_subdirectory(External)

add_executable (
	FEPR
	Sources/Main.cpp
	Sources/Error.h
	Sources/Error.cpp
	Sources/Transform.h
	Sources/Camera.h
	Sources/Mesh.h
	Sources/Mesh.cpp
	Sources/MeshLoader.h
	Sources/MeshLoader.cpp
	Sources/ShaderProgram.h
	Sources/ShaderProgram.cpp
	Sources/stb_image.h
)

# Copy the shader files in the binary location.

add_custom_command(TARGET FEPR
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:FEPR> ${CMAKE_CURRENT_SOURCE_DIR})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(FEPR LINK_PRIVATE glad)

target_link_libraries(FEPR LINK_PRIVATE glfw)

target_link_libraries(FEPR LINK_PRIVATE glm)
