# ---- Executable ----

project(PTree)

file(GLOB_RECURSE SOURCES "src/*.cpp")
file(GLOB_RECURSE HEADERS "include/*.hpp" "*.h" "*.inl")

add_library(ptree ${SOURCES} ${HEADERS})
target_include_directories(ptree PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(ptree PUBLIC glm)

