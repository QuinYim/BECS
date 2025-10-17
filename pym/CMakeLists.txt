#cmake_minimum_required(VERSION 3.29)
#project(operation)
#
#set(CMAKE_CXX_STANDARD 20)
#
#
#
## 查找 OpenSSL
#find_package(OpenSSL REQUIRED)
#
## 查找 GMP
#find_package(GMP REQUIRED)
#
## 查找 PBC
#find_path(PBC_INCLUDE_DIR pbc/pbc.h)
#find_library(PBC_LIBRARIES NAMES pbc)
#
#if (NOT PBC_INCLUDE_DIR OR NOT PBC_LIBRARIES)
#    message(FATAL_ERROR "Could not find PBC library")
#endif()
#
#include_directories(${PBC_INCLUDE_DIR} ${GMP_INCLUDE_DIR} ${OPENSSL_INCLUDE_DIR})
#
#add_executable(operation main.cpp)
#
## 链接库
#target_link_libraries(operation PRIVATE ${PBC_LIBRARIES} GMP::GMP OpenSSL::SSL OpenSSL::Crypto)
#
#
#

cmake_minimum_required(VERSION 3.29)
project(operation)

set(CMAKE_CXX_STANDARD 20)

# 设置 CMake 前缀路径以包含 /usr/local
# 这将帮助 CMake 在 /usr/local 目录下查找库和包
set(CMAKE_PREFIX_PATH "/usr/local" ${CMAKE_PREFIX_PATH})

# 设置 OpenSSL 的根目录
set(OPENSSL_ROOT_DIR "/usr/local/openssl")

# 查找 OpenSSL
find_package(OpenSSL REQUIRED)
if (NOT OpenSSL_FOUND)
    message(FATAL_ERROR "OpenSSL not found. Please set OPENSSL_ROOT_DIR.")
endif()

# 手动查找 GMP 的头文件和库文件
find_path(GMP_INCLUDE_DIR gmp.h
        PATHS /usr/local/include /usr/include
)

find_library(GMP_LIBRARIES gmp
        PATHS /usr/local/lib /usr/lib/x86_64-linux-gnu
)

if (NOT GMP_INCLUDE_DIR OR NOT GMP_LIBRARIES)
    message(FATAL_ERROR "Could not find GMP library. Please ensure GMP is installed correctly.")
endif()

# 查找 PBC 的头文件和库文件
find_path(PBC_INCLUDE_DIR pbc/pbc.h
        PATHS /usr/include /usr/local/include
)

find_library(PBC_LIBRARIES pbc
        PATHS /usr/lib /usr/local/lib
)

if (NOT PBC_INCLUDE_DIR OR NOT PBC_LIBRARIES)
    message(FATAL_ERROR "Could not find PBC library. Please install PBC development package.")
endif()

# 添加可执行文件
add_executable(operation main.cpp)

# 包含目录
target_include_directories(operation PRIVATE
        ${OPENSSL_INCLUDE_DIR}
        ${GMP_INCLUDE_DIR}
        ${PBC_INCLUDE_DIR}
)

# 链接库
target_link_libraries(operation PRIVATE
        OpenSSL::SSL
        OpenSSL::Crypto
        ${GMP_LIBRARIES}
        ${PBC_LIBRARIES}
)

# 打印找到的库信息，帮助调试
message(STATUS "OpenSSL include dir: ${OPENSSL_INCLUDE_DIR}")
message(STATUS "OpenSSL libraries: ${OPENSSL_LIBRARIES}")
message(STATUS "GMP include dir: ${GMP_INCLUDE_DIR}")
message(STATUS "GMP libraries: ${GMP_LIBRARIES}")
message(STATUS "PBC include dir: ${PBC_INCLUDE_DIR}")
message(STATUS "PBC libraries: ${PBC_LIBRARIES}")
