cmake_minimum_required(VERSION 3.5)
project(A1_Robot VERSION 1.1.1)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3 -fopenmp -std=c++17 -ggdb")
set(OGRE_DIR "/home/basit/workspace/raisimLib/build/lib/OGRE/cmake/")

set(LL_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# ===== Path to RaiSim and RaiSimOgre ===== #
if($ENV{USER} STREQUAL "basit")
    find_package(raisim REQUIRED HINTS "/home/basit/workspace/raisimLib/build/lib/cmake/raisim")
    # find_package(raisim REQUIRED HINTS "/home/basit/workspace/raisimLib/raisim/linux/lib/cmake/raisim")
    find_package(raisimOgre REQUIRED HINTS "/home/basit/workspace/raisimOgre/build/cmake")
else()
    # find_package(raisim REQUIRED PATHS "/home/basit/workspace/raisimLib/cmake/raisim")
    # find_package(raisimOgre CONFIG PATHS "/home/basit/workspace/raisimOgre/cmake")
endif()


find_package(ifopt)
find_package(casadi REQUIRED)
find_package(Eigen3 REQUIRED HINTS ${Eigen3_HINT})
find_package(Boost REQUIRED COMPONENTS timer system thread)
find_library(LCM lcm HINTS /usr/local/lib)
find_package(yaml-cpp REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(IPOPT REQUIRED ipopt)

if (raisimOgre_FOUND)
    include_directories(PUBLIC ${LL_DIR}/raisimOgre_include)
    file(COPY ${LL_DIR}/rsc DESTINATION ${CMAKE_BINARY_DIR})
endif()


# include header directories
include_directories(
    ${Boost_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/include
    SYSTEM ${IPOPT_INCLUDE_DIRS}
    SYSTEM ${CASADI_INCLUDE_DIR}
    SYSTEM ${SNOPT_INCLUDE_DIRS}
    ${CODE_GEN_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIRS}
    /home/basit/workspace/libsnopt7_cpp)

# library directories
link_directories(
    ${LIBRARY_DIRS}
    ${IPOPT_LIBRARY_DIRS}
    /home/basit/workspace/libsnopt7_cpp)

# find all the header files
file(GLOB HEADER_FILES_HPP ${CMAKE_SOURCE_DIR}/include/*.hpp)
file(GLOB HEADER_FILES_H ${CMAKE_SOURCE_DIR}/include/*.h)


set(UT_INC ${CMAKE_CURRENT_SOURCE_DIR}/unitree_legged_sdk/include)
set(UT_LINK ${CMAKE_CURRENT_SOURCE_DIR}/unitree_legged_sdk/lib)
set(EXTRA_LIBS -pthread ${UT_LINK}/libunitree_legged_sdk_amd64.so)

# set casadi include directory
set(CASADI_INCLUDE_DIR /usr/local/include/casadi)

# find casadi library
find_library(CASADI_LIBRARY
    NAMES casadi
    HINTS ${CASADI_INCLUDE_DIR}/../lib $ENV{CASADI_PREFIX}/lib)
if(CASADI_LIBRARY)
    set(CASADI_LIBRARIES ${CASADI_LIBRARIES} ${CASADI_LIBRARY})
endif()



# ===== Dynamics Library ===== #

# Settings:
set(COMPILE_CORIOLIS FALSE) # Option to compile the coriolis terms. Not implemented in code though.

set(dyn_srcs 
    "${LL_DIR}/dynamics/A1_New_wFoot/A1_D_B_G.cpp"
    "${LL_DIR}/dynamics/A1_New_wFoot/A1_FK.cpp"
    "${LL_DIR}/dynamics/A1_New_wFoot/A1_FK.cpp"
    "${LL_DIR}/dynamics/A1_New_wFoot/A1_Jaco.cpp"
    "${LL_DIR}/dynamics/A1_New_wFoot/A1_Jaco_dot.cpp")
if (COMPILE_CORIOLIS)
    set(dyn_srcs ${dyn_srcs} "${LL_DIR}/dynamics/A1_New_wFoot/A1_C.cpp")
endif()
add_library(dyn_lib ${dyn_srcs})
target_include_directories(dyn_lib PUBLIC
    "${LL_DIR}/dynamics"
    "${LL_DIR}/dynamics/A1_New_wFoot"
)

# ===== Controller Library ===== #
add_library(ctrl_lib
    "${LL_DIR}/src/Parameters.cpp"
    "${LL_DIR}/src/RobotModel.cpp"
    "${LL_DIR}/src/LocoWrapper.cpp"
    "${LL_DIR}/src/LowLevelCtrl.cpp"
    "${LL_DIR}/src/VirtualConstraints.cpp"
    "${LL_DIR}/src/ContactEst.cpp"
    "${LL_DIR}/src/MotionPlanner.cpp"
    "${LL_DIR}/src/DataLog.cpp"
    "${LL_DIR}/src/MPC_dist.cpp"
    "${LL_DIR}/src/nmpc.cpp"
)
target_include_directories(ctrl_lib PUBLIC 
    ${LL_DIR}/include
    ${EIGEN3_INCLUDE_DIRS}
    ${LL_DIR}/optimization/iSWIFT/include
    ${LL_DIR}/optimization/iSWIFT/ldl/include
    ${LL_DIR}/optimization/iSWIFT/cpp_wrapper
)
target_link_libraries(ctrl_lib PRIVATE dyn_lib)

# ===== High Level Library (basic test only) ===== #
# add_library(hl_lib
# 	"${LL_DIR}/src/forceUpdate.cpp"
# )

# ===== QP SWIFT Library ===== #
add_library(optim_lib
    "${LL_DIR}/optimization/iSWIFT/src/Auxilary.c"
    "${LL_DIR}/optimization/iSWIFT/src/Prime.c"
    "${LL_DIR}/optimization/iSWIFT/src/timer.c"
    "${LL_DIR}/optimization/iSWIFT/ldl/src/ldl.c"
    "${LL_DIR}/optimization/iSWIFT/cpp_wrapper/iswift_qp.cpp"
)
target_include_directories(optim_lib PUBLIC 
    ${EIGEN3_INCLUDE_DIRS}
    ${LL_DIR}/optimization/iSWIFT/include
    ${LL_DIR}/optimization/iSWIFT/ldl/include
    ${LL_DIR}/optimization/iSWIFT/cpp_wrapper
)


# ===== Executables ===== #
if (raisimOgre_FOUND)
    add_compile_definitions(RUNNING_SIM) # adds the toe offset to be used in simulation (exp needs no offset)
    add_executable(Run_Sim src/A1_Sim.cpp)
    target_link_libraries(Run_Sim PRIVATE ctrl_lib optim_lib dyn_lib raisim::raisimOgre ifopt::ifopt_snopt ipopt -Llibsnopt7_cpp -Llibassimp ${CASADI_LIBRARIES} yaml-cpp)
    target_compile_definitions(Run_Sim PRIVATE "-DEXAMPLE_ROBOT_RESOURCE_DIR=${LL_DIR}/rsc/")
# else()
    add_executable(nmpc_exp src/nmpc_exp.cpp)
    target_link_libraries(nmpc_exp PRIVATE ctrl_lib optim_lib dyn_lib ${EXTRA_LIBS} ${LCM} ${Boost_LIBRARIES} ifopt::ifopt_snopt ipopt -Llibsnopt7_cpp -Llibassimp ${CASADI_LIBRARIES} yaml-cpp)
    # target_link_libraries(nmpc_exp PRIVATE c  trl_lib optim_lib dyn_lib ${EXTRA_LIBS} ${LCM} ${Boost_LIBRARIES})
    target_include_directories(nmpc_exp PRIVATE ${UT_INC})
    
    # add_executable(A1_Experiment_w_HL src/A1_Exp_w_HL.cpp)
    # target_link_libraries(A1_Experiment_w_HL PRIVATE ctrl_lib hl_lib optim_lib dyn_lib ${EXTRA_LIBS} ${LCM} ${Boost_LIBRARIES})
    # target_include_directories(A1_Experiment_w_HL PRIVATE ${UT_INC})
    
    add_executable(restore src/restoreConfig.cpp)
    target_link_libraries(restore PRIVATE ${EXTRA_LIBS} ${LCM} ${Boost_LIBRARIES})
    target_include_directories(restore PRIVATE ${UT_INC})
endif()

# EOF