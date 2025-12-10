cmake_minimum_required(VERSION 3.10)
project(simpleMuon)

# Geant4 with UI + VIS (and analysis)
find_package(Geant4 REQUIRED ui_all vis_all)

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

set(SOURCES
        main.cc
        src/DetectorConstruction.cc
        src/PrimaryGeneratorAction.cc
        src/ActionInitialization.cc
        src/RunAction.cc
        src/PlaneSD.cc
)

add_executable(simpleMuon ${SOURCES})
target_link_libraries(simpleMuon ${Geant4_LIBRARIES})

# Absolute path to macros dir
target_compile_definitions(simpleMuon PRIVATE
        MACROS_DIR="${PROJECT_SOURCE_DIR}/macros"
)
