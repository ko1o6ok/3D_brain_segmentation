cmake_minimum_required(VERSION 3.10.0)

project(vtk_workspace VERSION 0.1.0 LANGUAGES C CXX)

# Найти ITK
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

# Найти VTK с необходимыми компонентами
find_package(VTK REQUIRED COMPONENTS
    CommonCore
    CommonDataModel
    RenderingCore
    RenderingOpenGL2
    InteractionStyle
    RenderingVolume
    RenderingVolumeOpenGL2
    IOImage
    FiltersCore
    RenderingContext2D
    FiltersSources  # Добавлено для простого теста
)

include(${VTK_USE_FILE})

# Настройка компилятора
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Основная программа с volume rendering
add_executable(vtk_workspace main.cpp)
target_link_libraries(vtk_workspace ${VTK_LIBRARIES} ${ITK_LIBRARIES})

vtk_module_autoinit(TARGETS vtk_workspace  MODULES ${VTK_LIBRARIES} )

# Простая тестовая программа
# add_executable(vtk_simple simple_test.cpp)
# target_link_libraries(vtk_simple ${VTK_LIBRARIES})