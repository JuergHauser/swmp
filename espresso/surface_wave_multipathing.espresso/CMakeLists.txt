include(FetchContent)
project(espresso_swmp)
cmake_minimum_required(VERSION 3.11.0)

FetchContent_Declare(
  swmp
  GIT_REPOSITORY https://github.com/JuergHauser/swmp.git
  GIT_TAG main
)

FetchContent_MakeAvailable(swmp)

