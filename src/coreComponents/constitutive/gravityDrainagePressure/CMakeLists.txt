set( gravityDrainagePressureSources
     GravityDrainagePressureBase.cpp
     SimpleGravityDrainagePressure.cpp
     ThomasGasOilGravityDrainagePressure.cpp
   )

set( gravityDrainagePressureHeaders
     GravityDrainagePressureBase.hpp
     GravityDrainagePressureFields.hpp
     SimpleGravityDrainagePressure.hpp
     ThomasGasOilGravityDrainagePressure.hpp
   )

set( dependencyList
     constitutive_common
     constitutive_fluid
   )

blt_add_library(
  NAME                  gravityDrainagePressure
  SOURCES               ${gravityDrainagePressureSources}
  HEADERS               ${gravityDrainagePressureHeaders}
  DEPENDS_ON            ${dependencyList}
  OBJECT                TRUE
)

target_include_directories( gravityDrainagePressure
                            PUBLIC
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
                          )
