# SPDX-License-Identifier: LGPL-2.1-only

set( test_name testDualContinuumExchangeUpwinding )

blt_add_executable( NAME ${test_name}
                    SOURCES testDualContinuumExchangeUpwinding.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON gtest gtest_main multiPhysicsSolvers ${parallelDeps} )

geos_add_test( NAME ${test_name}
               COMMAND ${test_name} )
