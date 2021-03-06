itk_module_test()

#Set output path
set(TESTING_OUTPUT_PATH "${CMAKE_CURRENT_BINARY_DIR}/Output")
make_directory(${TESTING_OUTPUT_PATH})

set(DATA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Input)
set(BASELINE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/baseline)
  
set( _test
Test_SpeedToPath_GradientDescent_2D.cxx
Test_SpeedToPath_GradientDescent_3D.cxx
Test_SpeedToPath_IterateNeighborhood_2D.cxx
Test_SpeedToPath_IterateNeighborhood_ExtendedSeed_2D.cxx
Test_SpeedToPath_IterateNeighborhood_3D.cxx
Test_SpeedToPath_RegularStepGradientDescent_2D.cxx
Test_SpeedToPath_RegularStepGradientDescent_3D.cxx
)

CreateTestDriver( MinimalPathExtraction "${MinimalPathExtraction-Test_LIBRARIES}" "${_test}" )

#Add tests for the source images
#Noise-01
add_test(
  NAME MinimalPathExtraction_Noise_01_NeighborhoodIterate
  COMMAND MinimalPathExtractionTestDriver
    --compare ${TESTING_OUTPUT_PATH}/Noise-01-IterateNeighborhood.png ${BASELINE_DIR}/Noise-01-IterateNeighborhood.png
    Test_SpeedToPath_IterateNeighborhood_2D
    ${TESTING_OUTPUT_PATH}/Noise-01-IterateNeighborhood.png
    ${DATA_DIR}/Noise-01-Speed.mhd
    ${DATA_DIR}/Noise-01.path
    1.0
    1.0
)
add_test(
  NAME MinimalPathExtraction_Noise_01_GradientDescent
  COMMAND MinimalPathExtractionTestDriver
    --compare ${TESTING_OUTPUT_PATH}/Noise-01-GradientDescent.png ${BASELINE_DIR}/Noise-01-GradientDescent.png
    Test_SpeedToPath_GradientDescent_2D
    ${TESTING_OUTPUT_PATH}/Noise-01-GradientDescent.png
    ${DATA_DIR}/Noise-01-Speed.mhd
    ${DATA_DIR}/Noise-01.path
    1.0
    1000
)
add_test(
  NAME MinimalPathExtraction_Noise_01_RegularStepGradientDescent
  COMMAND MinimalPathExtractionTestDriver
    --compare ${TESTING_OUTPUT_PATH}/Noise-01-RegularStepGradientDescent.png ${BASELINE_DIR}/Noise-01-RegularStepGradientDescent.png
    Test_SpeedToPath_RegularStepGradientDescent_2D
    ${TESTING_OUTPUT_PATH}/Noise-01-RegularStepGradientDescent.png
    ${DATA_DIR}/Noise-01-Speed.mhd
    ${DATA_DIR}/Noise-01.path
    1.0
    1000
    1.0
    0.5
)
#Synthetic-01
add_test(
  NAME MinimalPathExtraction_Synthetic_01_NeighborhoodIterate
  COMMAND MinimalPathExtractionTestDriver
    --compare ${TESTING_OUTPUT_PATH}/Synthetic-01-IterateNeighborhood.png ${BASELINE_DIR}/Synthetic-01-IterateNeighborhood.png
    Test_SpeedToPath_IterateNeighborhood_2D
    ${TESTING_OUTPUT_PATH}/Synthetic-01-IterateNeighborhood.png
    ${DATA_DIR}/Synthetic-01-Speed.mhd
    ${DATA_DIR}/Synthetic-01.path
    1.0
    1.0
)
add_test(
  NAME MinimalPathExtraction_Synthetic_01_NeighborhoodIterate_ExtendedSeed
  COMMAND MinimalPathExtractionTestDriver
    --compare ${TESTING_OUTPUT_PATH}/Synthetic-01-IterateNeighborhood_extended.png ${BASELINE_DIR}/Synthetic-01-IterateNeighborhood_extended.png
    Test_SpeedToPath_IterateNeighborhood_ExtendedSeed_2D
    ${TESTING_OUTPUT_PATH}/Synthetic-01-IterateNeighborhood_extended.png
    ${DATA_DIR}/Synthetic-01-Speed.mhd
    ${DATA_DIR}/Synthetic-01_seed.png
    1.0
    1.0
)
add_test(
  NAME MinimalPathExtraction_Synthetic_01_GradientDescent
  COMMAND MinimalPathExtractionTestDriver
    --compare ${TESTING_OUTPUT_PATH}/Synthetic-01-GradientDescent.png ${BASELINE_DIR}/Synthetic-01-GradientDescent.png
    Test_SpeedToPath_GradientDescent_2D
    ${TESTING_OUTPUT_PATH}/Synthetic-01-GradientDescent.png
    ${DATA_DIR}/Synthetic-01-Speed.mhd
    ${DATA_DIR}/Synthetic-01.path
    1.0
    1000
)
add_test(
  NAME MinimalPathExtraction_Synthetic_01_RegularStepGradientDescent
  COMMAND MinimalPathExtractionTestDriver
    --compare ${TESTING_OUTPUT_PATH}/Synthetic-01-RegularStepGradientDescent.png ${BASELINE_DIR}/Synthetic-01-RegularStepGradientDescent.png
    Test_SpeedToPath_RegularStepGradientDescent_2D
    ${TESTING_OUTPUT_PATH}/Synthetic-01-RegularStepGradientDescent.png
    ${DATA_DIR}/Synthetic-01-Speed.mhd
    ${DATA_DIR}/Synthetic-01.path
    1.0
    1000
    1.0
    0.5
)
#Synthetic-02
add_test(
  NAME MinimalPathExtraction_Synthetic_02_NeighborhoodIterate
  COMMAND MinimalPathExtractionTestDriver
    --compare ${TESTING_OUTPUT_PATH}/Synthetic-02-IterateNeighborhood.png ${BASELINE_DIR}/Synthetic-02-IterateNeighborhood.png
    Test_SpeedToPath_IterateNeighborhood_2D
    ${TESTING_OUTPUT_PATH}/Synthetic-02-IterateNeighborhood.png
    ${DATA_DIR}/Synthetic-02-Speed.mhd
    ${DATA_DIR}/Synthetic-02.path
    1.0
    1.0
)
add_test(
  NAME MinimalPathExtraction_Synthetic_02_GradientDescent
  COMMAND MinimalPathExtractionTestDriver
  --compare ${TESTING_OUTPUT_PATH}/Synthetic-02-GradientDescent.png ${BASELINE_DIR}/Synthetic-02-GradientDescent.png
  Test_SpeedToPath_GradientDescent_2D
    ${TESTING_OUTPUT_PATH}/Synthetic-02-GradientDescent.png
    ${DATA_DIR}/Synthetic-02-Speed.mhd
    ${DATA_DIR}/Synthetic-02.path
    1.0
    3000
)
add_test(
  NAME MinimalPathExtraction_Synthetic_02_RegularStepGradientDescent
  COMMAND MinimalPathExtractionTestDriver
    --compare ${TESTING_OUTPUT_PATH}/Synthetic-02-RegularStepGradientDescent.png ${BASELINE_DIR}/Synthetic-02-RegularStepGradientDescent.png
    Test_SpeedToPath_RegularStepGradientDescent_2D
    ${TESTING_OUTPUT_PATH}/Synthetic-02-RegularStepGradientDescent.png
    ${DATA_DIR}/Synthetic-02-Speed.mhd
    ${DATA_DIR}/Synthetic-02.path
    1.0
    3000
    1.0
    0.5
)
#Synthetic-03
add_test(
  NAME MinimalPathExtraction_Synthetic_03_NeighborhoodIterate
  COMMAND MinimalPathExtractionTestDriver
    --compare ${TESTING_OUTPUT_PATH}/Synthetic-03-IterateNeighborhood.png ${BASELINE_DIR}/Synthetic-03-IterateNeighborhood.png
    Test_SpeedToPath_IterateNeighborhood_2D
    ${TESTING_OUTPUT_PATH}/Synthetic-03-IterateNeighborhood.png
    ${DATA_DIR}/Synthetic-03-Speed.mhd
    ${DATA_DIR}/Synthetic-03.path
    1.0
    1.0
)
#Synthetic-04
add_test(
  NAME MinimalPathExtraction_Synthetic_04_NeighborhoodIterate
  COMMAND MinimalPathExtractionTestDriver
    --compare ${TESTING_OUTPUT_PATH}/Synthetic-04-IterateNeighborhood.mhd ${BASELINE_DIR}/Synthetic-04-IterateNeighborhood.mhd
    Test_SpeedToPath_IterateNeighborhood_3D
    ${TESTING_OUTPUT_PATH}/Synthetic-04-IterateNeighborhood.mhd
    ${DATA_DIR}/Synthetic-04-Speed.mhd
    ${DATA_DIR}/Synthetic-04.path
    1.0
    1.0
)
add_test(
  NAME MinimalPathExtraction_Synthetic_04_GradientDescent
  COMMAND MinimalPathExtractionTestDriver
    --compare ${TESTING_OUTPUT_PATH}/Synthetic-04-GradientDescent.mhd ${BASELINE_DIR}/Synthetic-04-GradientDescent.mhd
    Test_SpeedToPath_GradientDescent_3D
    ${TESTING_OUTPUT_PATH}/Synthetic-04-GradientDescent.mhd
    ${DATA_DIR}/Synthetic-04-Speed.mhd
    ${DATA_DIR}/Synthetic-04.path
    1.0
    3000
)
add_test(
  NAME MinimalPathExtraction_Synthetic_04_RegularStepGradientDescent
  COMMAND MinimalPathExtractionTestDriver
    --compare ${TESTING_OUTPUT_PATH}/Synthetic-04-RegularStepGradientDescent.mhd ${BASELINE_DIR}/Synthetic-04-RegularStepGradientDescent.mhd
    Test_SpeedToPath_RegularStepGradientDescent_3D
    ${TESTING_OUTPUT_PATH}/Synthetic-04-RegularStepGradientDescent.mhd
    ${DATA_DIR}/Synthetic-04-Speed.mhd
    ${DATA_DIR}/Synthetic-04.path
    1.0
    3000
    1.0
    0.5
)
#Real-DSA-01-01
add_test(
  NAME MinimalPathExtraction_Real_DSA_01_01_NeighborhoodIterate
  COMMAND MinimalPathExtractionTestDriver
    --compare ${TESTING_OUTPUT_PATH}/Real-DSA-01-01-IterateNeighborhood.png ${BASELINE_DIR}/Real-DSA-01-01-IterateNeighborhood.png
    Test_SpeedToPath_IterateNeighborhood_2D
    ${TESTING_OUTPUT_PATH}/Real-DSA-01-01-IterateNeighborhood.png
    ${DATA_DIR}/Real-DSA-01-Speed-01.mhd
    ${DATA_DIR}/Real-DSA-01.path
    1.0
    0.05
)
add_test(
  NAME MinimalPathExtraction_Real_DSA_01_01_RegularStepGradientDescent
  COMMAND MinimalPathExtractionTestDriver
    --compare ${TESTING_OUTPUT_PATH}/Real-DSA-01-01-RegularStepGradientDescent.png ${BASELINE_DIR}/Real-DSA-01-01-RegularStepGradientDescent.png
    Test_SpeedToPath_RegularStepGradientDescent_2D
    ${TESTING_OUTPUT_PATH}/Real-DSA-01-01-RegularStepGradientDescent.png
    ${DATA_DIR}/Real-DSA-01-Speed-01.mhd
    ${DATA_DIR}/Real-DSA-01.path
    3.0
    4000
    0.25
    0.5
)
#Real-DSA-01-02
add_test(
  NAME MinimalPathExtraction_Real_DSA_01_02_NeighborhoodIterate
  COMMAND MinimalPathExtractionTestDriver
    --compare ${TESTING_OUTPUT_PATH}/Real-DSA-01-02-IterateNeighborhood.png ${BASELINE_DIR}/Real-DSA-01-02-IterateNeighborhood.png
    Test_SpeedToPath_IterateNeighborhood_2D
    ${TESTING_OUTPUT_PATH}/Real-DSA-01-02-IterateNeighborhood.png
    ${DATA_DIR}/Real-DSA-01-Speed-02.mhd
    ${DATA_DIR}/Real-DSA-01.path
    1.0
    0.05
)
add_test(
  NAME MinimalPathExtraction_Real_DSA_01_02_RegularStepGradientDescent
  COMMAND MinimalPathExtractionTestDriver
    --compare ${TESTING_OUTPUT_PATH}/Real-DSA-01-02-RegularStepGradientDescent.png ${BASELINE_DIR}/Real-DSA-01-02-RegularStepGradientDescent.png
    Test_SpeedToPath_RegularStepGradientDescent_2D
    ${TESTING_OUTPUT_PATH}/Real-DSA-01-02-RegularStepGradientDescent.png
    ${DATA_DIR}/Real-DSA-01-Speed-02.mhd
    ${DATA_DIR}/Real-DSA-01.path
    3.0
    4000
    0.25
    0.5
)
