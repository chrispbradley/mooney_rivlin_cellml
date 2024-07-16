add_executable(${ExampleName} ${ExampleName}.F90)

target_link_libraries(${ExampleName} PRIVATE OpenCMISS::libOpenCMISS_Fortran)

install(TARGETS ${ExampleName} DESTINATION bin)

configure_file(visualise.cmgui visualise.cmgui)

configure_file(inputs/mooney_rivlin.xml ${PROJECT_BINARY_DIR}/mooney_rivlin.xml)
