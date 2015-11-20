
# Turn on more warnings
if(MSVC)
  # Force to always compile with W4
  if(CMAKE_C_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W4")
  endif()
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
endif()

include_directories("../include/")
<<<<<<< HEAD

add_library (optv SHARED tracking_frame_buf.c calibration.c parameters.c lsqadj.c ray_tracing.c trafo.c vec_utils.c image_processing.c multimed.c imgcoord.c epi.c orientation.c sortgrid.c segmentation.c)

=======
add_library (optv SHARED tracking_frame_buf.c calibration.c parameters.c lsqadj.c ray_tracing.c trafo.c vec_utils.c image_processing.c multimed.c imgcoord.c epi.c sortgrid.c correspondences.c btree.c)
>>>>>>> d810a88... updated the calls to the renewed versions of img_coord, pixel_to_metric from check_corresp


if(UNIX)
  target_link_libraries(optv m)
endif(UNIX)

install(TARGETS optv DESTINATION lib)
