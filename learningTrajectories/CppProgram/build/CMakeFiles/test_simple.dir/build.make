# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/build

# Include any dependencies generated for this target.
include CMakeFiles/test_simple.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/test_simple.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/test_simple.dir/flags.make

CMakeFiles/test_simple.dir/test_simple.cpp.o: CMakeFiles/test_simple.dir/flags.make
CMakeFiles/test_simple.dir/test_simple.cpp.o: ../test_simple.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/test_simple.dir/test_simple.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/test_simple.dir/test_simple.cpp.o -c /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/test_simple.cpp

CMakeFiles/test_simple.dir/test_simple.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_simple.dir/test_simple.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/test_simple.cpp > CMakeFiles/test_simple.dir/test_simple.cpp.i

CMakeFiles/test_simple.dir/test_simple.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_simple.dir/test_simple.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/test_simple.cpp -o CMakeFiles/test_simple.dir/test_simple.cpp.s

CMakeFiles/test_simple.dir/test_simple.cpp.o.requires:
.PHONY : CMakeFiles/test_simple.dir/test_simple.cpp.o.requires

CMakeFiles/test_simple.dir/test_simple.cpp.o.provides: CMakeFiles/test_simple.dir/test_simple.cpp.o.requires
	$(MAKE) -f CMakeFiles/test_simple.dir/build.make CMakeFiles/test_simple.dir/test_simple.cpp.o.provides.build
.PHONY : CMakeFiles/test_simple.dir/test_simple.cpp.o.provides

CMakeFiles/test_simple.dir/test_simple.cpp.o.provides.build: CMakeFiles/test_simple.dir/test_simple.cpp.o

CMakeFiles/test_simple.dir/cartesianClient.cpp.o: CMakeFiles/test_simple.dir/flags.make
CMakeFiles/test_simple.dir/cartesianClient.cpp.o: ../cartesianClient.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/test_simple.dir/cartesianClient.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/test_simple.dir/cartesianClient.cpp.o -c /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/cartesianClient.cpp

CMakeFiles/test_simple.dir/cartesianClient.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_simple.dir/cartesianClient.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/cartesianClient.cpp > CMakeFiles/test_simple.dir/cartesianClient.cpp.i

CMakeFiles/test_simple.dir/cartesianClient.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_simple.dir/cartesianClient.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/cartesianClient.cpp -o CMakeFiles/test_simple.dir/cartesianClient.cpp.s

CMakeFiles/test_simple.dir/cartesianClient.cpp.o.requires:
.PHONY : CMakeFiles/test_simple.dir/cartesianClient.cpp.o.requires

CMakeFiles/test_simple.dir/cartesianClient.cpp.o.provides: CMakeFiles/test_simple.dir/cartesianClient.cpp.o.requires
	$(MAKE) -f CMakeFiles/test_simple.dir/build.make CMakeFiles/test_simple.dir/cartesianClient.cpp.o.provides.build
.PHONY : CMakeFiles/test_simple.dir/cartesianClient.cpp.o.provides

CMakeFiles/test_simple.dir/cartesianClient.cpp.o.provides.build: CMakeFiles/test_simple.dir/cartesianClient.cpp.o

# Object files for target test_simple
test_simple_OBJECTS = \
"CMakeFiles/test_simple.dir/test_simple.cpp.o" \
"CMakeFiles/test_simple.dir/cartesianClient.cpp.o"

# External object files for target test_simple
test_simple_EXTERNAL_OBJECTS =

bin/test_simple: CMakeFiles/test_simple.dir/test_simple.cpp.o
bin/test_simple: CMakeFiles/test_simple.dir/cartesianClient.cpp.o
bin/test_simple: CMakeFiles/test_simple.dir/build.make
bin/test_simple: /home/odermy/Software/src/yarp/build/lib/libYARP_OS.so
bin/test_simple: /home/odermy/Software/src/yarp/build/lib/libYARP_sig.so
bin/test_simple: /home/odermy/Software/src/yarp/build/lib/libYARP_math.so
bin/test_simple: /home/odermy/Software/src/yarp/build/lib/libYARP_dev.so
bin/test_simple: /home/odermy/Software/src/yarp/build/lib/libYARP_init.so
bin/test_simple: /home/odermy/Software/src/yarp/build/lib/libYARP_name.so
bin/test_simple: lib/libfake_geomagic_driver.a
bin/test_simple: /home/odermy/Software/src/yarp/build/lib/libYARP_math.so
bin/test_simple: /home/odermy/Software/src/yarp/build/lib/libYARP_dev.so
bin/test_simple: /home/odermy/Software/src/yarp/build/lib/libYARP_sig.so
bin/test_simple: /home/odermy/Software/src/yarp/build/lib/libYARP_init.so
bin/test_simple: /home/odermy/Software/src/yarp/build/lib/libYARP_name.so
bin/test_simple: /home/odermy/Software/src/yarp/build/lib/libYARP_OS.so
bin/test_simple: CMakeFiles/test_simple.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable bin/test_simple"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_simple.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/test_simple.dir/build: bin/test_simple
.PHONY : CMakeFiles/test_simple.dir/build

CMakeFiles/test_simple.dir/requires: CMakeFiles/test_simple.dir/test_simple.cpp.o.requires
CMakeFiles/test_simple.dir/requires: CMakeFiles/test_simple.dir/cartesianClient.cpp.o.requires
.PHONY : CMakeFiles/test_simple.dir/requires

CMakeFiles/test_simple.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/test_simple.dir/cmake_clean.cmake
.PHONY : CMakeFiles/test_simple.dir/clean

CMakeFiles/test_simple.dir/depend:
	cd /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/build /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/build /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/build/CMakeFiles/test_simple.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/test_simple.dir/depend
