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
include CMakeFiles/recordTrajectories.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/recordTrajectories.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/recordTrajectories.dir/flags.make

CMakeFiles/recordTrajectories.dir/recordTrajectories.cpp.o: CMakeFiles/recordTrajectories.dir/flags.make
CMakeFiles/recordTrajectories.dir/recordTrajectories.cpp.o: ../recordTrajectories.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/recordTrajectories.dir/recordTrajectories.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/recordTrajectories.dir/recordTrajectories.cpp.o -c /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/recordTrajectories.cpp

CMakeFiles/recordTrajectories.dir/recordTrajectories.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/recordTrajectories.dir/recordTrajectories.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/recordTrajectories.cpp > CMakeFiles/recordTrajectories.dir/recordTrajectories.cpp.i

CMakeFiles/recordTrajectories.dir/recordTrajectories.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/recordTrajectories.dir/recordTrajectories.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/recordTrajectories.cpp -o CMakeFiles/recordTrajectories.dir/recordTrajectories.cpp.s

CMakeFiles/recordTrajectories.dir/recordTrajectories.cpp.o.requires:
.PHONY : CMakeFiles/recordTrajectories.dir/recordTrajectories.cpp.o.requires

CMakeFiles/recordTrajectories.dir/recordTrajectories.cpp.o.provides: CMakeFiles/recordTrajectories.dir/recordTrajectories.cpp.o.requires
	$(MAKE) -f CMakeFiles/recordTrajectories.dir/build.make CMakeFiles/recordTrajectories.dir/recordTrajectories.cpp.o.provides.build
.PHONY : CMakeFiles/recordTrajectories.dir/recordTrajectories.cpp.o.provides

CMakeFiles/recordTrajectories.dir/recordTrajectories.cpp.o.provides.build: CMakeFiles/recordTrajectories.dir/recordTrajectories.cpp.o

CMakeFiles/recordTrajectories.dir/cartesianClient.cpp.o: CMakeFiles/recordTrajectories.dir/flags.make
CMakeFiles/recordTrajectories.dir/cartesianClient.cpp.o: ../cartesianClient.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/recordTrajectories.dir/cartesianClient.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/recordTrajectories.dir/cartesianClient.cpp.o -c /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/cartesianClient.cpp

CMakeFiles/recordTrajectories.dir/cartesianClient.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/recordTrajectories.dir/cartesianClient.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/cartesianClient.cpp > CMakeFiles/recordTrajectories.dir/cartesianClient.cpp.i

CMakeFiles/recordTrajectories.dir/cartesianClient.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/recordTrajectories.dir/cartesianClient.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/cartesianClient.cpp -o CMakeFiles/recordTrajectories.dir/cartesianClient.cpp.s

CMakeFiles/recordTrajectories.dir/cartesianClient.cpp.o.requires:
.PHONY : CMakeFiles/recordTrajectories.dir/cartesianClient.cpp.o.requires

CMakeFiles/recordTrajectories.dir/cartesianClient.cpp.o.provides: CMakeFiles/recordTrajectories.dir/cartesianClient.cpp.o.requires
	$(MAKE) -f CMakeFiles/recordTrajectories.dir/build.make CMakeFiles/recordTrajectories.dir/cartesianClient.cpp.o.provides.build
.PHONY : CMakeFiles/recordTrajectories.dir/cartesianClient.cpp.o.provides

CMakeFiles/recordTrajectories.dir/cartesianClient.cpp.o.provides.build: CMakeFiles/recordTrajectories.dir/cartesianClient.cpp.o

# Object files for target recordTrajectories
recordTrajectories_OBJECTS = \
"CMakeFiles/recordTrajectories.dir/recordTrajectories.cpp.o" \
"CMakeFiles/recordTrajectories.dir/cartesianClient.cpp.o"

# External object files for target recordTrajectories
recordTrajectories_EXTERNAL_OBJECTS =

bin/recordTrajectories: CMakeFiles/recordTrajectories.dir/recordTrajectories.cpp.o
bin/recordTrajectories: CMakeFiles/recordTrajectories.dir/cartesianClient.cpp.o
bin/recordTrajectories: CMakeFiles/recordTrajectories.dir/build.make
bin/recordTrajectories: /home/odermy/Software/src/yarp/build/lib/libYARP_OS.so
bin/recordTrajectories: /home/odermy/Software/src/yarp/build/lib/libYARP_sig.so
bin/recordTrajectories: /home/odermy/Software/src/yarp/build/lib/libYARP_math.so
bin/recordTrajectories: /home/odermy/Software/src/yarp/build/lib/libYARP_dev.so
bin/recordTrajectories: /home/odermy/Software/src/yarp/build/lib/libYARP_init.so
bin/recordTrajectories: /home/odermy/Software/src/yarp/build/lib/libYARP_name.so
bin/recordTrajectories: lib/libfake_geomagic_driver.a
bin/recordTrajectories: /home/odermy/Software/src/yarp/build/lib/libYARP_math.so
bin/recordTrajectories: /home/odermy/Software/src/yarp/build/lib/libYARP_dev.so
bin/recordTrajectories: /home/odermy/Software/src/yarp/build/lib/libYARP_sig.so
bin/recordTrajectories: /home/odermy/Software/src/yarp/build/lib/libYARP_init.so
bin/recordTrajectories: /home/odermy/Software/src/yarp/build/lib/libYARP_name.so
bin/recordTrajectories: /home/odermy/Software/src/yarp/build/lib/libYARP_OS.so
bin/recordTrajectories: CMakeFiles/recordTrajectories.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable bin/recordTrajectories"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/recordTrajectories.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/recordTrajectories.dir/build: bin/recordTrajectories
.PHONY : CMakeFiles/recordTrajectories.dir/build

CMakeFiles/recordTrajectories.dir/requires: CMakeFiles/recordTrajectories.dir/recordTrajectories.cpp.o.requires
CMakeFiles/recordTrajectories.dir/requires: CMakeFiles/recordTrajectories.dir/cartesianClient.cpp.o.requires
.PHONY : CMakeFiles/recordTrajectories.dir/requires

CMakeFiles/recordTrajectories.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/recordTrajectories.dir/cmake_clean.cmake
.PHONY : CMakeFiles/recordTrajectories.dir/clean

CMakeFiles/recordTrajectories.dir/depend:
	cd /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/build /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/build /home/odermy/Desktop/meeting_test/icub-learning-trajectories/learningTrajectories/CppProgram/build/CMakeFiles/recordTrajectories.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/recordTrajectories.dir/depend
