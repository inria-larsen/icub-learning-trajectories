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
CMAKE_SOURCE_DIR = /home/odermy/Desktop/learningTrajectories/CppProgram

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/odermy/Desktop/learningTrajectories/CppProgram/build

# Include any dependencies generated for this target.
include CMakeFiles/replayTrajectories.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/replayTrajectories.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/replayTrajectories.dir/flags.make

CMakeFiles/replayTrajectories.dir/replayTrajectories.cpp.o: CMakeFiles/replayTrajectories.dir/flags.make
CMakeFiles/replayTrajectories.dir/replayTrajectories.cpp.o: ../replayTrajectories.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/odermy/Desktop/learningTrajectories/CppProgram/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/replayTrajectories.dir/replayTrajectories.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/replayTrajectories.dir/replayTrajectories.cpp.o -c /home/odermy/Desktop/learningTrajectories/CppProgram/replayTrajectories.cpp

CMakeFiles/replayTrajectories.dir/replayTrajectories.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/replayTrajectories.dir/replayTrajectories.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/odermy/Desktop/learningTrajectories/CppProgram/replayTrajectories.cpp > CMakeFiles/replayTrajectories.dir/replayTrajectories.cpp.i

CMakeFiles/replayTrajectories.dir/replayTrajectories.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/replayTrajectories.dir/replayTrajectories.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/odermy/Desktop/learningTrajectories/CppProgram/replayTrajectories.cpp -o CMakeFiles/replayTrajectories.dir/replayTrajectories.cpp.s

CMakeFiles/replayTrajectories.dir/replayTrajectories.cpp.o.requires:
.PHONY : CMakeFiles/replayTrajectories.dir/replayTrajectories.cpp.o.requires

CMakeFiles/replayTrajectories.dir/replayTrajectories.cpp.o.provides: CMakeFiles/replayTrajectories.dir/replayTrajectories.cpp.o.requires
	$(MAKE) -f CMakeFiles/replayTrajectories.dir/build.make CMakeFiles/replayTrajectories.dir/replayTrajectories.cpp.o.provides.build
.PHONY : CMakeFiles/replayTrajectories.dir/replayTrajectories.cpp.o.provides

CMakeFiles/replayTrajectories.dir/replayTrajectories.cpp.o.provides.build: CMakeFiles/replayTrajectories.dir/replayTrajectories.cpp.o

CMakeFiles/replayTrajectories.dir/cartesianClient.cpp.o: CMakeFiles/replayTrajectories.dir/flags.make
CMakeFiles/replayTrajectories.dir/cartesianClient.cpp.o: ../cartesianClient.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/odermy/Desktop/learningTrajectories/CppProgram/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/replayTrajectories.dir/cartesianClient.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/replayTrajectories.dir/cartesianClient.cpp.o -c /home/odermy/Desktop/learningTrajectories/CppProgram/cartesianClient.cpp

CMakeFiles/replayTrajectories.dir/cartesianClient.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/replayTrajectories.dir/cartesianClient.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/odermy/Desktop/learningTrajectories/CppProgram/cartesianClient.cpp > CMakeFiles/replayTrajectories.dir/cartesianClient.cpp.i

CMakeFiles/replayTrajectories.dir/cartesianClient.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/replayTrajectories.dir/cartesianClient.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/odermy/Desktop/learningTrajectories/CppProgram/cartesianClient.cpp -o CMakeFiles/replayTrajectories.dir/cartesianClient.cpp.s

CMakeFiles/replayTrajectories.dir/cartesianClient.cpp.o.requires:
.PHONY : CMakeFiles/replayTrajectories.dir/cartesianClient.cpp.o.requires

CMakeFiles/replayTrajectories.dir/cartesianClient.cpp.o.provides: CMakeFiles/replayTrajectories.dir/cartesianClient.cpp.o.requires
	$(MAKE) -f CMakeFiles/replayTrajectories.dir/build.make CMakeFiles/replayTrajectories.dir/cartesianClient.cpp.o.provides.build
.PHONY : CMakeFiles/replayTrajectories.dir/cartesianClient.cpp.o.provides

CMakeFiles/replayTrajectories.dir/cartesianClient.cpp.o.provides.build: CMakeFiles/replayTrajectories.dir/cartesianClient.cpp.o

# Object files for target replayTrajectories
replayTrajectories_OBJECTS = \
"CMakeFiles/replayTrajectories.dir/replayTrajectories.cpp.o" \
"CMakeFiles/replayTrajectories.dir/cartesianClient.cpp.o"

# External object files for target replayTrajectories
replayTrajectories_EXTERNAL_OBJECTS =

bin/replayTrajectories: CMakeFiles/replayTrajectories.dir/replayTrajectories.cpp.o
bin/replayTrajectories: CMakeFiles/replayTrajectories.dir/cartesianClient.cpp.o
bin/replayTrajectories: CMakeFiles/replayTrajectories.dir/build.make
bin/replayTrajectories: /home/odermy/Software/src/yarp/build/lib/libYARP_OS.so
bin/replayTrajectories: /home/odermy/Software/src/yarp/build/lib/libYARP_sig.so
bin/replayTrajectories: /home/odermy/Software/src/yarp/build/lib/libYARP_math.so
bin/replayTrajectories: /home/odermy/Software/src/yarp/build/lib/libYARP_dev.so
bin/replayTrajectories: /home/odermy/Software/src/yarp/build/lib/libYARP_init.so
bin/replayTrajectories: /home/odermy/Software/src/yarp/build/lib/libYARP_name.so
bin/replayTrajectories: /home/odermy/Software/src/yarp/build/lib/libYARP_sig.so
bin/replayTrajectories: /home/odermy/Software/src/yarp/build/lib/libYARP_OS.so
bin/replayTrajectories: CMakeFiles/replayTrajectories.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable bin/replayTrajectories"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/replayTrajectories.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/replayTrajectories.dir/build: bin/replayTrajectories
.PHONY : CMakeFiles/replayTrajectories.dir/build

CMakeFiles/replayTrajectories.dir/requires: CMakeFiles/replayTrajectories.dir/replayTrajectories.cpp.o.requires
CMakeFiles/replayTrajectories.dir/requires: CMakeFiles/replayTrajectories.dir/cartesianClient.cpp.o.requires
.PHONY : CMakeFiles/replayTrajectories.dir/requires

CMakeFiles/replayTrajectories.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/replayTrajectories.dir/cmake_clean.cmake
.PHONY : CMakeFiles/replayTrajectories.dir/clean

CMakeFiles/replayTrajectories.dir/depend:
	cd /home/odermy/Desktop/learningTrajectories/CppProgram/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/odermy/Desktop/learningTrajectories/CppProgram /home/odermy/Desktop/learningTrajectories/CppProgram /home/odermy/Desktop/learningTrajectories/CppProgram/build /home/odermy/Desktop/learningTrajectories/CppProgram/build /home/odermy/Desktop/learningTrajectories/CppProgram/build/CMakeFiles/replayTrajectories.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/replayTrajectories.dir/depend

