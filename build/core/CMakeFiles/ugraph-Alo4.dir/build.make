# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/gflychee/vldb2018/kmpc

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/gflychee/vldb2018/kmpc/build

# Include any dependencies generated for this target.
include core/CMakeFiles/ugraph-Alo4.dir/depend.make

# Include the progress variables for this target.
include core/CMakeFiles/ugraph-Alo4.dir/progress.make

# Include the compile flags for this target's objects.
include core/CMakeFiles/ugraph-Alo4.dir/flags.make

core/CMakeFiles/ugraph-Alo4.dir/Alo4.cpp.o: core/CMakeFiles/ugraph-Alo4.dir/flags.make
core/CMakeFiles/ugraph-Alo4.dir/Alo4.cpp.o: ../core/Alo4.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/gflychee/vldb2018/kmpc/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object core/CMakeFiles/ugraph-Alo4.dir/Alo4.cpp.o"
	cd /home/gflychee/vldb2018/kmpc/build/core && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ugraph-Alo4.dir/Alo4.cpp.o -c /home/gflychee/vldb2018/kmpc/core/Alo4.cpp

core/CMakeFiles/ugraph-Alo4.dir/Alo4.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ugraph-Alo4.dir/Alo4.cpp.i"
	cd /home/gflychee/vldb2018/kmpc/build/core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/gflychee/vldb2018/kmpc/core/Alo4.cpp > CMakeFiles/ugraph-Alo4.dir/Alo4.cpp.i

core/CMakeFiles/ugraph-Alo4.dir/Alo4.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ugraph-Alo4.dir/Alo4.cpp.s"
	cd /home/gflychee/vldb2018/kmpc/build/core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/gflychee/vldb2018/kmpc/core/Alo4.cpp -o CMakeFiles/ugraph-Alo4.dir/Alo4.cpp.s

core/CMakeFiles/ugraph-Alo4.dir/Alo4.cpp.o.requires:

.PHONY : core/CMakeFiles/ugraph-Alo4.dir/Alo4.cpp.o.requires

core/CMakeFiles/ugraph-Alo4.dir/Alo4.cpp.o.provides: core/CMakeFiles/ugraph-Alo4.dir/Alo4.cpp.o.requires
	$(MAKE) -f core/CMakeFiles/ugraph-Alo4.dir/build.make core/CMakeFiles/ugraph-Alo4.dir/Alo4.cpp.o.provides.build
.PHONY : core/CMakeFiles/ugraph-Alo4.dir/Alo4.cpp.o.provides

core/CMakeFiles/ugraph-Alo4.dir/Alo4.cpp.o.provides.build: core/CMakeFiles/ugraph-Alo4.dir/Alo4.cpp.o


# Object files for target ugraph-Alo4
ugraph__Alo4_OBJECTS = \
"CMakeFiles/ugraph-Alo4.dir/Alo4.cpp.o"

# External object files for target ugraph-Alo4
ugraph__Alo4_EXTERNAL_OBJECTS =

core/ugraph-Alo4: core/CMakeFiles/ugraph-Alo4.dir/Alo4.cpp.o
core/ugraph-Alo4: core/CMakeFiles/ugraph-Alo4.dir/build.make
core/ugraph-Alo4: core/libugraph.a
core/ugraph-Alo4: /usr/lib/libboost_system.so
core/ugraph-Alo4: /usr/lib/libboost_date_time.so
core/ugraph-Alo4: /usr/lib/libboost_program_options.so
core/ugraph-Alo4: /usr/lib/libboost_iostreams.so
core/ugraph-Alo4: /usr/lib/x86_64-linux-gnu/libbz2.so
core/ugraph-Alo4: core/CMakeFiles/ugraph-Alo4.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/gflychee/vldb2018/kmpc/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ugraph-Alo4"
	cd /home/gflychee/vldb2018/kmpc/build/core && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ugraph-Alo4.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
core/CMakeFiles/ugraph-Alo4.dir/build: core/ugraph-Alo4

.PHONY : core/CMakeFiles/ugraph-Alo4.dir/build

core/CMakeFiles/ugraph-Alo4.dir/requires: core/CMakeFiles/ugraph-Alo4.dir/Alo4.cpp.o.requires

.PHONY : core/CMakeFiles/ugraph-Alo4.dir/requires

core/CMakeFiles/ugraph-Alo4.dir/clean:
	cd /home/gflychee/vldb2018/kmpc/build/core && $(CMAKE_COMMAND) -P CMakeFiles/ugraph-Alo4.dir/cmake_clean.cmake
.PHONY : core/CMakeFiles/ugraph-Alo4.dir/clean

core/CMakeFiles/ugraph-Alo4.dir/depend:
	cd /home/gflychee/vldb2018/kmpc/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/gflychee/vldb2018/kmpc /home/gflychee/vldb2018/kmpc/core /home/gflychee/vldb2018/kmpc/build /home/gflychee/vldb2018/kmpc/build/core /home/gflychee/vldb2018/kmpc/build/core/CMakeFiles/ugraph-Alo4.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : core/CMakeFiles/ugraph-Alo4.dir/depend

