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
CMAKE_SOURCE_DIR = /tmp/stuff/masbcpp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /tmp/stuff/masbcpp

# Include any dependencies generated for this target.
include CMakeFiles/compute_ma.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/compute_ma.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/compute_ma.dir/flags.make

CMakeFiles/compute_ma.dir/src/compute_ma.cpp.o: CMakeFiles/compute_ma.dir/flags.make
CMakeFiles/compute_ma.dir/src/compute_ma.cpp.o: src/compute_ma.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/tmp/stuff/masbcpp/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/compute_ma.dir/src/compute_ma.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/compute_ma.dir/src/compute_ma.cpp.o -c /tmp/stuff/masbcpp/src/compute_ma.cpp

CMakeFiles/compute_ma.dir/src/compute_ma.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/compute_ma.dir/src/compute_ma.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /tmp/stuff/masbcpp/src/compute_ma.cpp > CMakeFiles/compute_ma.dir/src/compute_ma.cpp.i

CMakeFiles/compute_ma.dir/src/compute_ma.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/compute_ma.dir/src/compute_ma.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /tmp/stuff/masbcpp/src/compute_ma.cpp -o CMakeFiles/compute_ma.dir/src/compute_ma.cpp.s

CMakeFiles/compute_ma.dir/src/compute_ma.cpp.o.requires:

.PHONY : CMakeFiles/compute_ma.dir/src/compute_ma.cpp.o.requires

CMakeFiles/compute_ma.dir/src/compute_ma.cpp.o.provides: CMakeFiles/compute_ma.dir/src/compute_ma.cpp.o.requires
	$(MAKE) -f CMakeFiles/compute_ma.dir/build.make CMakeFiles/compute_ma.dir/src/compute_ma.cpp.o.provides.build
.PHONY : CMakeFiles/compute_ma.dir/src/compute_ma.cpp.o.provides

CMakeFiles/compute_ma.dir/src/compute_ma.cpp.o.provides.build: CMakeFiles/compute_ma.dir/src/compute_ma.cpp.o


# Object files for target compute_ma
compute_ma_OBJECTS = \
"CMakeFiles/compute_ma.dir/src/compute_ma.cpp.o"

# External object files for target compute_ma
compute_ma_EXTERNAL_OBJECTS =

compute_ma: CMakeFiles/compute_ma.dir/src/compute_ma.cpp.o
compute_ma: CMakeFiles/compute_ma.dir/build.make
compute_ma: libmasbcpp.a
compute_ma: /usr/lib/libpcl_features.so
compute_ma: /usr/lib/libpcl_filters.so
compute_ma: /usr/lib/libpcl_sample_consensus.so
compute_ma: /usr/lib/libpcl_search.so
compute_ma: /usr/lib/libpcl_kdtree.so
compute_ma: /usr/lib/libpcl_octree.so
compute_ma: /usr/lib/libpcl_common.so
compute_ma: CMakeFiles/compute_ma.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/tmp/stuff/masbcpp/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable compute_ma"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/compute_ma.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/compute_ma.dir/build: compute_ma

.PHONY : CMakeFiles/compute_ma.dir/build

CMakeFiles/compute_ma.dir/requires: CMakeFiles/compute_ma.dir/src/compute_ma.cpp.o.requires

.PHONY : CMakeFiles/compute_ma.dir/requires

CMakeFiles/compute_ma.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/compute_ma.dir/cmake_clean.cmake
.PHONY : CMakeFiles/compute_ma.dir/clean

CMakeFiles/compute_ma.dir/depend:
	cd /tmp/stuff/masbcpp && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /tmp/stuff/masbcpp /tmp/stuff/masbcpp /tmp/stuff/masbcpp /tmp/stuff/masbcpp /tmp/stuff/masbcpp/CMakeFiles/compute_ma.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/compute_ma.dir/depend

