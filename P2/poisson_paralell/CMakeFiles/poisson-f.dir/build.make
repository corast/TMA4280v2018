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
CMAKE_SOURCE_DIR = /home/sondenge/git/TMA4280v2018/P2/poisson

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/sondenge/git/TMA4280v2018/P2/poisson

# Include any dependencies generated for this target.
include CMakeFiles/poisson-f.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/poisson-f.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/poisson-f.dir/flags.make

CMakeFiles/poisson-f.dir/poisson.f90.o: CMakeFiles/poisson-f.dir/flags.make
CMakeFiles/poisson-f.dir/poisson.f90.o: poisson.f90
	$(CMAKE_COMMAND) -E cmake_progress_report /home/sondenge/git/TMA4280v2018/P2/poisson/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object CMakeFiles/poisson-f.dir/poisson.f90.o"
	/share/apps/software/Compiler/GCC/6.4.0-2.28/OpenMPI/2.1.2/bin/mpif90  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /home/sondenge/git/TMA4280v2018/P2/poisson/poisson.f90 -o CMakeFiles/poisson-f.dir/poisson.f90.o

CMakeFiles/poisson-f.dir/poisson.f90.o.requires:
.PHONY : CMakeFiles/poisson-f.dir/poisson.f90.o.requires

CMakeFiles/poisson-f.dir/poisson.f90.o.provides: CMakeFiles/poisson-f.dir/poisson.f90.o.requires
	$(MAKE) -f CMakeFiles/poisson-f.dir/build.make CMakeFiles/poisson-f.dir/poisson.f90.o.provides.build
.PHONY : CMakeFiles/poisson-f.dir/poisson.f90.o.provides

CMakeFiles/poisson-f.dir/poisson.f90.o.provides.build: CMakeFiles/poisson-f.dir/poisson.f90.o

# Object files for target poisson-f
poisson__f_OBJECTS = \
"CMakeFiles/poisson-f.dir/poisson.f90.o"

# External object files for target poisson-f
poisson__f_EXTERNAL_OBJECTS =

poisson-f: CMakeFiles/poisson-f.dir/poisson.f90.o
poisson-f: CMakeFiles/poisson-f.dir/build.make
poisson-f: libcommon.a
poisson-f: /share/apps/software/Compiler/GCC/6.4.0-2.28/OpenMPI/2.1.2/lib/libmpi_usempif08.so
poisson-f: /share/apps/software/Compiler/GCC/6.4.0-2.28/OpenMPI/2.1.2/lib/libmpi_usempi_ignore_tkr.so
poisson-f: /share/apps/software/Compiler/GCC/6.4.0-2.28/OpenMPI/2.1.2/lib/libmpi_mpifh.so
poisson-f: /share/apps/software/Compiler/GCC/6.4.0-2.28/OpenMPI/2.1.2/lib/libmpi.so
poisson-f: CMakeFiles/poisson-f.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking Fortran executable poisson-f"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/poisson-f.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/poisson-f.dir/build: poisson-f
.PHONY : CMakeFiles/poisson-f.dir/build

CMakeFiles/poisson-f.dir/requires: CMakeFiles/poisson-f.dir/poisson.f90.o.requires
.PHONY : CMakeFiles/poisson-f.dir/requires

CMakeFiles/poisson-f.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/poisson-f.dir/cmake_clean.cmake
.PHONY : CMakeFiles/poisson-f.dir/clean

CMakeFiles/poisson-f.dir/depend:
	cd /home/sondenge/git/TMA4280v2018/P2/poisson && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/sondenge/git/TMA4280v2018/P2/poisson /home/sondenge/git/TMA4280v2018/P2/poisson /home/sondenge/git/TMA4280v2018/P2/poisson /home/sondenge/git/TMA4280v2018/P2/poisson /home/sondenge/git/TMA4280v2018/P2/poisson/CMakeFiles/poisson-f.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/poisson-f.dir/depend

