# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/felix/CLionProjects/EPITA-BioInfo-2023/projet_bioinfo_wirth

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/felix/CLionProjects/EPITA-BioInfo-2023/projet_bioinfo_wirth/build

# Include any dependencies generated for this target.
include CMakeFiles/HMM-align.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/HMM-align.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/HMM-align.dir/flags.make

CMakeFiles/HMM-align.dir/src/hmm-align/main.cpp.o: CMakeFiles/HMM-align.dir/flags.make
CMakeFiles/HMM-align.dir/src/hmm-align/main.cpp.o: ../src/hmm-align/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/felix/CLionProjects/EPITA-BioInfo-2023/projet_bioinfo_wirth/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/HMM-align.dir/src/hmm-align/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/HMM-align.dir/src/hmm-align/main.cpp.o -c /home/felix/CLionProjects/EPITA-BioInfo-2023/projet_bioinfo_wirth/src/hmm-align/main.cpp

CMakeFiles/HMM-align.dir/src/hmm-align/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/HMM-align.dir/src/hmm-align/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/felix/CLionProjects/EPITA-BioInfo-2023/projet_bioinfo_wirth/src/hmm-align/main.cpp > CMakeFiles/HMM-align.dir/src/hmm-align/main.cpp.i

CMakeFiles/HMM-align.dir/src/hmm-align/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/HMM-align.dir/src/hmm-align/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/felix/CLionProjects/EPITA-BioInfo-2023/projet_bioinfo_wirth/src/hmm-align/main.cpp -o CMakeFiles/HMM-align.dir/src/hmm-align/main.cpp.s

CMakeFiles/HMM-align.dir/src/HMM/Fasta.cpp.o: CMakeFiles/HMM-align.dir/flags.make
CMakeFiles/HMM-align.dir/src/HMM/Fasta.cpp.o: ../src/HMM/Fasta.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/felix/CLionProjects/EPITA-BioInfo-2023/projet_bioinfo_wirth/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/HMM-align.dir/src/HMM/Fasta.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/HMM-align.dir/src/HMM/Fasta.cpp.o -c /home/felix/CLionProjects/EPITA-BioInfo-2023/projet_bioinfo_wirth/src/HMM/Fasta.cpp

CMakeFiles/HMM-align.dir/src/HMM/Fasta.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/HMM-align.dir/src/HMM/Fasta.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/felix/CLionProjects/EPITA-BioInfo-2023/projet_bioinfo_wirth/src/HMM/Fasta.cpp > CMakeFiles/HMM-align.dir/src/HMM/Fasta.cpp.i

CMakeFiles/HMM-align.dir/src/HMM/Fasta.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/HMM-align.dir/src/HMM/Fasta.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/felix/CLionProjects/EPITA-BioInfo-2023/projet_bioinfo_wirth/src/HMM/Fasta.cpp -o CMakeFiles/HMM-align.dir/src/HMM/Fasta.cpp.s

CMakeFiles/HMM-align.dir/src/HMM/HMM.cpp.o: CMakeFiles/HMM-align.dir/flags.make
CMakeFiles/HMM-align.dir/src/HMM/HMM.cpp.o: ../src/HMM/HMM.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/felix/CLionProjects/EPITA-BioInfo-2023/projet_bioinfo_wirth/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/HMM-align.dir/src/HMM/HMM.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/HMM-align.dir/src/HMM/HMM.cpp.o -c /home/felix/CLionProjects/EPITA-BioInfo-2023/projet_bioinfo_wirth/src/HMM/HMM.cpp

CMakeFiles/HMM-align.dir/src/HMM/HMM.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/HMM-align.dir/src/HMM/HMM.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/felix/CLionProjects/EPITA-BioInfo-2023/projet_bioinfo_wirth/src/HMM/HMM.cpp > CMakeFiles/HMM-align.dir/src/HMM/HMM.cpp.i

CMakeFiles/HMM-align.dir/src/HMM/HMM.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/HMM-align.dir/src/HMM/HMM.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/felix/CLionProjects/EPITA-BioInfo-2023/projet_bioinfo_wirth/src/HMM/HMM.cpp -o CMakeFiles/HMM-align.dir/src/HMM/HMM.cpp.s

# Object files for target HMM-align
HMM__align_OBJECTS = \
"CMakeFiles/HMM-align.dir/src/hmm-align/main.cpp.o" \
"CMakeFiles/HMM-align.dir/src/HMM/Fasta.cpp.o" \
"CMakeFiles/HMM-align.dir/src/HMM/HMM.cpp.o"

# External object files for target HMM-align
HMM__align_EXTERNAL_OBJECTS =

HMM-align: CMakeFiles/HMM-align.dir/src/hmm-align/main.cpp.o
HMM-align: CMakeFiles/HMM-align.dir/src/HMM/Fasta.cpp.o
HMM-align: CMakeFiles/HMM-align.dir/src/HMM/HMM.cpp.o
HMM-align: CMakeFiles/HMM-align.dir/build.make
HMM-align: CMakeFiles/HMM-align.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/felix/CLionProjects/EPITA-BioInfo-2023/projet_bioinfo_wirth/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable HMM-align"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/HMM-align.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/HMM-align.dir/build: HMM-align

.PHONY : CMakeFiles/HMM-align.dir/build

CMakeFiles/HMM-align.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/HMM-align.dir/cmake_clean.cmake
.PHONY : CMakeFiles/HMM-align.dir/clean

CMakeFiles/HMM-align.dir/depend:
	cd /home/felix/CLionProjects/EPITA-BioInfo-2023/projet_bioinfo_wirth/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/felix/CLionProjects/EPITA-BioInfo-2023/projet_bioinfo_wirth /home/felix/CLionProjects/EPITA-BioInfo-2023/projet_bioinfo_wirth /home/felix/CLionProjects/EPITA-BioInfo-2023/projet_bioinfo_wirth/build /home/felix/CLionProjects/EPITA-BioInfo-2023/projet_bioinfo_wirth/build /home/felix/CLionProjects/EPITA-BioInfo-2023/projet_bioinfo_wirth/build/CMakeFiles/HMM-align.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/HMM-align.dir/depend

