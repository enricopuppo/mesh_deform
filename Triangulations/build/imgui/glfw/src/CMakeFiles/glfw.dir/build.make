# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.26.3/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.26.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build

# Include any dependencies generated for this target.
include imgui/glfw/src/CMakeFiles/glfw.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include imgui/glfw/src/CMakeFiles/glfw.dir/compiler_depend.make

# Include the progress variables for this target.
include imgui/glfw/src/CMakeFiles/glfw.dir/progress.make

# Include the compile flags for this target's objects.
include imgui/glfw/src/CMakeFiles/glfw.dir/flags.make

imgui/glfw/src/CMakeFiles/glfw.dir/context.c.o: imgui/glfw/src/CMakeFiles/glfw.dir/flags.make
imgui/glfw/src/CMakeFiles/glfw.dir/context.c.o: /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/context.c
imgui/glfw/src/CMakeFiles/glfw.dir/context.c.o: imgui/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object imgui/glfw/src/CMakeFiles/glfw.dir/context.c.o"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT imgui/glfw/src/CMakeFiles/glfw.dir/context.c.o -MF CMakeFiles/glfw.dir/context.c.o.d -o CMakeFiles/glfw.dir/context.c.o -c /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/context.c

imgui/glfw/src/CMakeFiles/glfw.dir/context.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/context.c.i"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/context.c > CMakeFiles/glfw.dir/context.c.i

imgui/glfw/src/CMakeFiles/glfw.dir/context.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/context.c.s"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/context.c -o CMakeFiles/glfw.dir/context.c.s

imgui/glfw/src/CMakeFiles/glfw.dir/init.c.o: imgui/glfw/src/CMakeFiles/glfw.dir/flags.make
imgui/glfw/src/CMakeFiles/glfw.dir/init.c.o: /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/init.c
imgui/glfw/src/CMakeFiles/glfw.dir/init.c.o: imgui/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object imgui/glfw/src/CMakeFiles/glfw.dir/init.c.o"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT imgui/glfw/src/CMakeFiles/glfw.dir/init.c.o -MF CMakeFiles/glfw.dir/init.c.o.d -o CMakeFiles/glfw.dir/init.c.o -c /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/init.c

imgui/glfw/src/CMakeFiles/glfw.dir/init.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/init.c.i"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/init.c > CMakeFiles/glfw.dir/init.c.i

imgui/glfw/src/CMakeFiles/glfw.dir/init.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/init.c.s"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/init.c -o CMakeFiles/glfw.dir/init.c.s

imgui/glfw/src/CMakeFiles/glfw.dir/input.c.o: imgui/glfw/src/CMakeFiles/glfw.dir/flags.make
imgui/glfw/src/CMakeFiles/glfw.dir/input.c.o: /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/input.c
imgui/glfw/src/CMakeFiles/glfw.dir/input.c.o: imgui/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object imgui/glfw/src/CMakeFiles/glfw.dir/input.c.o"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT imgui/glfw/src/CMakeFiles/glfw.dir/input.c.o -MF CMakeFiles/glfw.dir/input.c.o.d -o CMakeFiles/glfw.dir/input.c.o -c /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/input.c

imgui/glfw/src/CMakeFiles/glfw.dir/input.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/input.c.i"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/input.c > CMakeFiles/glfw.dir/input.c.i

imgui/glfw/src/CMakeFiles/glfw.dir/input.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/input.c.s"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/input.c -o CMakeFiles/glfw.dir/input.c.s

imgui/glfw/src/CMakeFiles/glfw.dir/monitor.c.o: imgui/glfw/src/CMakeFiles/glfw.dir/flags.make
imgui/glfw/src/CMakeFiles/glfw.dir/monitor.c.o: /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/monitor.c
imgui/glfw/src/CMakeFiles/glfw.dir/monitor.c.o: imgui/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object imgui/glfw/src/CMakeFiles/glfw.dir/monitor.c.o"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT imgui/glfw/src/CMakeFiles/glfw.dir/monitor.c.o -MF CMakeFiles/glfw.dir/monitor.c.o.d -o CMakeFiles/glfw.dir/monitor.c.o -c /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/monitor.c

imgui/glfw/src/CMakeFiles/glfw.dir/monitor.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/monitor.c.i"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/monitor.c > CMakeFiles/glfw.dir/monitor.c.i

imgui/glfw/src/CMakeFiles/glfw.dir/monitor.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/monitor.c.s"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/monitor.c -o CMakeFiles/glfw.dir/monitor.c.s

imgui/glfw/src/CMakeFiles/glfw.dir/vulkan.c.o: imgui/glfw/src/CMakeFiles/glfw.dir/flags.make
imgui/glfw/src/CMakeFiles/glfw.dir/vulkan.c.o: /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/vulkan.c
imgui/glfw/src/CMakeFiles/glfw.dir/vulkan.c.o: imgui/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object imgui/glfw/src/CMakeFiles/glfw.dir/vulkan.c.o"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT imgui/glfw/src/CMakeFiles/glfw.dir/vulkan.c.o -MF CMakeFiles/glfw.dir/vulkan.c.o.d -o CMakeFiles/glfw.dir/vulkan.c.o -c /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/vulkan.c

imgui/glfw/src/CMakeFiles/glfw.dir/vulkan.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/vulkan.c.i"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/vulkan.c > CMakeFiles/glfw.dir/vulkan.c.i

imgui/glfw/src/CMakeFiles/glfw.dir/vulkan.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/vulkan.c.s"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/vulkan.c -o CMakeFiles/glfw.dir/vulkan.c.s

imgui/glfw/src/CMakeFiles/glfw.dir/window.c.o: imgui/glfw/src/CMakeFiles/glfw.dir/flags.make
imgui/glfw/src/CMakeFiles/glfw.dir/window.c.o: /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/window.c
imgui/glfw/src/CMakeFiles/glfw.dir/window.c.o: imgui/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object imgui/glfw/src/CMakeFiles/glfw.dir/window.c.o"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT imgui/glfw/src/CMakeFiles/glfw.dir/window.c.o -MF CMakeFiles/glfw.dir/window.c.o.d -o CMakeFiles/glfw.dir/window.c.o -c /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/window.c

imgui/glfw/src/CMakeFiles/glfw.dir/window.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/window.c.i"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/window.c > CMakeFiles/glfw.dir/window.c.i

imgui/glfw/src/CMakeFiles/glfw.dir/window.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/window.c.s"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/window.c -o CMakeFiles/glfw.dir/window.c.s

imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_init.m.o: imgui/glfw/src/CMakeFiles/glfw.dir/flags.make
imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_init.m.o: /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/cocoa_init.m
imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_init.m.o: imgui/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building C object imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_init.m.o"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_init.m.o -MF CMakeFiles/glfw.dir/cocoa_init.m.o.d -o CMakeFiles/glfw.dir/cocoa_init.m.o -c /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/cocoa_init.m

imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_init.m.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/cocoa_init.m.i"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/cocoa_init.m > CMakeFiles/glfw.dir/cocoa_init.m.i

imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_init.m.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/cocoa_init.m.s"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/cocoa_init.m -o CMakeFiles/glfw.dir/cocoa_init.m.s

imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_joystick.m.o: imgui/glfw/src/CMakeFiles/glfw.dir/flags.make
imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_joystick.m.o: /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/cocoa_joystick.m
imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_joystick.m.o: imgui/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building C object imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_joystick.m.o"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_joystick.m.o -MF CMakeFiles/glfw.dir/cocoa_joystick.m.o.d -o CMakeFiles/glfw.dir/cocoa_joystick.m.o -c /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/cocoa_joystick.m

imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_joystick.m.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/cocoa_joystick.m.i"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/cocoa_joystick.m > CMakeFiles/glfw.dir/cocoa_joystick.m.i

imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_joystick.m.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/cocoa_joystick.m.s"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/cocoa_joystick.m -o CMakeFiles/glfw.dir/cocoa_joystick.m.s

imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_monitor.m.o: imgui/glfw/src/CMakeFiles/glfw.dir/flags.make
imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_monitor.m.o: /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/cocoa_monitor.m
imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_monitor.m.o: imgui/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building C object imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_monitor.m.o"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_monitor.m.o -MF CMakeFiles/glfw.dir/cocoa_monitor.m.o.d -o CMakeFiles/glfw.dir/cocoa_monitor.m.o -c /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/cocoa_monitor.m

imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_monitor.m.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/cocoa_monitor.m.i"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/cocoa_monitor.m > CMakeFiles/glfw.dir/cocoa_monitor.m.i

imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_monitor.m.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/cocoa_monitor.m.s"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/cocoa_monitor.m -o CMakeFiles/glfw.dir/cocoa_monitor.m.s

imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_window.m.o: imgui/glfw/src/CMakeFiles/glfw.dir/flags.make
imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_window.m.o: /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/cocoa_window.m
imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_window.m.o: imgui/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building C object imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_window.m.o"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_window.m.o -MF CMakeFiles/glfw.dir/cocoa_window.m.o.d -o CMakeFiles/glfw.dir/cocoa_window.m.o -c /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/cocoa_window.m

imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_window.m.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/cocoa_window.m.i"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/cocoa_window.m > CMakeFiles/glfw.dir/cocoa_window.m.i

imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_window.m.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/cocoa_window.m.s"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/cocoa_window.m -o CMakeFiles/glfw.dir/cocoa_window.m.s

imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_time.c.o: imgui/glfw/src/CMakeFiles/glfw.dir/flags.make
imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_time.c.o: /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/cocoa_time.c
imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_time.c.o: imgui/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building C object imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_time.c.o"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_time.c.o -MF CMakeFiles/glfw.dir/cocoa_time.c.o.d -o CMakeFiles/glfw.dir/cocoa_time.c.o -c /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/cocoa_time.c

imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_time.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/cocoa_time.c.i"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/cocoa_time.c > CMakeFiles/glfw.dir/cocoa_time.c.i

imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_time.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/cocoa_time.c.s"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/cocoa_time.c -o CMakeFiles/glfw.dir/cocoa_time.c.s

imgui/glfw/src/CMakeFiles/glfw.dir/posix_thread.c.o: imgui/glfw/src/CMakeFiles/glfw.dir/flags.make
imgui/glfw/src/CMakeFiles/glfw.dir/posix_thread.c.o: /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/posix_thread.c
imgui/glfw/src/CMakeFiles/glfw.dir/posix_thread.c.o: imgui/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building C object imgui/glfw/src/CMakeFiles/glfw.dir/posix_thread.c.o"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT imgui/glfw/src/CMakeFiles/glfw.dir/posix_thread.c.o -MF CMakeFiles/glfw.dir/posix_thread.c.o.d -o CMakeFiles/glfw.dir/posix_thread.c.o -c /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/posix_thread.c

imgui/glfw/src/CMakeFiles/glfw.dir/posix_thread.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/posix_thread.c.i"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/posix_thread.c > CMakeFiles/glfw.dir/posix_thread.c.i

imgui/glfw/src/CMakeFiles/glfw.dir/posix_thread.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/posix_thread.c.s"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/posix_thread.c -o CMakeFiles/glfw.dir/posix_thread.c.s

imgui/glfw/src/CMakeFiles/glfw.dir/nsgl_context.m.o: imgui/glfw/src/CMakeFiles/glfw.dir/flags.make
imgui/glfw/src/CMakeFiles/glfw.dir/nsgl_context.m.o: /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/nsgl_context.m
imgui/glfw/src/CMakeFiles/glfw.dir/nsgl_context.m.o: imgui/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building C object imgui/glfw/src/CMakeFiles/glfw.dir/nsgl_context.m.o"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT imgui/glfw/src/CMakeFiles/glfw.dir/nsgl_context.m.o -MF CMakeFiles/glfw.dir/nsgl_context.m.o.d -o CMakeFiles/glfw.dir/nsgl_context.m.o -c /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/nsgl_context.m

imgui/glfw/src/CMakeFiles/glfw.dir/nsgl_context.m.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/nsgl_context.m.i"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/nsgl_context.m > CMakeFiles/glfw.dir/nsgl_context.m.i

imgui/glfw/src/CMakeFiles/glfw.dir/nsgl_context.m.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/nsgl_context.m.s"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/nsgl_context.m -o CMakeFiles/glfw.dir/nsgl_context.m.s

imgui/glfw/src/CMakeFiles/glfw.dir/egl_context.c.o: imgui/glfw/src/CMakeFiles/glfw.dir/flags.make
imgui/glfw/src/CMakeFiles/glfw.dir/egl_context.c.o: /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/egl_context.c
imgui/glfw/src/CMakeFiles/glfw.dir/egl_context.c.o: imgui/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building C object imgui/glfw/src/CMakeFiles/glfw.dir/egl_context.c.o"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT imgui/glfw/src/CMakeFiles/glfw.dir/egl_context.c.o -MF CMakeFiles/glfw.dir/egl_context.c.o.d -o CMakeFiles/glfw.dir/egl_context.c.o -c /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/egl_context.c

imgui/glfw/src/CMakeFiles/glfw.dir/egl_context.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/egl_context.c.i"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/egl_context.c > CMakeFiles/glfw.dir/egl_context.c.i

imgui/glfw/src/CMakeFiles/glfw.dir/egl_context.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/egl_context.c.s"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/egl_context.c -o CMakeFiles/glfw.dir/egl_context.c.s

imgui/glfw/src/CMakeFiles/glfw.dir/osmesa_context.c.o: imgui/glfw/src/CMakeFiles/glfw.dir/flags.make
imgui/glfw/src/CMakeFiles/glfw.dir/osmesa_context.c.o: /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/osmesa_context.c
imgui/glfw/src/CMakeFiles/glfw.dir/osmesa_context.c.o: imgui/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Building C object imgui/glfw/src/CMakeFiles/glfw.dir/osmesa_context.c.o"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT imgui/glfw/src/CMakeFiles/glfw.dir/osmesa_context.c.o -MF CMakeFiles/glfw.dir/osmesa_context.c.o.d -o CMakeFiles/glfw.dir/osmesa_context.c.o -c /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/osmesa_context.c

imgui/glfw/src/CMakeFiles/glfw.dir/osmesa_context.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/osmesa_context.c.i"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/osmesa_context.c > CMakeFiles/glfw.dir/osmesa_context.c.i

imgui/glfw/src/CMakeFiles/glfw.dir/osmesa_context.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/osmesa_context.c.s"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src/osmesa_context.c -o CMakeFiles/glfw.dir/osmesa_context.c.s

# Object files for target glfw
glfw_OBJECTS = \
"CMakeFiles/glfw.dir/context.c.o" \
"CMakeFiles/glfw.dir/init.c.o" \
"CMakeFiles/glfw.dir/input.c.o" \
"CMakeFiles/glfw.dir/monitor.c.o" \
"CMakeFiles/glfw.dir/vulkan.c.o" \
"CMakeFiles/glfw.dir/window.c.o" \
"CMakeFiles/glfw.dir/cocoa_init.m.o" \
"CMakeFiles/glfw.dir/cocoa_joystick.m.o" \
"CMakeFiles/glfw.dir/cocoa_monitor.m.o" \
"CMakeFiles/glfw.dir/cocoa_window.m.o" \
"CMakeFiles/glfw.dir/cocoa_time.c.o" \
"CMakeFiles/glfw.dir/posix_thread.c.o" \
"CMakeFiles/glfw.dir/nsgl_context.m.o" \
"CMakeFiles/glfw.dir/egl_context.c.o" \
"CMakeFiles/glfw.dir/osmesa_context.c.o"

# External object files for target glfw
glfw_EXTERNAL_OBJECTS =

imgui/glfw/src/libglfw3.a: imgui/glfw/src/CMakeFiles/glfw.dir/context.c.o
imgui/glfw/src/libglfw3.a: imgui/glfw/src/CMakeFiles/glfw.dir/init.c.o
imgui/glfw/src/libglfw3.a: imgui/glfw/src/CMakeFiles/glfw.dir/input.c.o
imgui/glfw/src/libglfw3.a: imgui/glfw/src/CMakeFiles/glfw.dir/monitor.c.o
imgui/glfw/src/libglfw3.a: imgui/glfw/src/CMakeFiles/glfw.dir/vulkan.c.o
imgui/glfw/src/libglfw3.a: imgui/glfw/src/CMakeFiles/glfw.dir/window.c.o
imgui/glfw/src/libglfw3.a: imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_init.m.o
imgui/glfw/src/libglfw3.a: imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_joystick.m.o
imgui/glfw/src/libglfw3.a: imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_monitor.m.o
imgui/glfw/src/libglfw3.a: imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_window.m.o
imgui/glfw/src/libglfw3.a: imgui/glfw/src/CMakeFiles/glfw.dir/cocoa_time.c.o
imgui/glfw/src/libglfw3.a: imgui/glfw/src/CMakeFiles/glfw.dir/posix_thread.c.o
imgui/glfw/src/libglfw3.a: imgui/glfw/src/CMakeFiles/glfw.dir/nsgl_context.m.o
imgui/glfw/src/libglfw3.a: imgui/glfw/src/CMakeFiles/glfw.dir/egl_context.c.o
imgui/glfw/src/libglfw3.a: imgui/glfw/src/CMakeFiles/glfw.dir/osmesa_context.c.o
imgui/glfw/src/libglfw3.a: imgui/glfw/src/CMakeFiles/glfw.dir/build.make
imgui/glfw/src/libglfw3.a: imgui/glfw/src/CMakeFiles/glfw.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_16) "Linking C static library libglfw3.a"
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && $(CMAKE_COMMAND) -P CMakeFiles/glfw.dir/cmake_clean_target.cmake
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/glfw.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
imgui/glfw/src/CMakeFiles/glfw.dir/build: imgui/glfw/src/libglfw3.a
.PHONY : imgui/glfw/src/CMakeFiles/glfw.dir/build

imgui/glfw/src/CMakeFiles/glfw.dir/clean:
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src && $(CMAKE_COMMAND) -P CMakeFiles/glfw.dir/cmake_clean.cmake
.PHONY : imgui/glfw/src/CMakeFiles/glfw.dir/clean

imgui/glfw/src/CMakeFiles/glfw.dir/depend:
	cd /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations /Users/enricopuppo/Documents/GitHub/lib/cinolib/external/imgui/glfw/src /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src /Users/enricopuppo/Documents/GitHub/code/mesh_deform/Triangulations/build/imgui/glfw/src/CMakeFiles/glfw.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : imgui/glfw/src/CMakeFiles/glfw.dir/depend

