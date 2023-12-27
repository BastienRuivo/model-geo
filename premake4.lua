solution "gKit2light"
    configurations { "debug", "release" }

    includedirs { ".", "src/gKit" }
    
    gkit_dir = path.getabsolute(".")
    
    configuration "debug"
        targetdir "bin/debug"
        defines { "DEBUG" }
        if _PREMAKE_VERSION >="5.0" then
            symbols "on"
            cppdialect "c++17"
        else
            buildoptions { "-std=c++17" }
            flags { "Symbols" }
        end
    
    configuration "release"
        targetdir "bin/release"
--~ 		defines { "NDEBUG" }
--~ 		defines { "GK_RELEASE" }
        if _PREMAKE_VERSION >="5.0" then
            optimize "full"
            cppdialect "c++17"
        else
            buildoptions { "-std=c++17" }
            flags { "OptimizeSpeed" }
        end
        
    configuration "linux"
        buildoptions { "-mtune=native -march=native" }
        buildoptions { "-W -Wall -Wextra -O3 -Wsign-compare -Wno-unused-parameter -Wno-unused-function -Wno-unused-variable", "-pipe" }
        links { "GLEW", 
            "SDL2", 
            "SDL2_image", 
            "GL", 
            "imgui",
            "boost_system",
            "boost_filesystem",
            "boost_program_options",
            "boost_regex",
            "boost_thread",
            "boost_chrono",
            "boost_date_time",
            "boost_locale",
            "boost_iostreams",
            "boost_serialization",
            "boost_log",
            "boost_log_setup",
            "boost_atomic",
            "boost_timer",
        }
    
    configuration { "linux", "debug" }
        buildoptions { "-g"}
        linkoptions { "-g"}
        buildoptions { "-fopenmp" }
        linkoptions { "-fopenmp" }
        buildoptions { "-flto"}
        linkoptions { "-flto"}
    
    configuration { "linux", "release" }
        buildoptions { "-fopenmp" }
        linkoptions { "-fopenmp" }
        buildoptions { "-flto"}
        linkoptions { "-flto"}
    
if _PREMAKE_VERSION >="5.0" then
    configuration { "windows", "codeblocks" }
        location "build"
        debugdir "."
        
        buildoptions { "-U__STRICT_ANSI__"} -- pour definir M_PI
        defines { "WIN32", "_WIN32" }
        includedirs { "extern/mingw/include" }
        libdirs { "extern/mingw/lib" }
        links { "mingw32", "SDL2main", "SDL2", "SDL2_image", "opengl32", "glew32" }
end
    
if _PREMAKE_VERSION >="5.0" then
    configuration { "windows" }
        defines { "WIN32", "_USE_MATH_DEFINES", "_CRT_SECURE_NO_WARNINGS" }
        defines { "NOMINMAX" } -- allow std::min() and std::max() in vc++ :(((

    configuration { "windows", "vs*" }
        location "build"
        debugdir "."
        
        system "Windows"
        architecture "x64"
        disablewarnings { "4244", "4305" }
        flags { "MultiProcessorCompile", "NoMinimalRebuild" }
        
        includedirs { "extern/visual/include" }
        libdirs { "extern/visual/lib" }
        links { "opengl32", "glew32", "SDL2", "SDL2main", "SDL2_image" }
end
    
    configuration "macosx"
        frameworks= "-F /Library/Frameworks/"
        buildoptions { "-std=c++17 -Wno-deprecated-declarations" }
        defines { "GK_MACOS" }
        buildoptions { frameworks }
        linkoptions { frameworks .. " -framework OpenGL -framework SDL2 -framework SDL2_image" }
    
 -- description des fichiers communs
gkit_files = { gkit_dir .. "/src/gKit/*.cpp", gkit_dir .. "/src/gKit/*.h" }
imgui_files = { gkit_dir .. "/src/imgui/*.cpp", gkit_dir .. "/src/imgui/*.h" }

project("tpImplicite")
    language "C++"
    kind "ConsoleApp"
    targetdir "bin"
    files ( imgui_files )
    files ( gkit_files )
    files { "projets/tpImplicite/*.cpp" }


