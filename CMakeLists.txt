# CMakeList.txt : CMake project for SCF Solver, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("SCF Solver")

# Add source to this project's executable.
add_executable (CMakeTarget "main.cpp"  "src/Molecule/Molecule.cpp" "src/Molecule/Atom.cpp" "src/Basis/Basis.cpp" "src/Molecule/Atom.h" "src/Basis/Basis.h" "src/Molecule/Molecule.h" "src/Utils/Vector3.h" "src/Utils/Vector.cpp" "src/Basis/Shell.cpp" "src/Basis/Shell.h" "src/Scf/Scf.h" "src/Scf/Scf.cpp" "src/Parser/MoleculeParser.h" "src/Parser/BasisParser.h" "src/Parser/FileParser.h" "src/Parser/MoleculeParser.cpp")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET CMakeTarget PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
