file(REMOVE_RECURSE
  "libdef.a"
  "libdef.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/def.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
