file(REMOVE_RECURSE
  "liblef.a"
  "liblef.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/lef.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
