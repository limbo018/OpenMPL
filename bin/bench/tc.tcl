

set fileList [open "./a" r+]
puts "opened $fileList for read"

set newList [open "./out" w+]
puts "opened $newList for write"

while { [gets $fileList line] >=0 } {
  set fileName [lindex $line 0]
  layout create L1 $fileName
  set topCell [L1 topcell]
  puts $newList
  puts $fileName
  puts $topCell
  puts $newList "$fileName $topCell"
  layout delete L1
}

close $newList

exit
