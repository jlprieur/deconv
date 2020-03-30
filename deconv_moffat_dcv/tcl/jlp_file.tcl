#############################################################
# jlp_file.tcl
#
# JLP
# Version 06/03/2002
#############################################################
# This file is used as a package provider:
# package provide jlp_file 1.0
# But rather complex: see man pkg_mkIndex ...

#############################################################
# ReadFromFile
# To read the value of a keyword 
# Arguments:
# IN:
# param           Name of the parameter to read 
# OUT:
# value           Value of the parameter to read 
# status
#############################################################
proc dcv:ReadKeyFromFile {param value status} {

# To obtain address of value and status:
upvar $value my_value
upvar $status my_status

set my_status -1
set filename "dcv_param.dat"

set istat [file exists $filename]
  if {$istat != 1} {
#   puts "ReadKeyFromFile/Error: $filename does not exist"  
  return
  }

set istat [file readable $filename]
  if {$istat != 1} {
#   puts "ReadKeyFromFile/Error: $filename is not readable"  
  return
  }

# Removes trailing spaces:
set param [string trim $param]

set fp [open $filename r+]
  while {[gets $fp line] >= 0} {
     if [regexp "$param=" $line] {
        set eqsign_index [string first "=" $line] 
        if {$eqsign_index != -1} {
          set my_value [string range $line [expr $eqsign_index + 1] end]
          set my_status 0
          close $fp
          return
        }
     }
  }

# puts "ReadKeyFromFile/status = $my_status, error reading $param in $filename!"  
close $fp
return
}
#############################################################
# WriteKeyToFile
# To read the value of a keyword 
# Arguments:
# IN:
# param           Name of the keyword 
# value           Value of the keyword 
# OUT:
# status
#############################################################
proc dcv:WriteKeyToFile {param new_value status} {

# To obtain address of status:
upvar $status my_status

set my_status -1
set filename "dcv_param.dat"

  if { [file exists $filename] != 1} {
     puts "WriteKeyToFile/Error: $filename does not exist"  
     return
   }

  if { [file writable $filename] != 1} {
     puts "WriteKeyToFile/Error: $filename is not writable"  
     return
   }

set fp [open $filename r+]
  while {[gets $fp line] >= 0} {
     if [regexp $param $line] {
        set eqsign_index [string first "=" $line] 
        if {$eqsign_index != -1} {
          set old_value [string range $line [expr $eqsign_index + 1] end]
puts "Before: $line and new value: $new_value"
          regsub $old_value [string range $line [expr $eqsign_index + 1] end] \
                 $new_value [string range $line [expr $eqsign_index + 1] end] 
puts "After: $line"
          set my_status 0
          close $fp
          return
        }
     }
  }

puts "WriteKeyToFile/Error: Keyword not found in $filename"  
close $fp
return
}
#########################################################################
# dcv:FileSelect
# To select a file in current directory:
#########################################################################
proc dcv:FileSelect {data_filename} {
upvar $data_filename my_filename
  set tk_strictMotif 1
  set types {
    { {.F* Files} {.F*} }
    { {.f* Files} {.f*} }
    { {All Files} {*} }
  }
  set my_filename [tk_getOpenFile -title "Load data file" -filetypes $types]
# If "Cancel" get "" (empty string)
    puts "OK: filename = >$my_filename<"
} 
#######################################################################
# dcv:SaveToFile
#
#######################################################################
proc dcv:SaveToFile {w status} {
global dcv_met_items dcv_method 
global dcv_file_val dcv_interactive dcv_positive 
global dcv_alpha dcv_ss 

    upvar $status my_status
    set my_status 1
    dcv:Validate $w my_status
    if { $my_status != 0 } {return}

set filename "dcv_param.dat"

set fp [open $filename w]
puts $fp "##################################################################### "
puts $fp "# Parameters used by \"jlp_dcv\" "
puts $fp "# File automatically generated/updated by \"jlp_dcv.tcl\" "
puts $fp "# "
puts $fp "# JLP "
puts $fp "# Version 07/05/2002 "
puts $fp "##################################################################### "
puts $fp "# Input PSF: "
puts $fp "DCV_PSF=$dcv_file_val(psf)"
puts $fp "# Input image: "
puts $fp "DCV_INPUT=$dcv_file_val(in)"
puts $fp "# Output image: "
puts $fp "DCV_OUTPUT=$dcv_file_val(out)"
puts $fp "# Interactive: 1 if interactive, 0 if batch processing  "
puts $fp "INTERACTIVE=$dcv_interactive"
puts $fp "# Deconvolution method: tikho wiener maxen sqrtr ggauss mark" 
puts $fp "DCV_METHOD=$dcv_method"
puts $fp "# Alpha: parameter used for regularization" 
puts $fp "DCV_ALPHA=$dcv_alpha"
puts $fp "# ss: parameter used for sqrt method only" 
puts $fp "DCV_SS=$dcv_ss"

close $fp

set my_status 0
return
}
