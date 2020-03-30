#############################################################
# jlp_edit.tcl
# Contains:
#
#  dcv:Edit_Files
#  dcv:Edit_dcv_Param
#
#############################################################
proc dcv:Edit_Files {w} {

#
# [Hint] The -options switch sets the options of the subwidgets.
# [Hint] We set the label.width subwidget option of the Controls to 
#        be 16 so that their labels appear to be aligned.
#
    global dcv_file_list dcv_file_lbl dcv_file_val 

    set dcv_file_list {psf in out}
    set dcv_file_lbl(psf)     "Input PSF: "
    set dcv_file_lbl(in)   "Input image: "
    set dcv_file_lbl(out)  "Output image: "

# Default values:
    set dcv_file_val(psf)    "s1_psf.fits"
    set dcv_file_val(in)     "s1_raw.fits"
    set dcv_file_val(out)    "s1_dcv.fits"

# Read in parameter file if present:
    dcv:ReadKeyFromFile DCV_PSF dcv_file_val(psf) status
    dcv:ReadKeyFromFile DCV_INPUT dcv_file_val(in) status
    dcv:ReadKeyFromFile DCV_OUTPUT dcv_file_val(out) status

# Title:
label $w.msg1 -text "Files"
frame $w.files
pack $w.msg1 $w.files -side top -anchor n -padx .5c -expand yes

   foreach i $dcv_file_list {
#    puts "jlp_edit.tcl/ $i : $dcv_file_val($i)"
    frame $w.files.$i
    label $w.files.$i.lbl -width 16 -text $dcv_file_lbl($i) 
    entry $w.files.$i.entry -width 40 -relief sunken -bd 2 \
          -textvariable dcv_file_val($i)
    pack $w.files.$i.lbl $w.files.$i.entry -side left -anchor w
    pack $w.files.$i -side top -anchor w -in $w.files
    }

}
proc dcv:Edit_dcv_Param {w} {
global dcv_alpha dcv_ss 

# Load default values:
  dcv:Default_dcv_Param 

# Read in parameter file if present:
  set status 1
  dcv:ReadKeyFromFile DCV_ALPHA dcv_alpha status
  dcv:ReadKeyFromFile DCV_SS dcv_ss status

# Create edit fields:
  frame $w.wp -relief ridge -borderwidth 4 
  pack $w.wp -side top -pady 8 -padx 2 -expand yes
  frame $w.wp.wp0
  label $w.wp.wp0.lab -text "Parameters:"
  button $w.wp.wp0.but1 -text "Default values"  \
        -command dcv:Default_dcv_Param
  button $w.wp.wp0.but2 -text "Advice"  \
        -command "dcv:ReadLogFile $w"
  pack $w.wp.wp0 -side top -pady 2 -expand yes -in $w.wp
  pack $w.wp.wp0.lab $w.wp.wp0.but1 $w.wp.wp0.but2 -side left -in $w.wp.wp0 \
       -anchor n -padx .5c -expand yes

  frame $w.wp.wp1
  pack $w.wp.wp1 -side top -anchor n -pady 2 -expand yes -in $w.wp

# alpha:
  label $w.wp.wp1.wp_1 -text "alpha:" 
  entry $w.wp.wp1.ed1 -width 5 -relief sunken -bd 2 \
          -textvariable dcv_alpha 
  BalloonHelpCreate $w.wp.wp1.wp_1 \
      "alpha: parameter used for regularization"

# ss:
  label $w.wp.wp1.wp_2 -text " ss:" 
  entry $w.wp.wp1.ed2 -width 5 -relief sunken -bd 2 \
          -textvariable dcv_ss 
  BalloonHelpCreate $w.wp.wp1.wp_2 \
      "ss: parameter used only by sqrt method" 

  foreach i [list 1 2] {
    pack $w.wp.wp1.wp_$i $w.wp.wp1.ed$i -side left -in $w.wp.wp1 \
            -anchor n -pady 2 -padx 2 -expand yes
  }

}
proc dcv:Default_dcv_Param {} {
  global dcv_alpha dcv_ss 
  set dcv_alpha 1.0 
  set dcv_ss 0.1 
}
###############################################################
# Read advice from "jlp_dcv.log"
###############################################################
proc dcv:ReadLogFile {w} {
set filename "jlp_dcv.log"

set istat [file exists $filename]
  if {$istat != 1} {
   dcv:ErrorMsg $w.rlf1 "ReadKeyFromFile/Error" \
           "Sorry: $filename does not exist. \n\
            This file is created by jlp_dcv, when processing data."
  return
  }

# Read content of logfile into dcv_log_text:
set fp [open $filename r]
set dcv_log_text [read $fp]
close $fp

# Create window with text, if not already there:
if { [winfo exists $w.rdlog] == 0} {
    toplevel $w.rdlog
    wm title $w.rdlog "Advice from \"$filename\" "

# Text:
    frame $w.rdlog.text
    text $w.rdlog.text.main -width 80 -height 20 \
          -yscrollcommand {$w.rdlog.text.ysbar set}
    scrollbar $w.rdlog.text.ysbar -orient vertical \
          -command {$w.rdlog.text.main yview}
    pack $w.rdlog.text -side top
    pack $w.rdlog.text.main -side left -expand yes \
          -fill both -in $w.rdlog.text 
    pack $w.rdlog.text.ysbar -side right -fill y \
          -in $w.rdlog.text 

# Commands:
    frame $w.rdlog.cmd
    button $w.rdlog.cmd.ok -text OK \
         -command "destroy $w.rdlog" -default active
    button $w.rdlog.cmd.update -text Update \
         -command "dcv:UpdateAdvice $w.rdlog.text.main" 
    pack $w.rdlog.cmd -side top -pady 2
    pack $w.rdlog.cmd.ok $w.rdlog.cmd.update -padx 4 -side left 

# Default command is "OK":
    bind $w.rdlog <Return> "tkButtonInvoke $w.ok"
} else {
# Erase old text:
    $w.rdlog.text.main configure -state normal 
    $w.rdlog.text.main delete 1.0 end 
}

    $w.rdlog.text.main insert end $dcv_log_text "typewriter" 
# Configure as "read-only" to prevent any modification by the user:
    $w.rdlog.text.main configure -state disabled

return
}
###############################################################
# Update advice from "jlp_dcv.log"
###############################################################
proc dcv:UpdateAdvice {wl} {
# Read content of logfile into dcv_log_text:
set filename "jlp_dcv.log"
set fp [open $filename r]
set dcv_log_text [read $fp]
close $fp

# Erase old text:
    $wl configure -state normal 
    $wl delete 1.0 end 
# Write new text:
    $wl insert end $dcv_log_text "typewriter"
    $wl configure -state disabled

return
}
