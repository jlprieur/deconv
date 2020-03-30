#!/usr/bin/wish
#########################################################$
# Should be run with 
#     wish jlp_dcv.tcl        
#
##########################################################

##############################################################
# Load the script libraries.
#  Using packages would be good, but would require some more installation
#  overhead.  This will work without making things complex.
##############################################################

set exec_dir $env(EXEC)

# jlp_balloon should be first, since called by others: 
set my_dep {jlp_balloon.tcl jlp_tk_combo.tcl jlp_file.tcl jlp_edit.tcl}
foreach f $my_dep {
    eval source $f
puts "Evaluating $f"
}

###############################################################
proc dcv:CheckRadioButtons {w} {
global dcv_positive dcv_interactive 

# Default values:
   set dcv_interactive 1
   set dcv_positive 0
# From file if present:
   set value -1
   dcv:ReadKeyFromFile POSITIVE value status
   if { [lsearch -exact {0 1} $value] != -1} {set dcv_positive $value}
   set value -1
   dcv:ReadKeyFromFile INTERACTIVE value status
   if { [lsearch -exact {0 1} $value] != -1} {set dcv_interactive $value}
   set value -1
   dcv:ReadKeyFromFile TIME_INTERVAL value status
   if { [lsearch -exact {20 40} $value] != -1} {set dt $value}

frame $w.opt
frame $w.opt.left

# Radiobuttons:
label $w.options -justify left -text "Options:"
checkbutton $w.opt.left.b1 -text "Interactive"    \
            -variable dcv_interactive -relief flat
BalloonHelpCreate $w.opt.left.b1 \
      "Interactive mode (with display). Otherwise: batch only"
checkbutton $w.opt.left.b2 -text "Positive constraint" \
            -variable dcv_positive -relief flat
BalloonHelpCreate $w.opt.left.b2 \
      "When this flag is selected, \
       a positive constraint is applied \
       during minimization\n"

# Main frames:
pack $w.opt -side top -expand yes -pady 2 
pack $w.opt.left -side left -pady 2 -padx 2 \
     -expand yes -anchor w

# Fill left frame:
pack $w.options -side left -expand yes -in $w.opt.left \
        -pady 2 -padx 2 -fill x
pack $w.opt.left.b1 $w.opt.left.b2 -in $w.opt.left     \
        -side top -pady 1 -anchor w  -fill y

return
}

# showVars --
# Displays the values of one or more variables in a window, and
# updates the display whenever any of the variables changes.
#
# Arguments:
# w -           Name of new window to create for display.
# args -        Any number of names of variables.

proc showVars {w args} {
    catch {destroy $w}
    toplevel $w
    wm title $w "Variable values"
    label $w.title -text "Variable values:" -width 20 -anchor center \
            -font {Helvetica 18}
    pack $w.title -side top -fill x
    set len 1
    foreach i $args {
        if {[string length $i] > $len} {
            set len [string length $i]
        }
    }
    foreach i $args {
        frame $w.$i
        label $w.$i.name -text "$i: " -width [expr $len + 2] -anchor w
        label $w.$i.value -textvar $i -anchor w
        pack $w.$i.name -side left
        pack $w.$i.value -side left -expand 1 -fill x
        pack $w.$i -side top -anchor w -fill x
    }
    button $w.ok -text OK -command "destroy $w" -default active
    bind $w <Return> "tkButtonInvoke $w.ok"
    pack $w.ok -side bottom -pady 2
}
##########################################################
 proc dcv:Validate {w status} {
    global dcv_interactive dcv_positive 
    global dcv_file_list dcv_file_lbl dcv_file_val 
    global dcv_method dcv_met_items
    global dcv_alpha dcv_ss 

    upvar $status my_status 
    set my_status -1

    foreach i $dcv_file_list {
       set dcv_file_val($i) [string trim $dcv_file_val($i)] 
       if { [file isfile $dcv_file_val($i)] != 1 } {
        dcv:ErrorMsg $w.va1 "Validate/Error" \
             ">$dcv_file_val($i)< is not a file!"
        return
       }
    }
    set dcv_positive [string trim $dcv_positive] 
    if { [lsearch -exact {0 1} $dcv_positive] == -1} {
       dcv:ErrorMsg $w.va2 "Validate/Error" \
             "$dcv_positive not 0 nor 1!"
       return
    }

    set dcv_interactive [string trim $dcv_interactive] 
    if { [lsearch -exact {0 1} $dcv_interactive] == -1} {
       dcv:ErrorMsg $w.va2 "Validate/Error" "$dcv_interactive not 0 nor 1!"
       return
    }

set my_status 0

return
}
proc dcv:ErrorMsg {w title msg} {
    toplevel $w
    wm title $w $title 
    label $w.msg  -text $msg -justify left -anchor center \
            -font {Helvetica 12}
    pack $w.msg -side top -fill x

    button $w.ok -text OK -command "destroy $w" -default active
    bind $w <Return> "tkButtonInvoke $w.ok"
    pack $w.ok -side bottom -pady 2
}
proc dcv:exit_cmd {w} {
    set status 1
    dcv:SaveToFile $w status
    destroy $w
}


##############################################################
# Actual start:
##############################################################

wm withdraw .
set w .jlp_dcv
toplevel $w
wm title $w "jlp_dcv/configuration"
wm iconname $w "form"

# Create the tixLabelEntrys on the top of the dialog box
#
frame $w.top -border 1 -relief raised

# Edit "filenames" fields:
dcv:Edit_Files $w

# File selection (for future use...):
frame $w.file
set data_filename "" 
label $w.file.msg -text "Data file:" 
label $w.file.msg2 -textvariable data_filename -width 40 \
                   -relief sunken -bd 2 
button $w.file.select -text "Select File"  -underline 0 \
	-command "dcv:FileSelect data_filename"
# pack $w.file -side top -pady 2 
# pack $w.file.msg $w.file.msg2 $w.file.select -side left \
#           -padx 2 -pady 2 

# Method in ComboBox:
dcv:ComboBox_Method $w

# CheckButtons: 
dcv:CheckRadioButtons $w

# dcv parameters: (sig1, sig2, radmax)
dcv:Edit_dcv_Param $w

#########################################################
#
frame $w.rn
pack $w.rn -side top -pady 2 -expand yes 

label $w.rn.lbl -text "Deconvolution with dcv_deconv_2D: "  

set run_but [button $w.rn.run -text "Run"  -underline 0 \
	-command "dcv:Run $w"]

BalloonHelpCreate $w.rn.run \
" Run: Save config to file and run dcv_deconv_2D\
Stop: Stop dcv_deconv_2D"
#
pack $w.rn.lbl $w.rn.run -side left -fill x 

#########################################################
frame $w.t
set run_log [text $w.t.log -width 60 -height 10 \
        -borderwidth 2 -relief raised -setgrid true \
        -yscrollcommand {$w.t.scroll set}]
scrollbar $w.t.scroll -command {$w.t.log yview}
pack $w.t.scroll -side right -fill y
pack $w.t.log -side left -fill both -expand true
pack $w.t -side top -fill both -expand true

#########################################################
# Last line with Exit/Cancel/See Variables:
set status 1
frame $w.bt
pack $w.bt -side top -pady 2 -expand yes 

button $w.bt.exit -text "Exit" -underline 0 \
        -command "dcv:exit_cmd $w" -default active
BalloonHelpCreate $w.bt.exit \
      "Save to file (dcv_param.dat) and quit."
#
button $w.bt.save -text "Save"  -underline 0 \
	-command "dcv:SaveToFile $w status"
BalloonHelpCreate $w.bt.save \
      "Save to file (dcv_param.dat)."
#
button $w.bt.cancel -text "Cancel" -underline 0 \
        -command "destroy $w"
BalloonHelpCreate $w.bt.cancel \
      "Quit without saving to file (dcv_param.dat)."

pack $w.bt.exit $w.bt.save $w.bt.cancel -side left -fill x 

# Let's set some nice bindings for keyboard accelerators
#
    bind $w <Return> "tkButtonInvoke $w.bt.exit"
    bind $w <Alt-m> "focus $w.dcv_method"

pack $w.top -side top -fill both -expand yes
bind $w <Destroy> exit

########################################################################
# Start process
########################################################################
proc dcv:Run {w} {
    global dcv2D_process exec_dir 
    global run_log run_but
    set status 1
    dcv:SaveToFile $w status

# Bug here (why ?)
#    puts "EXEC=$env(EXEC)"

# The catch command guards against bogus commands.
# The variable dcv2D_process is set to an error message,
# or to the normal open return that is a file descriptor.
#
# The pipeline diverts error output from the command through the cat program.
# If you do not use cat like this, then the error output from the pipeline,
# if any, shows up as an error message when the pipeline is closed.
    if [catch {open "|$exec_dir/dcv_deconv_2D.exe |& cat"} dcv2D_process] {
         $run_log insert end "Error : $dcv2D_process \n"
    } else {
# JLP2002: I add this to avoid buffering problems which block the window:
#         fconfigure $dcv2D_process -buffering line
		 $run_log insert end "OK: $dcv2D_process\n"
# The run button is changed into a stop button after the program begins.
		 $run_but config -text Stop -command "dcv:Stop $w"
# Handle file events:
#		 fileevent $dcv2D_process readable "dcv:Log $w"
	    }
	     
	}

# Read and log output from the program
# The run:Log procedure is invoked whenever data can be read from the pipeline,
# and when end of file has been reached.
proc dcv:Log {w} {
	global dcv2D_process run_log
# Here EOF has been reached:
	if [eof $dcv2D_process] {
		dcv:Stop $w
# Here one line is read and inserted into the log:
	} else {
		gets $dcv2D_process line
		$run_log insert end $line\n
# The text widget's see operation is used to position the view on the
# text so the new line is visible to the user:
		$run_log see end
	}
}

# Stop the program and fix up the button
proc dcv:Stop {w} {
   global dcv2D_process run_but
   catch {close dcv2D_process}

# The button is restored to its run state 
   $run_but config -text "Run" -command "dcv:Run $w"

}
