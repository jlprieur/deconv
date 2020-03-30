###########################################################
# jlp_balloon
# From Effective TCL/TK Programming (Mark Harrison Michael McLennan)
# pages 257-260
#
# Version 01/03/2002
###########################################################
# Set it to 1 if user want help balloons
#           0 otherwise
set bh_info(active) 1

# Define "widget ddefault" options for Balloonhelp class: 
option add *Balloonhelp*background white widgetDefault
option add *Balloonhelp*foreground black widgetDefault

# If message larger than 12cm, wrap onto multiple left-justified lines:
option add *Balloonhelp.info.wrapLength 12c widgetDefault
option add *Balloonhelp.info.justify left widgetDefault
#option add *Balloonhelp.info.font \
#     -*-lucida-medium-r-normal-sans-*-120-* widgetDefault

# Create this frame:
toplevel .balloonhelp -class Balloonhelp \
   -background black -borderwidth 1 -relief flat

# Arrow bitmap...
#label .balloonhelp.arrow -anchor nw \
#   -bitmap @[file join $my_tcl_dir images arrow.xbm]
#pack .balloonhelp.arrow -side left -fill y

# Create label that will be further used for the help message:
label .balloonhelp.info 
pack .balloonhelp.info -side left -fill y

# Suppress window manager's decorative border:
wm overrideredirect .balloonhelp 1

# Forget it for now (but later will be called where and when needed)
wm withdraw .balloonhelp

###################################################
# BalloonHelpCreate 
###################################################
proc BalloonHelpCreate {win mesg} {
  global bh_info
# Load help message in a (global) array:
  set bh_info($win) $mesg

# When calling BalloonHelpPending,
# %W will be replaced by the name of the widget:
  bind $win <Enter> {BalloonHelpPending %W}
  bind $win <Leave> {BalloonHelpCancel}
}
###################################################
# BalloonHelpPending
# to pop up the balloon window after 1.5 seconds
###################################################
proc BalloonHelpPending {win} {
  global bh_info
# Cancel any other pending balloon:
  BalloonHelpCancel
# BalloonHelpShow will be called after 1.5 seconds later to pop up
# the balloon window:
  set bh_info(pending) [after 1500 [list BalloonHelpShow $win]]
# (the result is stored in the slot bh_info(pending), so that it can
# be used later in BalloonHelpCancel)
}
###################################################
# BalloonHelpCancel
# to cancel last pending balloon
###################################################
proc BalloonHelpCancel {} {
  global bh_info
  if {[info exists bh_info(pending)]} {
    after cancel $bh_info(pending)
    unset bh_info(pending)
  }
wm withdraw .balloonhelp
}
###################################################
# BalloonHelpShow
# to display the help message for a particular widget: 
###################################################
proc BalloonHelpShow {win} {
  global bh_info
  if {$bh_info(active)} {
    .balloonhelp.info configure -text $bh_info($win)
# rootx and rooty: coordinates of the upper-left corner of the widget:
    set x [expr [winfo rootx $win]+10]
    set y [expr [winfo rooty $win]+[winfo height $win]]
    wm geometry .balloonhelp +$x+$y
    wm deiconify .balloonhelp
# since the window is unmanaged, raise it to be able to see it:
    raise .balloonhelp
  }
unset bh_info(pending)
}
