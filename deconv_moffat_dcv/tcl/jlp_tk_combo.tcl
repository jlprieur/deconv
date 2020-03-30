#####################################################################
# Combo-box: 
# JLP
# Version 07/07/2002
#####################################################################
proc dcv:ComboBox_Method {w} {
global dcv_method 
global dcv_met_list dcv_met_items

set dcv_met_list {tikho wiener maxen sqrtr ggauss mark}
set dcv_met_items(tikho)   "Tikhonov"
set dcv_met_items(wiener)  "Wiener"
set dcv_met_items(maxen)   "Maximum entropy"
set dcv_met_items(sqrtr)   "sqrt(ss^2 + x^2)"
set dcv_met_items(ggauss)  "Generalized Gauss"
set dcv_met_items(mark)    "Markov"

# Set default values:
   set dcv_method "tikho" 
# If file is present:
   set value -1
   dcv:ReadKeyFromFile DCV_METHOD value status
   if { [lsearch -exact $dcv_met_list $value] != -1} {set dcv_method $value}

frame $w.m

# Menu to initialize dcv_method
  label $w.m.met_msg -width 16 -text "Method:"
  button $w.m.met -text "$dcv_met_items($dcv_method)" -width 20\
            -command "dcv:SelectMethod $w" 
  BalloonHelpCreate $w.m.met_msg \
      "Deconvolution method: currently six methods are available"


pack $w.m -side top -anchor w -pady 3 -padx 3
pack $w.m.met_msg $w.m.met -side left \
       -in $w.m -anchor w -pady 3 -padx 3

}
##################################
# Popup list to select method:
proc dcv:SelectMethod {w} {
global dcv_method dcv_met_list dcv_met_items

if { [winfo exists $w.m.met_popup] == 0} {
   listbox $w.m.met_popup -width 20 -height [llength $dcv_met_list]
     foreach i $dcv_met_list {
       $w.m.met_popup insert end $dcv_met_items($i) 
      }
    pack $w.m.met $w.m.met_popup -side top 
} else {
    pack $w.m.met_popup -after $w.m.met -side top 
}

bind $w.m.met_popup <Double-Button-1> {
       set value [$w.m.met_popup get active]
       pack forget $w.m.met_popup 
       foreach i $dcv_met_list {
         if { [string compare $dcv_met_items($i) $value] == 0 } {
            set dcv_method $i
            break
         }
       }
       $w.m.met configure -text $dcv_met_items($dcv_method)
       }  
}
