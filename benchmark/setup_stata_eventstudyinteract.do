version 18.0
cap which ftools
if _rc ssc install ftools, replace
cap which reghdfe
if _rc ssc install reghdfe, replace
cap which ivreg2
if _rc ssc install ivreg2, replace
cap which avar
if _rc ssc install avar, replace
cap which eventstudyinteract
if _rc net install eventstudyinteract, from("https://raw.githubusercontent.com/lsun20/EventStudyInteract/main/") replace
which eventstudyinteract
exit, clear