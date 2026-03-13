version 18.0
local input_path "$BENCH_INPUT"
local output_path "$BENCH_OUTPUT"

clear all
import delimited using `"`input_path'`", clear

scalar bench_nobs = _N

qui eventstudyinteract y g_* g0-g16, cohort(first_treat) control_cohort(never_treat) absorb(i.id i.t)
timer clear 1
timer on 1
qui eventstudyinteract y g_* g0-g16, cohort(first_treat) control_cohort(never_treat) absorb(i.id i.t)
timer off 1
timer list 1
scalar t1 = r(t1)

qui eventstudyinteract y g_* g0-g16, cohort(first_treat) control_cohort(never_treat) absorb(i.id i.id1 i.id2 i.t)
timer clear 2
timer on 2
qui eventstudyinteract y g_* g0-g16, cohort(first_treat) control_cohort(never_treat) absorb(i.id i.id1 i.id2 i.t)
timer off 2
timer list 2
scalar t2 = r(t2)

clear
input str28 engine str24 spec double seconds long nobs
"Stata eventstudyinteract" "id + t" . .
"Stata eventstudyinteract" "id + id1 + id2 + t" . .
end
replace seconds = t1 in 1
replace seconds = t2 in 2
replace nobs = `nobs' in 1/2
export delimited using `"`output_path'`", replace
exit, clear