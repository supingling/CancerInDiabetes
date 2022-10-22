
*********************
***R1 Diabetologia***
*********************

cd "/rfs/LRWE_Proj46/Shared/SL/analysis/Mortality"
set seed 8156248
use dmdb, clear
keep patid yearin agein out allcause cancerdeath gender ethnic imd bmi_grp smk_status
save dmdb0, replace

forvalues n = 1(1)500 {
	use dmdb0, clear
	distinct patid
	foreach var of varlist imd smk_status gender ethnic bmi_grp {
		replace `var' = 999 if `var' == .
	}

	preserve
	bsample, strata(imd smk_status gender ethnic bmi_grp)
	distinct patid
	save "db`n'", replace
	restore	
}

***Append bootstrap datasets***SL changed on 17/08 to be consistent with previous labels
cd "/rfs/LRWE_Proj46/Shared/SL/results/R2_Oct2022"
cls
clear
forvalues k = 1(1)500 {
	append using "DM_res`k'"
}
drop if bmi_grp != . & age != 72
drop lex_dur dur
tab out, m
replace bmi_grp = bmi_grp -1
foreach var of varlist gender ethnic imd smk_status bmi_grp {
	tab `var', m
	tab `var', m nolabel
}

reshape wide rate, i(boots age period gender ethnic imd smk_status bmi_grp ) j(out)
renames rate7 rate8 \ rall rcancer
gen ratio = 100*rcancer/rall
sum ratio
sort age gender ethnic imd smk_status bmi_grp period boots

foreach var of varlist rall rcancer ratio {
	by age gender ethnic imd smk_status bmi_grp period, sort : egen float mean_`var' = mean(`var')
}
by age gender ethnic imd smk_status bmi_grp period, sort : egen float p_lb = pctile(ratio), p(2.5)
by age gender ethnic imd smk_status bmi_grp period, sort : egen float p_ub = pctile(ratio), p(97.5)

/*
ratio rcancer rall if groups == 1								/*much smaller CIs*/
ratio rcancer rall, over(groups)
ratio rcancer rall if groups == 1, vce(bootstrap, reps(500))	/*bootstrap of bootstrap - what is ideal number for the two bootstraps?*/ 
*/
keep age period gender ethnic imd smk_status bmi_grp p_lb p_ub
duplicates drop
replace period = period - 0.5
save "prop_bootstrap", replace






