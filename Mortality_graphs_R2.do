/*
**#bootstrap proportion, run in HPC////
cd "/rfs/LRWE_Proj46/Shared/SL/results/R2_Oct2022"
cls
clear
forvalues k = 1(1)500 {
	append using "DM_res`k'"
}
drop if bmi_grp != . & age != 72
drop lex_dur
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

foreach var of varlist rall rcancer ratio {
	by age period gender ethnic imd smk_status bmi_grp, sort : egen float mean_`var' = mean(`var')
}
by age period gender ethnic imd smk_status bmi_grp, sort : egen float p_lb = pctile(ratio), p(2.5)
by age period gender ethnic imd smk_status bmi_grp, sort : egen float p_ub = pctile(ratio), p(97.5)

/*
ratio rcancer rall if groups == 1								/*much smaller CIs*/
ratio rcancer rall, over(groups)
ratio rcancer rall if groups == 1, vce(bootstrap, reps(500))	/*bootstrap of bootstrap - what is ideal number for the two bootstraps?*/ 
*/
keep age period gender ethnic imd smk_status bmi_grp p_lb p_ub
duplicates drop
replace period = period - 0.5
save "prop_bootstrap", replace
*/

global path "C:\Users\SupingLing\Desktop\UoL\KB\paper1_mortality"
cd "$path\submissions\Diabetologia\R2\Rerun_Oct2022"
**#clean and re-organise results database///
use DM_res_R2, clear
drop dur
gen out1 = "allcause" if out == 16
replace out1 = "cancerdeath" if out == 17
replace out1 = "colorectal"  if out == 18
replace out1 = "lung"  if out == 19
replace out1 = "liver"  if out == 20
replace out1 = "pancreas"  if out == 21
replace out1 = "gallbladder"  if out == 22
replace out1 = "breast"  if out == 23 
replace out1 = "endometrium"  if out == 24
replace out1 = "prostate" if out == 25
drop out

gen out =1 if out1 == "allcause" 
replace out = 2 if out1 == "cancerdeath" 
replace out =3  if out1 == "pancreas"  
replace out =4  if out1 == "liver"
replace out =5  if out1 == "gallbladder" 
replace out =6  if out1 == "endometrium" 
replace out =7  if out1 == "breast" 
replace out =8  if out1 == "prostate" 
replace out =9  if out1 == "lung"  
replace out =10  if out1 == "colorectal" 
sort out 
replace out1 = proper(out1)
sencode out1, replace
label define out1 1 "All-cause death" 2 "All cancer death", modify
drop out lex_dur
label define gender 1 "Men" 2 "Women"
label values gender gender
replace bmi_grp = bmi_grp -1
label define bmi 1 "18.5-25" 2 "25-30" 3 "30-35" 4 "≥35"
label values bmi_grp bmi
label define ethnic 1 "White" 2 "Non-White"
label values ethnic ethnic
label define imd 1 "Most affluent" 5 "Most deprived"
label values imd imd
label define smk 1 "Current smoker" 2 "Non-smoker" 3 "Ex-smoker"
label values smk_status smk
replace period = period - 0.5
save DM_res_for_graph, replace
**#export results for joinpoint regression///
preserve
foreach var of varlist rate lb ub {
	replace `var' = `var' * 100 if out1 ==5
}
gen se = (ln(ub) - ln(lb))/3.92
sort out1 age gender ethnic imd smk_status bmi_grp period
drop lb ub dur
export excel "$path\res.xls", firstrow(variable) replace
restore
**#For R2 response IMD all five quintiles: all-cause and all cancer mortality rates///
twoway 	///
(connected rate period  if imd==1, lcolor(blue%80) lwidth(vthin) mcolor(blue%80) msymbol(circle) msize(vsmall) sort)     /// 
(rcap lb ub period 		if imd==1, lcolor(blue%80) lwidth(vthin) mcolor(blue%80) msize(6-pt) sort )   ///
(connected rate period  if imd==2, lcolor(blue%50) lwidth(vthin) mcolor(blue%50) msymbol(diamond) msize(vsmall) sort)     ///
(rcap lb ub period 		if imd==2, lcolor(blue%50) lwidth(vthin) mcolor(blue%50) msize(6-pt)  sort )   ///
(connected rate period  if imd==3, lcolor(blue%20) lwidth(vthin) mcolor(blue%20) msymbol(triangle) msize(vsmall)sort)     ///
(rcap lb ub period 		if imd==3, lcolor(blue%20) lwidth(vthin) mcolor(blue%20) msize(6-pt)  sort)   ///
(connected rate period  if imd==4, lcolor(red%20) lwidth(vthin) mcolor(red%20) msymbol(square) msize(vsmall)sort)     ///
(rcap lb ub period 		if imd==4, lcolor(red%20) lwidth(vthin) mcolor(red%20) msize(6-pt)  sort )   ///
(connected rate period  if imd==5, lcolor(red%60) lwidth(vthin) mcolor(red%60) msymbol(circle_hollow) msize(vsmall)sort)     ///
(rcap lb ub period 		if imd==5, lcolor(red%60) lwidth(vthin) mcolor(red%60) msize(6-pt)  sort )   ///
if out1 < 3, by(out1, rows(1) legend(on) note("")) ///
legend(order(1 "The least deprived" 3 "2nd quintile" 5 "3rd quintile" 7 "4th quintile" 9 "The most deprived") ///
size(vsmall) color(black) margin(tiny) nobox region(fcolor(none) lcolor(none)) rows(1)) ///
ytitle("All cancer mortality rate (per 1000 person-years)", size(small) margin(zero)) ///
xtitle("Calendar year", size(small) margin(zero) alignment(middle)) xscale(outergap(2)) ///
xsize(11.75) ysize(8.25) by(, graphregion(margin(zero) fcolor(white) ifcolor(white)) ///
plotregion(margin(zero) fcolor(white) ifcolor(white)) ///
yrescale xrescale iyaxes ixaxes iytick ixtick iylabel ixlabel ixtitle iytitle) 	///
subtitle("") graphregion(margin(zero)) plotregion(margin(zero))	///
xlabel(1998(2)2018, labsize(small) angle(forty_five))  ///
ylabel(1 2 5 10 20 50, nogrid labsize(small) angle(hor)) yscale(log) name(rate_IMD2, replace)


twoway 	///
(connected rate period  if imd==1, lcolor(blue%80) lwidth(vthin) mcolor(blue%80) msymbol(circle) msize(vsmall) sort)     /// 
(connected rate period  if imd==2, lcolor(blue%50) lwidth(vthin) mcolor(blue%50) msymbol(diamond) msize(vsmall) sort)     ///
(connected rate period  if imd==3, lcolor(blue%20) lwidth(vthin) mcolor(blue%20) msymbol(triangle) msize(vsmall)sort)     ///
(connected rate period  if imd==4, lcolor(red%20) lwidth(vthin) mcolor(red%20) msymbol(square) msize(vsmall)sort)     ///
(connected rate period  if imd==5, lcolor(red%60) lwidth(vthin) mcolor(red%60) msymbol(circle_hollow) msize(vsmall)sort)     ///
if out1 < 3, by(out1, rows(1) legend(off) note("")) ///
ytitle("All cancer mortality rate (per 1000 person-years)", size(small) margin(zero)) ///
xtitle("Calendar year", size(small) margin(zero) alignment(middle)) xscale(outergap(2)) ///
xsize(11.75) ysize(8.25) by(, graphregion(margin(zero) fcolor(white) ifcolor(white)) ///
plotregion(margin(zero) fcolor(white) ifcolor(white)) ///
yrescale xrescale iyaxes ixaxes iytick ixtick iylabel ixlabel ixtitle iytitle) 	///
subtitle("") graphregion(margin(zero)) plotregion(margin(zero))	///
xlabel(1998(2)2018, labsize(small) angle(forty_five))  ///
ylabel(1 2 5 10 20 50, nogrid labsize(small) angle(hor)) yscale(log) name(rate_IMD1, replace)

graph combine rate_IMD1 rate_IMD2, xsize(8.25) ysize(8.25) cols(1) ///
imargin(zero) graphregion(margin(tiny) fcolor(white) ifcolor(white)) ///
plotregion(margin(tiny) fcolor(white) ifcolor(white)) name(rate_IMD, replace)
graph export "IMD_response.svg", as(svg) replace

**#For R2 response raw data: all-cause and all cancer mortality rates///

use general_t2dm, clear
recode age (1/5=1)(6/7=2)(8/9=3)(10/11=4), gen(age_grp)
label define age_grp 1 "35-60" 2  "60-69" 3 "70-79" 4 ">=80"
label values age_grp age_grp
collapse (sum) pop allcause cancerdeath colorectal lung liver pancreas gallbladder endometrium breast prostate, by(yr age_grp t2dm)
keep if t2dm ==1
local j = 1
foreach var of varlist allcause cancerdeath pancreas liver gallbladder endometrium breast prostate lung colorectal {
	gen rate`j' = `var'/pop * 1000
	local ++j
}
keep age_grp yr rate*
reshape long rate, i(age_grp yr) j(out1)
rename rate rate0
drop if yr == 2019
tempfile t2dm
save `t2dm', replace

/*R1 results, model included year * age only*/
use "C:\Users\SupingLing\Desktop\UoL\KB\paper1_mortality\submissions\Diabetologia\R1\rerun-Aug2022\DM_res_for_graph.dta", clear
keep if gender ==. & ethnic ==. & imd ==. & smk_status ==. & bmi_grp ==.
keep age period rate out1 
rename rate rate1
gen age_grp = 1 if age == 55
replace age_grp =2 if age == 65
replace age_grp =3 if age == 75
replace age_grp =4 if age == 85
rename period yr
drop if yr == 2019
merge 1:1 yr age_grp out1 using `t2dm', nogen
save `t2dm', replace

use "DM_res_for_graph.dta", clear
keep if gender ==. & ethnic ==. & imd ==. & smk_status ==. & bmi_grp ==.
keep age period rate out1 
rename rate rate2
gen age_grp = 1 if age == 55
replace age_grp =2 if age == 65
replace age_grp =3 if age == 75
replace age_grp =4 if age == 85
rename period yr
drop if yr == 2019
merge 1:1 yr age_grp out1 using `t2dm', nogen
keep if out1 < 3
replace rate0 = . if rate0 == 0

twoway ///
(connected rate2 yr if age_grp == 1, sort mcolor(blue) msymbol(circle_hollow) lcolor(blue) lpattern(longdash)) ///
(connected rate1 yr if age_grp == 1, sort mcolor(blue) msymbol(circle_hollow) lcolor(blue) lpattern(vshortdash)) ///
(connected rate0 yr  if age_grp == 1, sort mcolor(blue) msymbol(circle) lcolor(blue) lpattern(solid)) ///
(connected rate2 yr if age_grp == 2, sort mcolor(red) msymbol(diamond_hollow) lcolor(red) lpattern(longdash)) ///
(connected rate1 yr if age_grp == 2, sort mcolor(red) msymbol(diamond_hollow) lcolor(red) lpattern(vshortdash)) ///
(connected rate0 yr  if age_grp == 2, sort mcolor(red) msymbol(diamond) lcolor(red) lpattern(solid)) ///
(connected rate2 yr if age_grp == 3, sort mcolor(green) msymbol(square_hollow) lcolor(green) lpattern(longdash)) ///
(connected rate1 yr if age_grp == 3, sort mcolor(green) msymbol(square_hollow) lcolor(green) lpattern(vshortdash)) ///
(connected rate0 yr  if age_grp == 3, sort mcolor(green) msymbol(square) lcolor(green) lpattern(solid)) ///
(connected rate2 yr if age_grp == 4, sort mcolor(orange) msymbol(triangle_hollow) lcolor(orange) lpattern(longdash)) ///
(connected rate1 yr if age_grp == 4, sort mcolor(orange) msymbol(triangle_hollow) lcolor(orange) lpattern(vshortdash)) ///
(connected rate0 yr  if age_grp == 4, sort mcolor(orange) msymbol(triangle) lcolor(orange) lpattern(solid)) ///
, by(out1, note("") legend(on) yrescale xrescale iyaxes ixaxes iylabel ixlabel) ///
legend(order(1 "Modelled (age*year) - 55 years old" 2 "Modelled adjusting for duration - 55 years old" 3 "Raw data - 35-59 years old"  ///
			 4 "Modelled (age*year) - 65 years old" 5 "Modelled adjusting for duration - 65 years old" 6 "Raw data - 60-69 years old"  ///
			 7 "Modelled (age*year) - 75 years old" 8 "Modelled adjusting for duration - 75 years old" 9 "Raw data - 70-79 years old"  ///
			 10 "Modelled (age*year) - 85 years old" 11 "Modelled adjusting for duration - 85 years old"  12 "Raw data - 80-89 years old") ///
			 col(3) size(vsmall)) ///
xlabel(1998(2)2018, ang(45) grid) ///
xtitle("Calendar year") yscale(log) ///
ylabel(1 2 5 10 20 50 100, ang(h)) ///
ytitle("Mortality rate, per 1000 person-years")
graph export "rawdata_response.svg", as(svg) replace

**#/////proportion of all cancer out of all-cause death at each age
use DM_res_for_graph, clear
keep if out1<3
reshape wide rate lb ub, i(age period gender ethnic imd smk_status bmi_grp) j(out1)
gen prop = rate2 / rate1 * 100
merge 1:1 age period gender ethnic imd smk_status bmi_grp using prop_bootstrap, nogen keepusing(p_lb p_ub)

foreach var of varlist prop p_lb p_ub {
	replace `var' = . if p_ub>100
}
sort age gender ethnic imd smk_status bmi_grp period

**# proportion tables in the ESM
preserve
foreach var of varlist rate1-p_ub {
    tostring `var', gen (`var'1) format(%9.1f) force
}
gen allcause = rate11 + " (" + lb11 + ", " + ub11 + ")"
gen allcancer =  rate21 + " (" + lb21 + ", " + ub21 + ")"
gen proportion =  prop1 + " (" + p_lb1 + ", " + p_ub1 + ")"
keep age gender ethnic imd smk_status bmi_grp period allcancer allcause proportion
export excel age period all* proportion using "table.xlsx" if age!=72, sheet("TableS1") sheetreplace firstrow(variables)

local j = 2
foreach c in gender ethnic imd bmi_grp smk_status {
	export excel `c' period all* proportion using "table.xlsx" if `c'!=., sheet("TableS`j'") sheetreplace firstrow(variables)
	local ++j
}
restore

drop if imd>1 & imd <5
replace imd = 2 if imd ==5
label define imd 1 "Most affluent" 5 "Most deprived", modify
label values imd imd

foreach var of varlist gender ethnic imd smk_status bmi_grp {
	tab `var', m
	tab `var', m nolabel
}
label drop bmi smk imd gender ethnic

gen group = .
gen strata = ""
foreach var of varlist gender ethnic imd smk_status bmi_grp {
	replace group = `var' if `var'!=.
	replace strata = "`var'"  if `var'!=.
}
replace group = 1 if age == 55
replace group = 2 if age == 65
replace group = 3 if age == 75
replace group = 4 if age == 85
replace strata = "1age" if age!=72
replace strata = "6smk_status" if strata == "smk_status"
replace strata = "5bmi_grp" if strata == "bmi_grp"
replace strata = "2gender" if strata == "gender"
replace strata = "3ethnic" if strata == "ethnic"
replace strata = "4imd" if strata == "imd"
gen strata1 = "Age, years"					if strata == "1age"
replace strata1 = "Smoking status"  if strata == "6smk_status"
replace strata1 = "BMI, kg/m2"		if strata == "5bmi_grp" 
replace strata1 = "Gender"			if strata == "2gender" 
replace strata1 = "Ethnicity"		if strata == "3ethnic" 
replace strata1 = "Deprivation"		if strata == "4imd" 
sort strata
sencode strata1, replace
/*ethnicity: non-white don't plot*/
drop if ethnic == 2

tab strata, m
tab group, m

sort strata1 group period
bysort strata1: gen per = _n
replace per = per + 3 if per > 21
replace per = per + 3 if per > 45
replace per = per + 3 if per > 69
label define period ///
1  "1998" 25 "1998" 49 "1998" 73 "1998" ///
2  "1999" 26 "1999" 50 "1999" 74 "1999" ///
3  "2000" 27 "2000" 51 "2000" 75 "2000" ///
4  "2001" 28 "2001" 52 "2001" 76 "2001" ///
5  "2002" 29 "2002" 53 "2002" 77 "2002" ///
6  "2003" 30 "2003" 54 "2003" 78 "2003" ///
7  "2004" 31 "2004" 55 "2004" 79 "2004" ///
8  "2005" 32 "2005" 56 "2005" 80 "2005" ///
9  "2006" 33 "2006" 57 "2006" 81 "2006" ///
10 "2007" 34 "2007" 58 "2007" 82 "2007" ///
11 "2008" 35 "2008" 59 "2008" 83 "2008" ///
12 "2009" 36 "2009" 60 "2009" 84 "2009" ///
13 "2010" 37 "2010" 61 "2010" 85 "2010" ///
14 "2011" 38 "2011" 62 "2011" 86 "2011" ///
15 "2012" 39 "2012" 63 "2012" 87 "2012" ///
16 "2013" 40 "2013" 64 "2013" 88 "2013" ///
17 "2014" 41 "2014" 65 "2014" 89 "2014" ///
18 "2015" 42 "2015" 66 "2015" 90 "2015" ///
19 "2016" 43 "2016" 67 "2016" 91 "2016" ///
20 "2017" 44 "2017" 68 "2017" 92 "2017" ///
21 "2018" 45 "2018" 69 "2018" 93 "2018", modify 
label values per period

local j =20	
local margin = "medsmall"
local margin1 = "medsmall"
forvalues strata=1(1)6 {
	sum p_ub if strata1 == `strata' , d
	local u = round(`r(max)'+4.9,10)
twoway 	///
	(bar prop period, lcolor(blue%`j') lwidth(vthin) fcolor(blue%`j') sort)     ///
	(rcap p_lb p_ub period, lcolor(blue%60) lwidth(vthin) sort )   ///
	if strata1 == `strata' & group == 1, legend(off)  ///
	xlabel(1998(2)2018, labsize(small) tlength(1) angle(forty_five)) ///
	xtitle("Calendar year", size(small) margin(tiny)) ///
	graphregion(margin(`margin') fcolor(white) ifcolor(white)) ///
	plotregion(margin(`margin1') fcolor(white) ifcolor(white)) ///
	ytitle("Proportion (%)", size(small)) ///
	ylabel(0(10)`u', labsize(small) angle(hor) nogrid) name(prop`strata'1, replace) nodraw

twoway 	///
	(bar prop period, lcolor(red%`j') lwidth(vthin) fcolor(red%`j') sort)     ///
	(rcap p_lb p_ub period, lcolor(red%60) lwidth(vthin) sort )   ///
	if strata1 == `strata' & group == 2, legend(off)  ///
	xlabel(1998(2)2018, labsize(small) tlength(1) angle(forty_five)) ///
	xtitle("Calendar year", size(small) margin(tiny)) ///
	ytitle("Proportion (%)", size(small)) ///
	graphregion(margin(`margin') fcolor(white) ifcolor(white)) ///
	plotregion(margin(`margin1') fcolor(white) ifcolor(white)) ///
	ylabel(0(10)`u', labsize(small) angle(hor) nogrid) name(prop`strata'2, replace)	nodraw	

cap twoway 	///
	(bar prop period, lcolor(green%`j') lwidth(vthin) fcolor(green%`j') sort)     ///
	(rcap p_lb p_ub period, lcolor(green%60) lwidth(vthin) sort )   ///
	if strata1 == `strata' & group == 3, legend(off)  ///
	xlabel(1998(2)2018, labsize(small) tlength(1) angle(forty_five)) ///
	xtitle("Calendar year", size(small) margin(tiny)) ///
	ytitle("Proportion (%)", size(small)) ///
	graphregion(margin(`margin') fcolor(white) ifcolor(white)) ///
	plotregion(margin(`margin1') fcolor(white) ifcolor(white)) ///
	ylabel(0(10)`u', labsize(small) angle(hor) nogrid) name(prop`strata'3, replace)	nodraw			

cap twoway 	///
	(bar prop period, lcolor(orange%`j') lwidth(vthin) fcolor(orange%`j') sort)     ///
	(rcap p_lb p_ub period, lcolor(orange%60) lwidth(vthin) sort )   ///
	if strata1 == `strata' & group == 4, legend(off)  ///
	xlabel(1998(2)2018, labsize(small) tlength(1) angle(forty_five)) ///
	xtitle("Calendar year", size(small) margin(tiny)) ///
	ytitle("Proportion (%)", size(small)) ///
	graphregion(margin(`margin') fcolor(white) ifcolor(white)) ///
	plotregion(margin(`margin1') fcolor(white) ifcolor(white)) ///
	ylabel(0(10)`u', labsize(small) angle(hor) nogrid) name(prop`strata'4, replace)	nodraw				
}


graph combine ///
prop11 prop12 prop13 prop14 ///
prop21 prop22   ///
prop31 		   ///
prop41 prop42   ///
prop51 prop52 prop53 prop54 ///
prop61 prop62 prop63  ///
, rows(6) holes(7 11 12) imargin(zero) xsize(8.25) ysize(11.75) ///
graphregion(margin(medium) fcolor(white) ifcolor(white)) ///
plotregion(margin(tiny) fcolor(white) ifcolor(white)) ///
name(prop, replace)
graph export "F2.svg", as(svg) replace

/*for legend*/
local j = 20
twoway ///
(bar prop per if group ==1, lcolor(blue%`j') lwidth(vthin) fcolor(blue%`j') sort) ///
(rcap p_lb p_ub per if group == 1, msize(4-pt)  lcolor(blue%80) lwidth(thin) sort) ///
(bar prop per if group ==2, lcolor(red%`j') lwidth(vthin) fcolor(red%`j') sort) ///
(rcap p_lb p_ub per if group == 2, msize(4-pt) lcolor(red%80) lwidth(thin) sort) ///
(bar prop per if group ==3, lcolor(green%`j') lwidth(vthin) fcolor(green%`j') sort) ///
(rcap p_lb p_ub per if group == 3, msize(4-pt) lcolor(green%80) lwidth(thin) sort) ///
(bar prop per if group ==4, lcolor(orange%`j') lwidth(vthin) fcolor(orange%`j') sort) ///
(rcap p_lb p_ub per if group == 4, msize(4-pt) lcolor(orange%80) lwidth(thin) sort) ///
if strata1 ==1, ytitle("Proportion (%)") ytitle(, size(vsmall) margin(tiny)) ///
xlabel(1(2)87, labsize(vsmall) valuelabel nogrid ang(45)) ///
by(, legend(on size(*0.3) ring(0)) note("") graphregion(margin(zero) fcolor(white) ifcolor(white)) ///
plotregion(margin(zero) fcolor(white) ifcolor(white))) xsize(8.25) ysize(1.958) ///
subtitle("") xtitle("Calendar year", size(vsmall)) /// 
legend(title("Age (years)", size(vsmall) pos(9)) order(1 "55" 3 "65" 5 "75" 7 "85") ///
textfirst symplacement(center) rows(2) size(vsmall) symxsize(5) keygap(1))  ///
by(strata1, rows(6) yrescale iyaxes ixaxes iytick ixtick iylabel ixlabel) ///
ylabel(0(10)70, nogrid labsize(vsmall) angle(hor)) name(age, replace) 

twoway ///
(bar prop per if group ==1, lcolor(blue%`j') lwidth(vthin) fcolor(blue%`j') sort) ///
(rcap p_lb p_ub per if group == 1, msize(4-pt)  lcolor(blue%80) lwidth(thin) sort) ///
(bar prop per if group ==2, lcolor(red%`j') lwidth(vthin) fcolor(red%`j') sort) ///
(rcap p_lb p_ub per if group == 2, msize(4-pt) lcolor(red%80) lwidth(thin) sort) ///
(bar prop per if group ==3, lcolor(green%`j') lwidth(vthin) fcolor(green%`j') sort) ///
(rcap p_lb p_ub per if group == 3, msize(4-pt) lcolor(green%80) lwidth(thin) sort) ///
(bar prop per if group ==4, lcolor(orange%`j') lwidth(vthin) fcolor(orange%`j') sort) ///
(rcap p_lb p_ub per if group == 4, msize(4-pt) lcolor(orange%80) lwidth(thin) sort) ///
if strata1 ==1, ytitle("Proportion (%)") ytitle(, size(vsmall) margin(tiny)) ///
xlabel(1(2)87, labsize(vsmall) valuelabel nogrid ang(45)) ///
by(, legend(on size(*0.3) bplacement(northeast) ring(0)) note("") graphregion(margin(zero) fcolor(white) ifcolor(white)) ///
plotregion(margin(zero) fcolor(white) ifcolor(white))) xsize(8.25) ysize(1.958) ///
subtitle("") xtitle("Calendar year", size(vsmall)) /// 
legend(title("Gender", size(vsmall) pos(9)) order(1 "Men" 3 "Women") ///
notextfirst symplacement(center) rows(2) size(vsmall) symxsize(5) keygap(1))  ///
by(strata1, rows(6) yrescale iyaxes ixaxes iytick ixtick iylabel ixlabel) ///
ylabel(0(10)70, nogrid labsize(vsmall) angle(hor)) name(sex, replace) 

twoway ///
(bar prop per if group ==1, lcolor(blue%`j') lwidth(vthin) fcolor(blue%`j') sort) ///
(rcap p_lb p_ub per if group == 1, msize(4-pt)  lcolor(blue%80) lwidth(thin) sort) ///
(bar prop per if group ==2, lcolor(red%`j') lwidth(vthin) fcolor(red%`j') sort) ///
(rcap p_lb p_ub per if group == 2, msize(4-pt) lcolor(red%80) lwidth(thin) sort) ///
(bar prop per if group ==3, lcolor(green%`j') lwidth(vthin) fcolor(green%`j') sort) ///
(rcap p_lb p_ub per if group == 3, msize(4-pt) lcolor(green%80) lwidth(thin) sort) ///
(bar prop per if group ==4, lcolor(orange%`j') lwidth(vthin) fcolor(orange%`j') sort) ///
(rcap p_lb p_ub per if group == 4, msize(4-pt) lcolor(orange%80) lwidth(thin) sort) ///
if strata1 ==1, ytitle("Proportion (%)") ytitle(, size(vsmall) margin(tiny)) ///
xlabel(1(2)87, labsize(vsmall) valuelabel nogrid ang(45)) ///
by(, legend(on size(*0.3) bplacement(northeast) ring(0)) note("") graphregion(margin(zero) fcolor(white) ifcolor(white)) ///
plotregion(margin(zero) fcolor(white) ifcolor(white))) xsize(8.25) ysize(1.958) ///
subtitle("") xtitle("Calendar year", size(vsmall)) /// 
legend(title("Ethnicity", size(vsmall) pos(9)) order(1 "White" 3 "Non-White") ///
notextfirst symplacement(center) rows(2) size(vsmall) symxsize(5) keygap(1))  ///
by(strata1, rows(6) yrescale iyaxes ixaxes iytick ixtick iylabel ixlabel) ///
ylabel(0(10)70, nogrid labsize(vsmall) angle(hor)) name(ethnic, replace) 

twoway ///
(bar prop per if group ==1, lcolor(blue%`j') lwidth(vthin) fcolor(blue%`j') sort) ///
(rcap p_lb p_ub per if group == 1, msize(4-pt)  lcolor(blue%80) lwidth(thin) sort) ///
(bar prop per if group ==2, lcolor(red%`j') lwidth(vthin) fcolor(red%`j') sort) ///
(rcap p_lb p_ub per if group == 2, msize(4-pt) lcolor(red%80) lwidth(thin) sort) ///
(bar prop per if group ==3, lcolor(green%`j') lwidth(vthin) fcolor(green%`j') sort) ///
(rcap p_lb p_ub per if group == 3, msize(4-pt) lcolor(green%80) lwidth(thin) sort) ///
(bar prop per if group ==4, lcolor(orange%`j') lwidth(vthin) fcolor(orange%`j') sort) ///
(rcap p_lb p_ub per if group == 4, msize(4-pt) lcolor(orange%80) lwidth(thin) sort) ///
if strata1 ==1, ytitle("Proportion (%)") ytitle(, size(vsmall) margin(tiny)) ///
xlabel(1(2)87, labsize(vsmall) valuelabel nogrid ang(45)) ///
by(, legend(on size(*0.3) bplacement(northeast) ring(0)) note("") graphregion(margin(zero) fcolor(white) ifcolor(white)) ///
plotregion(margin(zero) fcolor(white) ifcolor(white))) xsize(8.25) ysize(1.958) ///
subtitle("") xtitle("Calendar year", size(vsmall)) /// 
legend(title("Deprivation", size(vsmall) pos(9)) order(1 "Least deprived" 3 "Most deprived") ///
notextfirst symplacement(center) rows(2) size(vsmall) symxsize(5) keygap(1))  ///
by(strata1, rows(6) yrescale iyaxes ixaxes iytick ixtick iylabel ixlabel) ///
ylabel(0(10)70, nogrid labsize(vsmall) angle(hor)) name(dep, replace) 

twoway ///
(bar prop per if group ==1, lcolor(blue%`j') lwidth(vthin) fcolor(blue%`j') sort) ///
(rcap p_lb p_ub per if group == 1, msize(4-pt)  lcolor(blue%80) lwidth(thin) sort) ///
(bar prop per if group ==2, lcolor(red%`j') lwidth(vthin) fcolor(red%`j') sort) ///
(rcap p_lb p_ub per if group == 2, msize(4-pt) lcolor(red%80) lwidth(thin) sort) ///
(bar prop per if group ==3, lcolor(green%`j') lwidth(vthin) fcolor(green%`j') sort) ///
(rcap p_lb p_ub per if group == 3, msize(4-pt) lcolor(green%80) lwidth(thin) sort) ///
(bar prop per if group ==4, lcolor(orange%`j') lwidth(vthin) fcolor(orange%`j') sort) ///
(rcap p_lb p_ub per if group == 4, msize(4-pt) lcolor(orange%80) lwidth(thin) sort) ///
if strata1 ==1, ytitle("Proportion (%)") ytitle(, size(vsmall) margin(tiny)) ///
xlabel(1(2)87, labsize(vsmall) valuelabel nogrid ang(45)) ///
by(, legend(on size(*0.3) bplacement(northeast) ring(0)) note("") graphregion(margin(zero) fcolor(white) ifcolor(white)) ///
plotregion(margin(zero) fcolor(white) ifcolor(white))) xsize(8.25) ysize(1.958) ///
subtitle("") xtitle("Calendar year", size(vsmall)) /// 
legend(title("Body Mass Index (Kg/m2)", size(vsmall) pos(9)) order(1 "18.5-25" 3 "25-30" 5 "30-35" 7 "≥35") ///
notextfirst symplacement(center) rows(2) size(vsmall) symxsize(5) keygap(1))  ///
by(strata1, rows(6) yrescale iyaxes ixaxes iytick ixtick iylabel ixlabel) ///
ylabel(0(10)70, nogrid labsize(vsmall) angle(hor)) name(bmi, replace) 

twoway ///
(bar prop per if group ==1, lcolor(blue%`j') lwidth(vthin) fcolor(blue%`j') sort) ///
(rcap p_lb p_ub per if group == 1, msize(4-pt)  lcolor(blue%80) lwidth(thin) sort) ///
(bar prop per if group ==2, lcolor(red%`j') lwidth(vthin) fcolor(red%`j') sort) ///
(rcap p_lb p_ub per if group == 2, msize(4-pt) lcolor(red%80) lwidth(thin) sort) ///
(bar prop per if group ==3, lcolor(green%`j') lwidth(vthin) fcolor(green%`j') sort) ///
(rcap p_lb p_ub per if group == 3, msize(4-pt) lcolor(green%80) lwidth(thin) sort) ///
(bar prop per if group ==4, lcolor(orange%`j') lwidth(vthin) fcolor(orange%`j') sort) ///
(rcap p_lb p_ub per if group == 4, msize(4-pt) lcolor(orange%80) lwidth(thin) sort) ///
if strata1 ==1, ytitle("Proportion (%)") ytitle(, size(vsmall) margin(tiny)) ///
xlabel(1(2)87, labsize(vsmall) valuelabel nogrid ang(45)) ///
by(, legend(on size(*0.3) bplacement(northeast) ring(0)) note("") graphregion(margin(zero) fcolor(white) ifcolor(white)) ///
plotregion(margin(zero) fcolor(white) ifcolor(white))) xsize(8.25) ysize(1.958) ///
subtitle("") xtitle("Calendar year", size(vsmall)) /// 
legend(title("Smoking status", size(vsmall) pos(9)) order(1 "Current smoker" 3 "Non-smoker" 5 "Ex-smoker") ///
notextfirst symplacement(center) rows(2) size(vsmall) symxsize(5) keygap(1))  ///
by(strata1, rows(6) yrescale iyaxes ixaxes iytick ixtick iylabel ixlabel) ///
ylabel(0(10)70, nogrid labsize(vsmall) angle(hor)) name(smk, replace) 	
foreach t in age sex ethnic dep bmi smk {
	graph export "`t'.svg", as(svg) replace name(`t')
}
restore

**#F1, F3, F4
use DM_res_for_graph, clear
drop if imd>1 & imd <5
replace imd = 2 if imd ==5
label define imd 1 "Most affluent" 5 "Most deprived", modify
label values imd imd

foreach var of varlist gender ethnic imd smk_status bmi_grp {
	tab `var', m
	tab `var', m nolabel
}
label drop bmi smk imd gender ethnic

gen group = .
gen strata = ""
foreach var of varlist gender ethnic imd smk_status bmi_grp {
	replace group = `var' if `var'!=.
	replace strata = "`var'"  if `var'!=.
}
replace group = 1 if age == 55
replace group = 2 if age == 65
replace group = 3 if age == 75
replace group = 4 if age == 85
replace strata = "1age" if age!=72
replace strata = "6smk_status" if strata == "smk_status"
replace strata = "5bmi_grp" if strata == "bmi_grp"
replace strata = "2gender" if strata == "gender"
replace strata = "3ethnic" if strata == "ethnic"
replace strata = "4imd" if strata == "imd"
gen strata1 = "Age, years"					if strata == "1age"
replace strata1 = "Smoking status"  if strata == "6smk_status"
replace strata1 = "BMI, kg/m2"		if strata == "5bmi_grp" 
replace strata1 = "Gender"			if strata == "2gender" 
replace strata1 = "Ethnicity"		if strata == "3ethnic" 
replace strata1 = "Deprivation"		if strata == "4imd" 
sort strata
sencode strata1, replace

tab strata, m
tab group, m

**#Figure 1: all-cause and all cancer mortality rates///
forvalues strata = 1(1)6 {
	twoway 	///
	(connected rate period  if group==1, lcolor(blue%60) lwidth(vthin) mcolor(blue%60) msymbol(circle) msize(vsmall) sort)     /// 
	(rcap lb ub period 		if group==1, lcolor(blue%60) lwidth(vthin) mcolor(blue%60) msize(6-pt) sort )   ///
	(connected rate period  if group==2, lcolor(red%60) lwidth(vthin) mcolor(red%60) msymbol(diamond) msize(vsmall) sort)     ///
	(rcap lb ub period 		if group==2, lcolor(red%60) lwidth(vthin) mcolor(red%60) msize(6-pt)  sort )   ///
	(connected rate period  if group==3, lcolor(green%60) lwidth(vthin) mcolor(green%60) msymbol(triangle) msize(vsmall)sort)     ///
	(rcap lb ub period 		if group==3, lcolor(green%60) lwidth(vthin) mcolor(green%60) msize(6-pt)  sort)   ///
	(connected rate period  if group==4, lcolor(orange) lwidth(vthin) mcolor(orange) msymbol(square) msize(vsmall)sort)     ///
	(rcap lb ub period 		if group==4, lcolor(orange) lwidth(vthin) mcolor(orange) msize(6-pt)  sort )   ///
	if out1==1 & strata1 ==`strata', by(out1 strata1, rows(1) legend(off) note("")) ///
	ytitle("All-cause mortality rate (per 1000 person-years)", size(small)) ytitle(, margin(zero)) ///
	xtitle("Calendar year", size (small)) xtitle(, margin(zero)) ///
	xsize(1.958) ysize(2.75) by(, graphregion(margin(zero) fcolor(white) ifcolor(white)) ///
	plotregion(margin(zero) fcolor(white) ifcolor(white)) ///
	yrescale xrescale iyaxes ixaxes iytick ixtick iylabel ixlabel) 	///
	subtitle("") graphregion(margin(zero)) plotregion(margin(zero))	///
	xlabel(1998(2)2018, labsize(small) angle(forty_five))  ///
	ylabel(3 5 10 20 50 100 150, nogrid labsize(small) angle(hor)) yscale(log) ///
	name(rate1`strata', replace) nodraw
}

forvalues strata = 1(1)6 {
	twoway 	///
	(connected rate period  if group==1, lcolor(blue%60) lwidth(vthin) mcolor(blue%60) msymbol(circle) msize(vsmall) sort)     /// 
	(rcap lb ub period 		if group==1, lcolor(blue%60) lwidth(vthin) mcolor(blue%60) msize(6-pt) sort )   ///
	(connected rate period  if group==2, lcolor(red%60) lwidth(vthin) mcolor(red%60) msymbol(diamond) msize(vsmall) sort)     ///
	(rcap lb ub period 		if group==2, lcolor(red%60) lwidth(vthin) mcolor(red%60) msize(6-pt)  sort )   ///
	(connected rate period  if group==3, lcolor(green%60) lwidth(vthin) mcolor(green%60) msymbol(triangle) msize(vsmall)sort)     ///
	(rcap lb ub period 		if group==3, lcolor(green%60) lwidth(vthin) mcolor(green%60) msize(6-pt)  sort)   ///
	(connected rate period  if group==4, lcolor(orange) lwidth(vthin) mcolor(orange) msymbol(square) msize(vsmall)sort)     ///
	(rcap lb ub period 		if group==4, lcolor(orange) lwidth(vthin) mcolor(orange) msize(6-pt)  sort )   ///
	if out1==2 & strata1 ==`strata', by(out1 strata1, rows(1) legend(off) note("")) ///
	ytitle("All cancer mortality rate (per 1000 person-years)", size(small) margin(zero)) ///
	xtitle("Calendar year", size (small) margin(zero)) ///
	xsize(1.958) ysize(2.75) by(, graphregion(margin(zero) fcolor(white) ifcolor(white)) ///
	plotregion(margin(zero) fcolor(white) ifcolor(white)) ///
	yrescale xrescale iyaxes ixaxes iytick ixtick iylabel ixlabel	) 	///
	subtitle("") graphregion(margin(zero)) plotregion(margin(zero))	///
	xlabel(1998(2)2018, labsize(small) angle(forty_five))  ///
	ylabel(1 2 5 10 20 50, nogrid labsize(small) angle(hor)) yscale(log) ///
	name(rate2`strata', replace) nodraw
}

graph combine ///
rate11 rate12 rate13 rate14 rate15 rate16 ///
rate21 rate22 rate23 rate24 rate25 rate26 ///
, rows(2) xsize(11.75) ysize(5.5) ///
imargin(zero) graphregion(margin(tiny) fcolor(white) ifcolor(white)) ///
plotregion(margin(tiny) fcolor(white) ifcolor(white)) name(rate, replace)
graph export "F1.svg", as(svg) replace
**#Figure 3 & 4: cancer-specific mortality rates, where possible
drop if out1<3
drop if (out1==5 | out1 == 6) & strata1>1
foreach o in rate lb ub{
	replace `o' = `o' * 100  /*per 100,000 person-years*/
}	
drop if lb < 0.01
///Diabetologia requires to show all x and y axis so I draw each graph individually///
local lint = 60
local mint = 50
local size = "vsmall"
local margin = "vsmall"
forvalues out = 3(1)10 {
	forvalues strata = 1(1)6 {
		sum lb if out1 == `out' & strata1 == `strata'
		
		if (`r(N)' != 0) {
			
			if (`r(min)' >= 0.01 & `r(min)' < 0.1) {
			local min = round(`r(min)'-0.005,0.01)
			}
			if (`r(min)' >= 0.1 & `r(min)' < 1) {
			local min = round(`r(min)'-0.05,0.1)
			}
			if (`r(min)' >= 1 & `r(min)' < 5) {
			local min = round(`r(min)'-0.5,1)
			}
			if (`r(min)' >= 5 & `r(min)' < 10) {
				local min = 5
			}
			if (`r(min)' >=10) {
				local min = round(`r(min)'-5,10)
			}
			di `min' `r(min)'
			sum ub if out1 == `out' & strata1 == `strata'
			local max = round(`r(max)'+50, 100)
			di `r(max)' `max'
			local title : label (out1) `out'
			niceloglabels `min' `max', local(ylab) style(125)
			
			twoway 	///
			(connected rate period  if group==1, lcolor(blue%`lint') lwidth(vthin) mcolor(blue%`mint') msymbol(circle) msize(vsmall) sort)     /// 
			(rcap lb ub period 		if group==1, lcolor(blue%`lint') lwidth(vthin) mcolor(blue%`mint') msize(6-pt) sort )   ///
			(connected rate period  if group==2, lcolor(red%`lint') lwidth(vthin) mcolor(red%`mint') msymbol(diamond) msize(vsmall) sort)     ///
			(rcap lb ub period 		if group==2, lcolor(red%`lint') lwidth(vthin) mcolor(red%`mint') msize(6-pt)  sort )   ///
			(connected rate period  if group==3, lcolor(green%`lint') lwidth(vthin) mcolor(green%`mint') msymbol(triangle) msize(vsmall)sort)     ///
			(rcap lb ub period 		if group==3, lcolor(green%`lint') lwidth(vthin) mcolor(green%`mint') msize(6-pt)  sort)   ///
			(connected rate period  if group==4, lcolor(orange%`lint') lwidth(vthin) mcolor(orange%`mint') msymbol(square) msize(vsmall)sort)     ///
			(rcap lb ub period 		if group==4, lcolor(orange%`lint') lwidth(vthin) mcolor(orange%`mint') msize(6-pt)  sort )   ///
			if out1 == `out' & strata1== `strata', ///
			graphregion(margin(`margin') fcolor(white) ifcolor(white) lcolor(white)) ///
			plotregion(margin(`margin') fcolor(white) ifcolor(white) lcolor(white))	///
			xlabel(1998(2)2018, labsize(`size') tlength(1) nogrid angle(forty_five)) ///
			ylabel(`min' `ylab' `max', nogrid labsize(`size') angle(hor))  ///
			ytitle("`title' cancer mortality rate", size(`size')) xtitle("Calendar year", size(`size')) ///
			xsize(1.958) ysize(2.065) yscale(log) legend(off)	///
			name(common`out'`strata', replace)
			}

	}
}
graph combine ///
common71 common73 common74 common75 common76 ///
common81 common83 common84 common85 common86 ///
common91 common92  common93 common94 common95 common96 ///
common101 common102  common93 common104 common105 common106 ///
, graphregion(margin(`margin') fcolor(white) ifcolor(white) lcolor(white)) ///
plotregion(margin(`margin') fcolor(white) ifcolor(white) lcolor(white))	///
holes(2 8) xsize(11.75) ysize(8.25) col(6) name(common, replace)
graph export "F3.svg", as(svg) replace name("common")

graph combine ///
common31 common32 common33 common34 common35 common36 ///
common41 common42 common43 common44 common45 common46 ///
, graphregion(margin(`margin') fcolor(white)) plotregion(margin(`margin') fcolor(white))	///
xsize(11.75) ysize(5.5) col(6) name(dmrelated, replace)
graph export "F4.svg", as(svg) replace name("dmrelated")
graph close _all

**#Table 2: joinpoint regression results
import excel "joint.xlsx", sheet("APCs") firstrow clear
destring APC95*, replace
foreach var of varlist APC APC95*  {
	tostring `var', force format(%4.1f) gen(`var'1)
}
tostring PValue, force format(%4.3f) gen(Pvalue)
replace Pvalue = "<0.001" if Pvalue == "0.000"
tostring SegmentStart, replace
tostring SegmentEnd, replace
gen period = SegmentStart + "-" + SegmentEnd
gen apc = APC1 + " (" + APC95LCL1 + ", " + APC95UCL1 + ")"
drop Model SegmentStart SegmentEnd APC APC1 APC95LCL APC95LCL1 APC95UCL APC95UCL1 APCSignificant TestStatistic PValue
order Pvalue, last
reshape wide  period apc Pvalue, i(age gender ethnic imd smk_status bmi_grp out1) j(Segment)
tostring age, gen(age1)
replace bmi_grp = "≥35.0" if bmi_grp == "?35"
replace bmi_grp = "18.5-24.9" if bmi_grp == "18.5-25"
replace bmi_grp = "25.0-29.9" if bmi_grp == "25-30"
replace bmi_grp = "30.0-34.9" if bmi_grp == "30-35"
replace ethnic = "Y-Others" if ethnic == "Non-White"
gen strata = ""
gen Group = age1 if age!=72
replace strata = "age" if age!=72
foreach var of varlist gender ethnic imd smk_status bmi_grp {
	replace Group = `var' if age == 72 & `var'!=""
	replace strata = "`var'" if age == 72 & `var'!=""
}
drop age* gender ethnic imd smk_status bmi_grp
sort out1 strata Group
order strata Group, after(out1)
tempfile apc
save `apc', replace

import excel "joint.xlsx", sheet("AAPCs") firstrow clear
foreach var of varlist AAPC*  {
	tostring `var', force format(%4.1f) gen(`var'1)
}
tostring PValue, force format(%4.3f) gen(Pvalue)
replace Pvalue = "<0.001" if Pvalue == "0.000"
gen aapc = AAPC1 + " (" + AAPCCILow1 + ", " + AAPCCIHigh1 + ")"
tostring age, gen(age1)
replace bmi_grp = "≥35.0" if bmi_grp == "?35"
replace bmi_grp = "18.5-24.9" if bmi_grp == "18.5-25"
replace bmi_grp = "25.0-29.9" if bmi_grp == "25-30"
replace bmi_grp = "30.0-34.9" if bmi_grp == "30-35"
replace ethnic = "Y-Others" if ethnic == "Non-White"
gen strata = ""
gen Group = age1 if age!=72
replace strata = "age" if age!=72
foreach var of varlist gender ethnic imd smk_status bmi_grp {
	replace Group = `var' if age == 72 & `var'!=""
	replace strata = "`var'" if age == 72 & `var'!=""
}
tab strata, m
tab Group, m
keep out1 strata Group aapc Pvalue
merge 1:1 out1 strata Group using `apc', nogen
gen strata1 = 1 if strata == "age"
replace strata1 = 2 if strata == "gender"
replace strata1 = 3 if strata == "ethnic"
replace strata1 = 4 if strata == "imd"
replace strata1 = 5 if strata == "bmi_grp"
replace strata1 = 6 if strata == "smk_status"
label define strata 1 "Age" 2 "Gender" 3 "Ethnicity" 4 "Deprivation" 5 "BMI group" 6 "Smoking status"
label values strata1 strata
sort out1 strata1 Group
order aapc, last
keep out1 *apc* period* Group
foreach var of varlist apc* period* {
	replace `var' = "–" if `var' == ""
}
export excel using "table.xlsx", sheet("table2") sheetreplace firstrow(variables)

**#graphical abstract
use "DM_res_for_graph.dta" , clear
keep if out1 ==2
keep if period == 1998 | period == 2008 | period == 2018
gen per = 1 if period == 1998
replace per = 2 if period == 2008
replace per = 3 if period == 2018
label define per 1 "1998" 2 "2008" 3 "2018"
label values per per
sort gender ethnic imd smk_status bmi_grp age
label define ethnic 2 "Others", modify
label define imd 1 "Least deprived", modify
drop if imd>1 & imd<5

local int = 50
local int2 = 20
twoway ///
(bar rate per if age == 55, lcolor(black%`int') lwidth(vthin)  fcolor(green%`int') sort ) ///
(rcap lb ub per if age == 55, lcolor(green) sort msize(4-pt)) ///
(bar rate per if age == 65, lcolor(black%`int') lwidth(vthin) fcolor(green%`int2') sort ) ///
(rcap lb ub per if age == 65, lcolor(green) sort msize(4-pt)) ///
(bar rate per if age == 75, lcolor(black%`int') lwidth(vthin) fcolor(red%`int2') sort ) ///
(rcap lb ub per if age == 75, lcolor(red) sort msize(4-pt)) ///
(bar rate per if age == 85, lcolor(black%`int') lwidth(vthin) fcolor(red%`int') sort) ///
(rcap lb ub per if age == 85, lcolor(red) sort msize(4-pt)) ///
if age !=72, ///
by(age, rows(1) note("") legend(off)) ///
xlabel(1(1)3, ang(45) valuelabel) xtitle("") ///
ylabel(1 2 5 10 20 50, ang(hor)) yscale(log) ///
xsize(5.875) ysize(2.75) name(age, replace)

twoway ///
(bar rate per if gender == 1, lcolor(black%`int') lwidth(vthin)  fcolor(red%`int2') sort ) ///
(rcap lb ub per if gender == 1, lcolor(red) sort msize(4-pt)) ///
(bar rate per if gender == 2, lcolor(black%`int') lwidth(vthin)  fcolor(red%`int') sort) ///
(rcap lb ub per if gender == 2, lcolor(red) sort msize(4-pt)) ///
if gender !=., ///
by(gender, rows(1) note("") legend(off)) ///
xlabel(1(1)3, ang(45) valuelabel) xtitle("") ///
ylabel(1 2 5 10 20 50, ang(hor)) yscale(log) ///
xsize(5.875) ysize(2.75) name(gender, replace)

twoway ///
(bar rate per if ethnic == 1, lcolor(black%`int') lwidth(vthin) fcolor(red%`int') sort ) ///
(rcap lb ub per if ethnic == 1, lcolor(red) sort msize(4-pt)) ///
(bar rate per if ethnic == 2, lcolor(black%`int') lwidth(vthin)  fcolor(green%`int') sort) ///
(rcap lb ub per if ethnic == 2, lcolor(green) sort msize(4-pt)) ///
if ethnic !=., ///
by(ethnic, rows(1) note("") legend(off)) ///
xlabel(1(1)3, ang(45) valuelabel) xtitle("") ///
ylabel(1 2 5 10 20 50, ang(hor)) yscale(log) ///
xsize(5.875) ysize(2.75) name(ethnic, replace)

twoway ///
(bar rate per if imd == 1, lcolor(black%`int') lwidth(vthin) fcolor(red%`int') sort ) ///
(rcap lb ub per if imd == 1, lcolor(red) sort msize(4-pt)) ///
(bar rate per if imd == 5, lcolor(black%`int') lwidth(vthin) fcolor(red%`int2') sort) ///
(rcap lb ub per if imd == 5, lcolor(red) sort msize(4-pt)) ///
if imd !=., ///
by(imd, rows(1) note("") legend(off)) ///
xlabel(1(1)3, ang(45) valuelabel) xtitle("") ///
ylabel(1 2 5 10 20 50, ang(hor)) yscale(log) ///
xsize(5.875) ysize(2.75) name(imd, replace)

twoway ///
(bar rate per if smk_status == 1, lcolor(black%`int') lwidth(vthin) fcolor(red%`int') sort ) ///
(rcap lb ub per if smk_status == 1, lcolor(red) sort msize(4-pt)) ///
(bar rate per if smk_status == 2, lcolor(black%`int') lwidth(vthin) fcolor(green%`int') sort) ///
(rcap lb ub per if smk_status == 2, lcolor(green) sort msize(4-pt)) ///
(bar rate per if smk_status == 3, lcolor(black%`int') lwidth(vthin) fcolor(red%`int2') sort) ///
(rcap lb ub per if smk_status == 3, lcolor(red) sort msize(4-pt)) ///
if smk_status !=., ///
by(smk_status, rows(1) note("") legend(off)) ///
xlabel(1(1)3, ang(45) valuelabel) xtitle("") ///
ylabel(1 2 5 10 20 50, ang(hor)) yscale(log) ///
xsize(5.875) ysize(2.75) name(smk_status, replace)

twoway ///
(bar rate per if bmi_grp == 1, lcolor(black%`int') lwidth(vthin) fcolor(red%`int2') sort ) ///
(rcap lb ub per if bmi_grp == 1, lcolor(red) sort msize(4-pt)) ///
(bar rate per if bmi_grp == 2, lcolor(black%`int') lwidth(vthin) fcolor(blue%`int2') sort ) ///
(rcap lb ub per if bmi_grp == 2, lcolor(blue) sort msize(4-pt)) ///
(bar rate per if bmi_grp ==3,lcolor(black%`int') lwidth(vthin)  fcolor(blue%`int2') sort ) ///
(rcap lb ub per if bmi_grp == 3, lcolor(blue) sort msize(4-pt)) ///
(bar rate per if bmi_grp == 4, lcolor(black%`int') lwidth(vthin)  fcolor(red%`int') sort) ///
(rcap lb ub per if bmi_grp == 4, lcolor(red) sort msize(4-pt)) ///
if bmi_grp !=., ///
by(bmi_grp, rows(1) note("") legend(off)) ///
xlabel(1(1)3, ang(45) valuelabel) xtitle("") ///
ylabel(1 2 5 10 20 50, ang(hor)) yscale(log) ///
xsize(5.875) ysize(2.75) name(bmi_grp, replace)
graph combine age gender ethnic imd bmi_grp smk_status , ///
row(3) xsize(11.75) ysize(8.25)
graph export "graphical_abstract2_oct2022.svg", as(svg) replace


import excel "C:\Users\SupingLing\Desktop\UoL\KB\paper1_mortality\submissions\Diabetologia\R2\Rerun_Oct2022\joint.xlsx", sheet("AAPCs") firstrow clear
keep if out1 == "All cancer death"
replace bmi_grp = ">35" if bmi_grp == "?35"
replace ethnic = "Others" if ethnic == "Non-White"
replace imd = "Least deprived" if imd == "Most affluent"
gen age1 = 1 if age == 55
replace age1 = 2 if age == 65
replace age1 = 3 if age == 75
replace age1 = 4 if age ==85
label define age1 1 "55" 2 "65" 3 "75" 4 "85"
label values age1 age1
drop if strlen(imd) == 1
sencode gender, replace
sencode ethnic, replace
sencode imd, replace
sencode bmi_grp, replace
sencode smk_status, replace

foreach var of varlist AAPC* {
	tostring `var', gen(`var'1) force format(%4.1f)
}
gen aapc1 = AAPC1 + " (" + AAPCCILow1 + ", " + AAPCCIHigh1 + ")"



local int = 50
twoway ///
(bar AAPC age1 , lcolor(none) fcolor(blue%`int') barwidth(0.8) sort) ///
(rcap AAPCCILow AAPCCIHigh age1, lcolor(blue) msize(4-pt) sort) ///
(scatter AAPC age1, sort mcolor(none) msymbol(point) mlabel(AAPC1) mlabcolor(black) mlabposition(12)) ///
if age !=72, ylabel(-6(2)8, ang(hor) nogrid) ///
yline(0, lwidth(medium) lpattern(solid) lcolor(black))  ///
ytitle("Annual change from 1998 to 2018 (%)", size(medsmall)) ///
xlabel(1(1)4, valuelabel) xtitle("") title("Age, years", color(black) size(medsmall)) ///
legend(off) ///
xsize(4.125) ysize(3.92) name(age, replace)

twoway ///
(bar AAPC gender, lcolor(none) fcolor(blue%`int') barwidth(0.8) sort) ///
(rcap AAPCCILow AAPCCIHigh gender, lcolor(blue) msize(4-pt) sort) ///
(scatter AAPC gender, sort mcolor(none) msymbol(point) mlabel(AAPC1) mlabcolor(black) mlabposition(12)) ///
if gender !=., ylabel(-6(2)8, ang(hor) nogrid) ///
yline(0, lwidth(medium) lpattern(solid) lcolor(black))  ///
ytitle("Annual change from 1998 to 2018 (%)", size(medsmall)) ///
xlabel(1(1)4, valuelabel) xtitle("") title("Gender", color(black) size(medsmall)) ///
legend(off) ///
xsize(4.125) ysize(3.92) name(gender, replace)

twoway ///
(bar AAPC ethnic, lcolor(none) fcolor(blue%`int') barwidth(0.8) sort) ///
(rcap AAPCCILow AAPCCIHigh ethnic, lcolor(blue) msize(4-pt) sort) ///
(scatter AAPC ethnic, sort mcolor(none) msymbol(point) mlabel(AAPC1) mlabcolor(black) mlabposition(12)) ///
if ethnic !=., ylabel(-6(2)8, ang(hor) nogrid) ///
yline(0, lwidth(medium) lpattern(solid) lcolor(black))  ///
ytitle("Annual change from 1998 to 2018 (%)", size(medsmall)) ///
xlabel(1(1)4, valuelabel) xtitle("") title("Ethnicity", color(black) size(medsmall)) ///
legend(off) ///
xsize(4.125) ysize(3.92) name(ethnic, replace)

twoway ///
(bar AAPC imd, lcolor(none) fcolor(blue%`int') barwidth(0.8) sort) ///
(rcap AAPCCILow AAPCCIHigh imd, lcolor(blue) msize(4-pt) sort) ///
(scatter AAPC imd, sort mcolor(none) msymbol(point) mlabel(AAPC1) mlabcolor(black) mlabposition(12)) ///
if imd !=., ylabel(-6(2)8, ang(hor) nogrid) ///
yline(0, lwidth(medium) lpattern(solid) lcolor(black))  ///
ytitle("Annual change from 1998 to 2018 (%)", size(medsmall)) ///
xlabel(1(1)4, valuelabel) xtitle("") title("Deprivation", color(black) size(medsmall)) ///
legend(off) ///
xsize(4.125) ysize(3.92) name(imd, replace)

twoway ///
(bar AAPC smk_status, lcolor(none) fcolor(blue%`int') barwidth(0.8) sort) ///
(rcap AAPCCILow AAPCCIHigh smk_status, lcolor(blue) msize(4-pt) sort) ///
(scatter AAPC smk_status, sort mcolor(none) msymbol(point) mlabel(AAPC1) mlabcolor(black) mlabposition(12)) ///
if smk_status !=., ylabel(-6(2)8, ang(hor) nogrid) ///
yline(0, lwidth(medium) lpattern(solid) lcolor(black))  ///
ytitle("Annual change from 1998 to 2018 (%)", size(medsmall)) ///
xlabel(1(1)4, valuelabel) xtitle("") title("Smoking status", color(black) size(medsmall)) ///
legend(off) ///
xsize(4.125) ysize(3.92) name(smk_status, replace)

twoway ///
(bar AAPC bmi_grp, lcolor(none) fcolor(blue%`int') barwidth(0.8) sort) ///
(rcap AAPCCILow AAPCCIHigh bmi_grp, lcolor(blue) msize(4-pt) sort) ///
(scatter AAPC bmi_grp, sort mcolor(none) msymbol(point) mlabel(AAPC1) mlabcolor(black) mlabposition(12)) ///
if bmi_grp !=., ylabel(-6(2)8, ang(hor) nogrid) ///
yline(0, lwidth(medium) lpattern(solid) lcolor(black))  ///
ytitle("Annual change from 1998 to 2018 (%)", size(medsmall)) ///
xlabel(1(1)4, valuelabel) xtitle("") title("Body mass index (Kg/m2)", color(black) size(medsmall)) ///
legend(off) ///
xsize(4.125) ysize(3.92) name(bmi_grp, replace)
graph combine age gender ethnic imd smk_status bmi_grp, ///
row(3) xsize(8.25) ysize(11.75)
graph export "graphical_abstract_oct2022.svg", as(svg) replace



