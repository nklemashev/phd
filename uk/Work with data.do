// Merge data
cd "D:\\PhD\\Work with BCC2008 data\\FES Data Stata\\"
use data75, clear
forvalues i=76(1)99 {
	append using "data`i'"
}
save "data75_99"

// Work with BBC2008
cd "D:\\PhD\\Work with BCC2008 data\\"
use "BBC2008", clear
// No prices on insgr -> drop it
drop insgr
// Negative income -> drop it
drop if hhincome < 0
// All commodity names
local names bread cereals biscuits beef lamb pork bacon poul_oth fish butter oil_fats cheese eggs milkfres milkprod tea coffee softdrin sug_pres swe_choc potatoes oth_vegs fruit oth_food canteen oth_snac beer winespir cigs othtob rent mortgage ratesetc repair diy coal electric gas oil_oth furnit furnish elec_app othhheq consuma petcare postage telephon domservs fees_sub menouter womouter kidouter othcloth footwear chemgood p_servs motorveh maintmot pet_oil tax_ins railfare busfares othtrav aud_vis rec_toys book_new garden tvlicen entertai
// Generate consumption data
foreach var in `names' {
	gen x_`var' = `var' * totexp / p_`var'
}

// Lorenz curve data
sort year hhincome
// Total income, by years
by year: egen tIncome = total(hhincome)
// Total expenditure by years
by year: egen tExpend = total(totexp)
// Income shares
gen iShare = hhincome / tIncome
// Cumulative sum of income shares
by year: gen cuSumIShares = sum(iShare)
// Population ranks by income, by years
by year: egen hhRank = rank(hhincome)
// Max population rank, by years
by year: egen hhMaxRank = max(hhRank)
// Cumulative sum of population shares
gen cuSumPShares = hhRank / hhMaxRank
// Unique ranks
gen I = 1 if hhRank[_n] != hhRank[_n + 1]

// Lorenz curve figures
// Straight line
by year: egen hhMinRank = min(hhRank)
gen sLine = 0 if hhRank == hhMinRank & I
replace sLine = 1 if hhRank == hhMaxRank & I
// Lorenz
line cuSumIShares sLine cuSumPShares if year == 75 & I, legend(off)
// Lorenz auto
forvalues i = 75(1)99 {
	line cuSumIShares sLine cuSumPShares if year == `i' & I, legend(off) xtitle("Cumulative share of people from lowest to highest incomes") ytitle("Cumulative share of income earned")
	graph export "D:\\PhD\\Work with BCC2008 data\\Lorenz\\Lorenz`i'.eps", as(eps) preview(on) replace
}

// Prepare data for non-parametric analysis
// All data
local names bread cereals biscuits beef lamb pork bacon poul_oth fish butter oil_fats cheese eggs milkfres milkprod tea coffee softdrin sug_pres swe_choc potatoes oth_vegs fruit oth_food canteen oth_snac beer winespir cigs othtob rent mortgage ratesetc repair diy coal electric gas oil_oth furnit furnish elec_app othhheq consuma petcare postage telephon domservs fees_sub menouter womouter kidouter othcloth footwear chemgood p_servs motorveh maintmot pet_oil tax_ins railfare busfares othtrav aud_vis rec_toys book_new garden tvlicen entertai
foreach var in `names' {
	by year: egen xx_`var' = total(x_`var')
}
gen IU = 1 if year[_n] != year[_n + 1]
outsheet year xx_* using "D:\\PhD\\Work with BCC2008 data\\Ten parts by income\\Part0.csv" if IU == 1, comma nolabel noquote replace
drop IU xx_*

// In several parts by income
local names bread cereals biscuits beef lamb pork bacon poul_oth fish butter oil_fats cheese eggs milkfres milkprod tea coffee softdrin sug_pres swe_choc potatoes oth_vegs fruit oth_food canteen oth_snac beer winespir cigs othtob rent mortgage ratesetc repair diy coal electric gas oil_oth furnit furnish elec_app othhheq consuma petcare postage telephon domservs fees_sub menouter womouter kidouter othcloth footwear chemgood p_servs motorveh maintmot pet_oil tax_ins railfare busfares othtrav aud_vis rec_toys book_new garden tvlicen entertai
local parts 100
forvalues i = 1(1)`parts' {
	display "*******************************************"
	display "***********     i = `i'     ***************"
	display "*******************************************"
	local lb (`i'-1)/`parts'
	local ub `i' / `parts'
	gen I`i' = 1 if cuSumPShares > `lb' & cuSumPShares <= `ub'
	gen IU = 1 if I`i'[_n] != I`i'[_n + 1] & I`i'[_n] != .
	foreach var in `names' {
		by year: egen xx_`var' = total(x_`var') if I`i' == 1
	}
	outsheet year xx_* using "D:\\PhD\\Work with BCC2008 data\\Hundred parts by income\\Part`i'.csv" if IU == 1, comma nolabel noquote replace
	drop I`i' IU xx_*
}

// Extract data for Jini index
gen IU = 1 if cuSumPShares[_n] != cuSumPShares[_n + 1]
forvalues yr = 75(1)99 {
	outsheet cuSumPShares cuSumIShares using "D:\\PhD\\Work with BCC2008 data\\Jini\\Year`yr'.csv" if year == `yr' & IU == 1, comma nolabel noquote replace
}
drop IU

// IncomePerCapita
gen incPerCapita = .
forvalues i = 75(1)99 {
	qui sum year if year == `i'
	local num r(N)
	replace incPerCapita = `num' / tIncome if year == `i'
}

// Percentiles of income
bysort year: egen p10 = pctile(hhincome), p(10)
bysort year: egen p50 = pctile(hhincome), p(50)
bysort year: egen p90 = pctile(hhincome), p(90)
bysort year: egen p99 = pctile(hhincome), p(99)
gen p50p10 = p50 / p10
gen p90p50 = p90 / p50
gen p99p90 = p99 / p90

//  By centiles
gen IU = 0
local parts 100
forvalues i = 1(1)`parts' {
	display "*******************************************"
	display "***********     i = `i'     ***************"
	display "*******************************************"
	local lb (`i'-1)/`parts'
	local ub `i' / `parts'
	gen I`i' = 1 if cuSumPShares > `lb' & cuSumPShares <= `ub'
	replace IU = 1 if I`i'[_n] != I`i'[_n + 1] & I`i'[_n] != .
}

// Expenditures by centiles
gen tExpendByPerc = .
local parts 100
forvalues i = 1(1)`parts' {
	display "*******************************************"
	display "***********     i = `i'     ***************"
	display "*******************************************"
	by year: egen tmp = total(totexp) if I`i' == 1
    replace tExpendByPerc = tmp if I`i' == 1
	drop tmp
}
gen tExpendByPercShare = tExpendByPerc / tExpend

// Incomes by centiles
gen tIncomeByPerc = .
local parts 100
forvalues i = 1(1)`parts' {
	display "*******************************************"
	display "***********     i = `i'     ***************"
	display "*******************************************"
	by year: egen tmp = total(hhincome) if I`i' == 1
    replace tIncomeByPerc = tmp if I`i' == 1
	drop tmp
}
gen tIncomeByPercShare = tIncomeByPerc / tIncome
