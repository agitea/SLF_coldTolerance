Repository to accompany

# Cold tolerance strategy and lower temperature thresholds of *Lycorma delicatula* egg masses ![SLF_nymph1-01](https://github.com/agitea/SLF_coldTolerance/assets/73284944/c8b8e2ee-0240-401d-8511-497af8635027)

*Icon of a 1-3 instar spotted lanternfly nymph designed by the Invasive Species Center*

### Anna J. Turbelin, Brent J. Sinclair, John Rost, Amanda D. Roe


*code written by Anna J Turbelin*

## SCRIPTS

1_SLF_coldTolerance_PATemp_F1S1 - script to plot temperature data from Blue Marsh lake where egg masses where collected and look at differences in autumn field acclimatization temperature between 2021 and 2022.

2_SLF_coldTolerance_disturbedEgg - script to look at the effect of removing individual eggs from egg masses on hatch.

3_SLF_coldTolerance_SCP_F2 - script to plot supercooling points of spotted lanternfly eggs.

4_SLF_coldTolerance_binomialDataset - script to create binomial dataset from percentage hatch values. Including binomial dataset with theoretical boundaries.

5_SLF_coldTolerance_tempTimeTreatmentHatch_F3 - script to summarise basic statistics on hatch for each treatments including Kruskal-Wallis rank sum test comparing haatch in individual treatments to hatch in control.

6_SLF_coldTolerance_GLMmodelSelection - script to assess which GLM to use to assess the effect of time and temperature on hatch and calculate LLT10, LLT50 and LLT90. Requires running script 4_SLF_coldTolerance_binomialDataset.

7_SLF_coldTolerance_clusteredEffect - script to assess the effect of clustered egg massess on hatch rate. Requires running script 4_SLF_coldTolerance_binomialDataset.

8_SLF_coldTolerance_LLT - script to calculate lower lethal temperatures (LLT10, LLT50, LLT90). Requires running script 4_SLF_coldTolerance_binomialDataset and adjusted_doseP_fun.

9_SLF_CFIA_obs_AHCCD_Temperature - script to calculate the reoccurence of cold events in cities where spotted lanternfly has been observed in Canada: i) daily minimum below -27.7 C, ii) 10 or more consecutive days with maximum temperature below -15 C, between 1973-2022 

adjusted_doseP_fun - function adapted from the dose.p() function from the MASS package (Venables and Ripley, 2002) to calculate LLTs.

## KEY FILES
*in data folder*

Dataset S1. Field and laboratory temperature exposure of Lycorma delicatula egg masses used in YR2 experiments. Climate data for Blue Marsh Lake Pennsylvania from NOAA (2023) online https://www.weather.gov/wrh

Dataset S2. Field temperature data from Blue Marsh Lake Pennsylvania where Lycorma delicatula egg masses were collected in YR1 (November 2021) and YR2 (December 2022). Climate data from NOAA (2023) online https://www.weather.gov/wrh

Dataset S3. Handling effects on survival. Results from testing the effect of removing eggs from egg masses and rearing individuals; includes hatch success of each treatment.

Dataset S4. Supercooling point (SCP) measurements. Results from exposing sets of eggs to decreasing temperatures until all eggs experience a SCP, eggs are returned to +6 °C to complete chill requirement and then eggs are warmed and checked for survival. 

Dataset S5.  Lower Lethal Limits assessment. Results from exposing groups of egg masses to set temperatures to identify lower lethal limits and hatch success assessment. 

Dataset S6.  Adjusted and Homogenized Canadian Climate Data (AHCCD) for stations near which L. delicatula was observed. Data from Environment and Climate Change Canada (Available online: https://climatedata.ca/download/#ahccd-download )

## OUTPUT
*in output folder*
