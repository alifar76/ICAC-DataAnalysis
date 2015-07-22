# ICAC-DataAnalysis

## src

A simple shell script with instructions to run analysis for 7 year outcome data to generate a multiply, rarefied OTU table.

Script assumes working within detachable hard drive as folder paths are hard-coded. 



## SensTraj

#### *Development of patient trajectories*
1. Create a **proportion of sensitivities** variable at each time point, which is a numeric value indicating the “degree” of sensitization. For example, if you were tested for 6 different allergens, and allergic to 3 (based on either IgE>=0.35 or a positive skin test), you have a value of 0.5 for that specific year. 
2. Use **latent class mixture models** that model the degree of sensitization over time. 
3. The modeling system places participants into groups based on their similar degree of sensitization across time.
4. Compare several iterations to determine the number of groups that provide the best fit for the model (3 groups in our case). 
 
#### *Notes*
* There are some participants who are missing data for the 4 groups at age 7 (michd_4grp_yr7 variable) but have data for the trajectories (sens_traj_y7). This is because the trajectories were developed on data that was imputed, so anyone missing sensitization data had this data imputed, which leads to more completeness for this variable. The michd_4grp_yr7 variable does not include imputed data. 
 
* Likewise, there are participants with missing trajectory data, but who have been included in the 4-group analysis – this is because the trajectories were developed only on those participants who were part of the main URECA high-risk cohort. When URECA began, 560 babies were enrolled in the main study (families who had history of asthma, allergies, etc), and 49 were enrolled as a small comparison group (families who had NO history of allergy, asthma, etc). 
It has been decided that this group of participants (n=22 of the 49 have been sequenced) should be removed from the analysis. This will maintain a consistent sample across both the microbiome analysis and URECA's main 7 year outcome analysis. There is a new variable in the CSV that indicates which cohort the child is in (variable = “cohort”) and it's suggested to only keep those who are indicated as “Main URECA cohort.”
