# Heterogeneity-unpredictability-parasite-outbreaks
Data and code associated with the manuscript titled "Host heterogeneity drives unpredictability in parasite outbreaks".

Details of the experimental setups for this work can be found in the paper titled “Host heterogeneity drives unpredictability in parasite outbreaks”.

Files and variables

File: Transmission_experiment.R
Description: Code to analyse transmission experiment data.

File: Susceptibility_infectiousness.R
Description: Code to analyse susceptibility and infectiousness data.

File: Trib_model_3.0.nlogo
Description: Agent based model from the paper.

File: SI_analyses.R
Description: Code to conduct the analyses presented in the SI.

File: transmission_experiment.csv
Description: Data from the transmission experiment. 
Variables
•	treatment: experimental treatment
•	replicate_no: experimental population within the treatment (populations 1-5)
•	dissection_no: Dissection week (weeks 1-8)
•	id: individual larval ID
•	trophont: number of trophonts
•	gamont: number of gamonts
•	gut_gametocyst: number of gametocysts
•	par_presence: parasite presence/absence (0: absent; 1: present)

File: susceptibility.csv
Description: Data from the initial susceptibility assay.
Variables
•	colony: colony identity
•	group: experimental or control
•	day_removed: experimental day larva was removed from flour
•	day_dissected: experimental day larva was dissected
•	id: individual ID
•	trophont: number of trophonts
•	gamont: number of gamonts
•	gut_gametocyst: number of gametocysts
•	par_presence: parasite presence/absence (0: absent; 1: present)
•	experimental_repeat: experimental repeat number

File: infectiousness.csv
Description: Data from the initial infectiousness assay.
Variables
•	colony: colony identity 
•	day_removed: experimental day larva was removed from flour
•	id: individual ID
•	oocysts: number of oocysts on haemocytometer slide
•	par_presence: parasite presence/absence (0: absent; 1: present)
•	experimental_repeat: experimental repeat number

File: gametocyst_oocyst.csv
Description: Data to test for covariation between gametocysts and oocysts.
Variables
•	colony: colony identity
•	day_removed: experimental day larva was removed from flour
•	id: individual ID
•	gametocysts: number of gametocysts in frass
•	oocysts: number of oocysts on haemocytometer slide
•	par_presence: parasite presence/absence (0: absent; 1: present)
•	experimental_repeat: experimental repeat number

File: development.csv
Description: Data to analyse development rates.
Variables
•	colony: colony identity
•	group: exposed (experimental) vs unexposed (control)
•	day: experimental day
•	development: 0 if larva had not pupated, 1 if it had pupated
•	experimental_repeat: experimental repeat number
