# scRNA_data_modeling
Single cell data modeling
    
General workflow

Model for generating data: Wang et al. 2010. Biophysical Journal Volume 99 July 2010 29–39. https://doi.org/10.1016/j.bpj.2010.03.058

- Fixed points (True system)
- Vector field (True system)
- Solutions in the vector field (True system)
- Generate data over time 

Model for decoding data: Daniels, B. C., & Nemenman, I. (2015). Nature Communications, 6, 1?8. https://doi.org/10.1038/ncomms9133
Model for 2 species (genes): Sigmoidal model (6 parameters) + decaying (2 parameters) + interacions (4 parameters) --> 12 parameters

- Particle swarm set up: target function solutions_version.m
	+ Non-linear constraint related with the stability of the fixed points: nlc_version.m 
	+ Cost function: Linear combination of the different constraint functions
- Fixed points (Solutions system)
- Vector field (Solutions system)
- Solutions in the vector field (Solutions system)

Versions:

- WF_a: includes the model for generating data
- WF_b: without generated data
- WF_c: without generated data, nullclines in the constraint
- WF_d: without generated data, nullclines in the constraint, different sigmoid 