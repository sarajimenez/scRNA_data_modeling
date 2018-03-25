# scRNA_data_modeling
Single cell data modeling

List of the different versions

- Workflow_a.m
    * Model for generating data: Wang et al. 2010. Biophysical Journal Volume 99 July 2010 29–39. https://doi.org/10.1016/j.bpj.2010.03.058
    Fixed points (True system)
    Vector field (True system)
    Solutions in the vector field (True system)
    Generate data over time 
    * Model for decoding data: Daniels, B. C., & Nemenman, I. (2015). Nature Communications, 6, 1?8. https://doi.org/10.1038/ncomms9133
    Model for 2 species (genes): Sigmoidal model (2 parameters) + decaying (2 parameters) + interacions (4 parameters) --> 8 parameters
    Optimization set-up using particle swarm 
    
- Worflow_b.m
    * Model for decoding data: 
    Model for 2 species (genes): Sigmoidal model (6 parameters) + decaying (2 parameters) + interacions (4 parameters) --> 12 parameters

    
- Worflow_c.m
    * Model for decoding data: 
    Fixed points (Solutions system)
    Vector field (Solutions system)
    Solutions in the vector field (Solutions system)
    
- Worflow_d.m
    * particleswarm: 
    Eigenvalues constraint