# S1-Doppler-Anomaly
In this repository, the main routine "read_dop_nc.pro" is able to correct the Doppler Anomaly product from S1 OCN file from instrumental noise and lack of attitude information (see DANILO, Céline et BELLA, Gábor. Extracting Usable Geophysical Doppler Properties from Sentinel–1 for Coastal Monitoring. In : 2018 Doppler Oceanography from Space (DOfS). IEEE, 2018. p. 1-5.). This solution is applicable in the coastal region. The precision is supposed to be in the order of 10 Hz which corresponds to 0.5 m/s.

These routine are written in IDL. I have tested them with IDL version 7.1.1. 
All these routines need to be put in the same repository.

The data repository path has to be specified at line 22 of the "read_dop_nc.pro" file.
