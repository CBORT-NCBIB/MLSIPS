Maximum-Likelihood Single Input Polarization Sensitive Processing
=======

Repository Description
-----------
This code repository works to extract polarimetric properties from OFDI data acquired with a single input state. Systems **must** have polarization diverse detection and sufficient system PMD. This code also includes a mechanism for determining the PMD / spectral polarization spread of systems given a sufficient amount data, which can be used to determine the system-specific error on SIPS measurements.

Introduction Script Descriptions
-----------
### 1. StarterScript.m
Presents how to analyze tomogram data with SIPS processing given a system compensation structure and a
tomogram acquired with polarization diverse detection

_Required Inputs_
1. **System Compensation**: Your own, computed via SystemCompensationStarterScript, or an example 
2. **Spectrally binned tomogram**: Your own, computed with recstrTom, or an example

### 2. SystemCompensationStarterScript.m
Presents how to compute a system compensation structure given and a large sampling of Stokes vectors or tomograms
acquired with polarization diverse detection. Shows how to assess to overall spectral polarization spread of the system
and estimates overall reliability

_Required Inputs_
1. **Spectrally-binned Stokes Vectors**: Your own, computed via recstrTom or your prefered reconstruction method
