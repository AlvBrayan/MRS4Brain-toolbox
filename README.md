![image](https://github.com/AlvBrayan/MRS4Brain-toolbox/assets/137405058/426b8e15-a07f-4758-bc38-48535edbc710)

## Description:
The *MRS4Brain Toolbox* is developed by [*MRS4Brain*](https://www.epfl.ch/labs/mrs4brain/) research group @ CIBM MRI EPFL AIT and was designed to offer advanced functionalities for Bruker preclinical MRSI data, encompassing preprocessing, fitting, quantification, semi-automatic quality control, co-registration and segmentation of metabolic maps using anatomical images, all conveniently integrated within a single open-source graphical user interface (GUI). The development of this user-friendly toolbox aims to streamline the processing workflow and enhance the accessibility of MRSI for researchers in the preclinical field.

The *MRS4Brain Toolbox* is written in MATLAB 2023 or later (MathWorks, USA) and encompasses three distinct spectroscopy modalities (MRSI, single voxel MRS, dMRS), with our primary focus directed towards MRSI.

### Supported file formats:
-	raw Bruker Paravision 360 (version 1.1 (*fid*), 3.X (*fid_proc.64*)) data

### Supported methods:
-	1H-FID-MRSI, 1H-PRESS-MRSI
-	1H-SVS
-	1H-dMRS

## Functionalities:

-	The main functionalities are described in the document [MRS4Brain_Toolbox.pdf](https://github.com/AlvBrayan/MRS4Brain-toolbox/blob/main/MRS4Brain_Toolbox.pdf) , published in NMR in Biomedicine (Simicic D, Alves B, Mosso J, Briand G, Lê TP, van Heeswijk RB, Starčuková J, Lanz B, Klauser A, Strasser B, Bogner W, Cudalbu C. *Fast High-Resolution Metabolite Mapping in the rat Brain Using 1H-FID-MRSI at 14.1 T*, [DOI: 10.1002/nbm.5304](https://doi.org/10.1002/nbm.5304)) and in the presented abstract at ISMRM 2024 (Briand G., Alves B., Mosso J, Pierzchala K., Near J., Lanz B., Cudalbu C. *MRS4Brain Toolbox: an harmonized and accessible workflow for preclinical MRSI data processing*, ISMRM 2024 Singapore)
-	External software:
 
    -	The *MRS4Brain Toolbox* utilizes an additional software, ANTs, for registration and segmentation purposes. ANTs is designed for the UNIX environment. To facilitate this on non-UNIX systems, we employ another tool called Docker Desktop, enabling the containerization of Linux-built software. Details can be found here: [Registration_MRS4Brain_Toolbox.pdf](https://github.com/AlvBrayan/MRS4Brain-toolbox/blob/main/Registration_MRS4Brain_Toolbox.pdf)

    -	ANTs is combined with an in-house developed template (Wistar adult rat brain, acquisitions at 14.1T) based on the SIGMA atlas, customized to handle the dimensions of the rodent brain. 

    -	LCModel (version 6.3) is used for metabolite fitting. LCModel is now freely available and can be downloaded from here: http://s-provencher.com/lcmodel.shtml or [schorschinho/LCModel: A collection of compiled LCModel binaries for various operating systems. (github.com)](https://github.com/schorschinho/LCModel)

    -	The preprocessing steps in SVS are based on FID-A functions [CIC-methods/FID-A: Toolbox for simulation and processing of in-vivo magnetic resonance spectroscopy (MRS) data (github.com)](https://github.com/CIC-methods/FID-A)

-	Live Demos on how to use the toolbox can be found here: https://www.epfl.ch/labs/mrs4brain/links/live-demos/

-	Data sets are available here:

    - RAW data
    - Same RAW data + processing
    - Control Files for LCModel
    - Metabolite basis set (for AD=1.3ms) 

## Getting started:

### Prerequisites:

  - MATLAB version 2023a or later version (on all operating systems)
  - ANTs (on all operating systems)
  - LCModel version 6.3 (on all operating systems)
  - Docker (for Windows & MAC OS,  https://www.docker.com/products/docker-desktop/)
                
### Installation: 

- **On Windows & MAC OS**:
  1.	Download and install docker desktop (no account is required)
  2.	Start docker desktop
  3.	Run in Windows Powershell: `wsl --update` (if `wsl` is not installed, please run `wsl --install`)
  4.	Run in Windows Powershell: `docker pull antsx/ants:latest`
  5.	Check that the image is successfully installed on docker by running in Windows Powershell : `docker run --rm antsx/ants /opt/ants/bin/antsRegistration`. If the console returns something akin to what is seen in the printscreen below, ANTs has been successfully installed on your computer docker !

![image](https://github.com/AlvBrayan/MRS4Brain-toolbox/assets/137405058/c0046519-ebd1-4a5e-a031-0522b8f337ca)

- **On Linux** : No need for docker desktop. Just launch the *MRS4BRAIN_toolbox.mlapp*

- **For all operating systems**:
  1.	Launch your version of MATLAB
  2.	Add to path the folder containing the MRS4Brain_toolbox (Selected folders and subfolders)
  3.	Launch the MRS4BRAIN_toolbox command (either via your MATLAB console or with any commands you would use in your operating system)
  4.	After selecting the modality of your choice, click on parameters button (symbolized with the symbol below) and select the necessary folders for your task (LCModel location, Template masks, Config files, etc.)
 
![image](https://github.com/AlvBrayan/MRS4Brain-toolbox/assets/137405058/97b1871f-7ea0-4dd2-b289-24afc9373f20)

## Main Developers:

Brayan Alves

Guillaume Briand

Jessie Mosso 

For more details please check [Contributors file](https://github.com/AlvBrayan/MRS4Brain-toolbox/blob/main/CONTRIBUTING)

### Should you publish material that made use of MRS4Brain toolbox, please cite the following publication:

“*Fast high-resolution metabolite mapping in the rat brain using 1 H-FID-MRSI at 14.1T*”

Dunja Simicic, Brayan Alves, Jessie Mosso, Guillaume Briand, Thanh Phong Lê, Ruud B. van Heeswijk, Jana Starčuková, Bernard Lanz, Antoine Klauser, Bernhard Strasser, Wolfgang Bogner, Cristina Cudalbu

*Published in NMR in Biomedicine* : [DOI: 10.1002/nbm.5304](https://doi.org/10.1002/nbm.5304)

“*Noise-reduction techniques for 1H-FID-MRSI at 14.1T: Monte-Carlo validation & in vivo application*”

Brayan Alves, Dunja Simicic, Jessie Mosso, Thanh Phong Lê, Guillaume Briand, Wolfgang Bogner, Bernard Lanz, Bernhard Strasser, Antoine Klauser, Cristina Cudalbu

*Published in NMR in Biomedicine* : [DOI: 10.1002/nbm.5211](https://doi.org/10.1002/nbm.5211)

## Contact, Feedback, Suggestions: 

Thanks for your interest in our MRS4Brain toolbox. For any sort of questions, feedback, suggestions, or critique, please contact us brayan.alves@epfl.ch, cristina.cudalbu@epfl.ch. We also welcome your direct contributions to MRS4Brain toolbox here in the GitHub repository.

## Acknowledgments: 

Financial support was provided by the Swiss National Science Foundation (Project No. 310030_201218).

We acknowledge access to the facilities and expertise of the CIBM Center for Biomedical Imaging founded and supported by Lausanne University Hospital (CHUV), University of Lausanne (UNIL), Ecole polytechnique fédérale de Lausanne (EPFL), University of Geneva (UNIGE) and Geneva University Hospitals (HUG). 

## License:

MRS4Brain toolbox is licensed under [a non-commercial license](https://github.com/AlvBrayan/MRS4Brain-toolbox/blob/main/LICENSE.txt). If you need more information, please contact us brayan.alves@epfl.ch, cristina.cudalbu@epfl.ch 
