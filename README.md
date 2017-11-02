# Matching Patches for Pollen Identification through Spatially-aware Dictionary Learning and Coding

[project page](http://www.ics.uci.edu/~skong2/recurrentDepthSeg)

![](https://github.com/aimerykong/PatchMatchingForPollenIdentification/blob/master/figures/example_demo.png)


![](https://github.com/aimerykong/PatchMatchingForPollenIdentification/blob/master/figures/patchMatch_demo.png)

  critchfieldii             |     glauca          |     mariana
:-------------------------:|:-------------------------:|:-------------------------:  
![](https://github.com/aimerykong/PatchMatchingForPollenIdentification/blob/master/figures/patches_critchfieldii_K300L0.1_D0.1_E200_B2_globalContrastNorm.png)  |  ![](https://github.com/aimerykong/PatchMatchingForPollenIdentification/blob/master/figures/patches_glauca_K300L0.1_D0.1_E200_B2_globalContrastNorm.png)   |  ![](https://github.com/aimerykong/PatchMatchingForPollenIdentification/blob/master/figures/patches_mariana_K300L0.1_D0.1_E200_B2_globalContrastNorm.png)


This repository contains scripts for generating exemplar patch candidates in training dictionary, 
how to do viewpoint aligment, and the whole classification pipeline for our project with the following published paper.
 
    @inproceedings{kong2016spatially,
      title={Spatially aware dictionary learning and coding for fossil pollen identification},
      author={Kong, Shu and Punyasena, Surangi and Fowlkes, Charless},
      booktitle={Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition Workshops},
      year={2016}
    }


The libsvm toolbox is required. Please compile it under ``libs'' directory before running the code.

The dataset is large; to download the dataset, 
please go [google drive](https://drive.google.com/folderview?id=0BxeylfSgpk1Mdk1HeVhaaEdxMEk&usp=sharing).
To run the demos, please put the downloaded folder under proper directory (refer to addpath in the script).
Note that the fossilized pollen grains are not released for now. Please stay tuned.

Here stores the fossil pollen data that have not been released yet [google drive](https://drive.google.com/drive/folders/0B6uW-Khc9uCDTGk0MUFSekJscWM?usp=sharing)

Here are some useful comments.

1. Folder ``part1_trainval'' contains all the scripts for training, evaluating and visualizing the dictionary.
2. Folder ``figures'' contains the figures used in the paper and others of interest.
3. Folder ``part2_test'' contains the scripts for testing, running the model on a holdout testing set without annotation. This test set can be downloaded in google [google drive](https://drive.google.com/drive/folders/0B6uW-Khc9uCDTGk0MUFSekJscWM?usp=sharing) (permission required). The scripts are cleaner than those in ``part1_trainval''.
4. The two demo directories contains useful scripts for visualization.



For questions, please contact
 
 Shu Kong (Aimery) aimerykong AT gmail com

