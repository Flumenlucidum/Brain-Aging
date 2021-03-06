# Brain Aging

Investigation of Genetic Variants and Causal Biomarkers Associated with Brain Aging
Jangho Kim, Junhyeong Lee, Seunggeun Lee
medRxiv 2022.03.04.22271813; doi: https://doi.org/10.1101/2022.03.04.22271813

### Overview of the Analysis
<div>
          <p align="left">
          <img width = "600" height = "350" src = "https://github.com/Flumenlucidum/Brain-Aging/blob/main/images/method.jpg">
</div>

### 3D Convlutional Neural Network Model for Age Prediction 
<div>
          <p align="left">
          <img width = "600" height = "350" src = "https://github.com/Flumenlucidum/Brain-Aging/blob/main/images/3dcnn.jpg">
</div>

The model takes 3D structural MRI of brain and predicts age of the individual. The mean absolute error (MAE) of our model was
2.64 years. We calculated 'delta age' by subtracting individuals' true age from the predicted age

### Saliency Map derived by Integrated Gradients
<div>
          <p align="center">
          <img width = "600" height = "350" src = "https://github.com/Flumenlucidum/Brain-Aging/blob/main/images/final_plot.jpg">
</div>

Integrated Gradients showed that the prediction model highlighted the fornix and the lower part of the thalamus.

### Genome-wide Association Tests on Single Variants (SAIGE) and Gene Regions (SKAT-CommonRare)

<div>
          <p align="left">
          <img width = "600" height = "350" src = "https://github.com/Flumenlucidum/Brain-Aging/blob/main/images/final.jpg">
</div>
<div>
          <p align="left">
          <img width = "350" height = "350" src = "https://github.com/Flumenlucidum/Brain-Aging/blob/main/images/delta_0227.jpeg">
</div>

The genetic variants associated with delta age had links to carcinogenesis, immune response, and neuron survival.
In addition, the volume of the fornix and thalamus also had similar associated genetic variants. Delta age and the volume of the two regions had high genetic correlation in LD Score regression.

### Causal Biomarkers 
<div>
          <p align="left">
          <img width = "350" height = "350" src = "https://github.com/Flumenlucidum/Brain-Aging/blob/main/images/plot_Egger.jpeg">
</div>

Mendelian randomization is a causal inference method using the genotype as an instrument variable.
Eosinophil count, Eosinophil percentage, Neutrophil count, and Total protein had a significant causal effect on delta age.
They all increased delta age. 

<div>
          <p align="left">
          <img width = "350" height = "350" src = "https://github.com/Flumenlucidum/Brain-Aging/blob/main/images/23436_pl.png">
          <img width = "350" height = "350" src = "https://github.com/Flumenlucidum/Brain-Aging/blob/main/images/23436_kiv.png">
</div>

Nonlinear causal patterns were derived using piecewise MR and Kernel IV regression.
Three biomarkers had p-values less than 0.05, but they were insignificant after multiple testing correction.
The above plots are from causal estimates of total cholines.
