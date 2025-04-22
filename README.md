# AF2χ

## 🏠	 About
Code and data used for the AF2χ project ([Cagiada M., Thomasen F.E., et al., bioRxiv 2025](https://www.biorxiv.com)).  AF2χ is a software for predicting side-chain heterogeneity using AlphaFold2 and its internal side-chain representations. AF2χ outputs side-chain χ-angle distributions and a structural ensemble around the predicted AF2 structure.

## 🖥️	 AF2χ implementation

AF2χ is currently available for the Linux distribution of localColabFold, you can find the implementaiton and all the information on how to use it [here](https://github.com/matteo-cagiada/af2chi_localcolabfold).

----
## 🗂️	Repository layout

The contents of this repository allow you to reproduce the results and figures that appear in the AF2χ manuscript.
- `notebooks`: Folder containing all the Jupyter notebooks used to generate the analysis and figures in the manuscript.
- `data`: Folder containing all the data necessary to run the analysis notebook and reproduce the manuscript results.
- `figures`: collection of the output figures from the analysis notebooks.

**N.B.: Due to their large size, the MDatlas predictions are not included in the current repository. If you wish to reproduce the analyses in full, please contact one of the authors..**

---- 
## 📝 Reference this work

If you use our model please cite:

Cagiada, M., Thomasen, F.E., Ovchinnikov S., Deane C.M &  Lindorff-Larsen, K. (2025). AF2χ: Predicting protein side-chain rotamer distributions with AlphaFold2. In bioRxiv (p. 2024.05.21.595203). https://doi.org/10.1101/2024.05.21.595203

```text
@ARTICLE{Cagiada2025-ax,
  title    = "AF2χ: Predicting protein side-chain rotamer distributions with AlphaFold2",
  author   = "Cagiada, Matteo and Thomasen, F. Emil and Ovchinnikov, Sergey and Deane, Charlotte M. and Lindorff-Larsen, Kresten",
  journal  = "bioRxiv",
  pages    = "",
  month    =  ,
  year     =  ,
  language = "en"
```

Also if you use this localColab implementation remember to cite:

- Mirdita M, Schütze K, Moriwaki Y, Heo L, Ovchinnikov S and Steinegger M. ColabFold - Making protein folding accessible to all. Nature Methods (2022) doi:[ 10.1038/s41592-022-01488-1]( 10.1038/s41592-022-01488-1)
- If you’re using AlphaFold, please also cite:
  Jumper et al. "Highly accurate protein structure prediction with AlphaFold." Nature (2021) doi: [10.1038/s41586-021-03819-2](10.1038/s41586-021-03819-2)
- If you’re using AlphaFold-multimer, please also cite:
  Evans et al. "Protein complex prediction with AlphaFold-Multimer." BioRxiv (2022) doi: [10.1101/2021.10.04.463034v2](10.1101/2021.10.04.463034v2)

## 🙌 Acknowledgements  

The research was supported by the PRISM (Protein Interactions and Stability in Medicine and Genomics) centre funded by the Novo Nordisk Foundation (NNF18OC0033950, to K.L.-L.), a Novo Nordisk Foundation Postdoctoral Fellowship (NNF23OC0082912; to MC). 
We acknowledge access to computational resources via a grant from the Carlsberg Foundation (CF21-0392; to K.L.-L.).

----

## 📜 License  
This project is licensed under the **MIT License**. See [LICENSE](LICENSE) for details.

## 📬 Contact  
For questions or support with this repository, please use the GitHub issue tab or reach out to us via email:  

📧 Matteo Cagiada: [matteo.cagiada@bio.ku.dk](mailto:matteo.cagiada@bio.ku.dk)

📧 Emil Thomasen: [fe.thomasen@bio.ku.dk](mailto:fe.thomasen@bio.ku.dk)

📧 Kresten Lindorff-Larsen: [lindorffb@io.ku.dk](mailto:lindorff@bio.ku.dk)
