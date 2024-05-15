# DeepSlice_PrV_adjustment

### Background
DeepSlice is a tool used to align brain histology (coronal sections primarily) with relevant atlases. It is primarily used for mice and with the Allen Brain Atlas common coordinate framework (CCF). It is fully integrated with the [QUINT workflow](https://www.ebrains.eu/tools/quint-workflow), which is used for quantitative analysis of brainwide datasets. For more details, check their [GitHub](https://github.com/PolarBean/DeepSlice). Disclaimer: I am not affiliated with DeepSlice in anyway, but have used their [web tool](https://www.deepslice.com.au/), developed by [Michael Pegios](https://github.com/ThermoDev/). <br>

I work with prairie voles (*Microtus ochrogaster*), which only recently had an atlas published ([Gustison et al., 2024](https://doi.org/10.7554/eLife.87029.3)). I received a test version of the QuickNII tool from EBRAINS that uses the vole atlas, and wanted to integrate DeepSlice to speed my work along, as it seems to dramatically help in the early stages of the workflow. Instead of training my own version of DeepSlice on vole brains, I attempted to find a way to easily convert from ABA mouse coordinates to the vole atlas coordinate system. By running downsampled versions of the vole reference brain through DeepSlice web, I was able to establish a *rough* conversion of depth (along the AP axis) betweeen the two atlases. My code uses the predicted depths of sample sections from the mouse atlas to interpolate the corresponding coordinate information from the vole referenec brain. It also adjusts the position and orientation of the section to register it reasonably well, and then exports adjusted JSON and XML files that can be further adjusted in QuickNII, the next step of the QUINT workflow.<br>

### Overview
This GitHub contains a jupyter notebook and the files I used for an example walkthrough of my process. *It's by no means perfect so if anyone has ideas for improvements and are a better coder than me, please let me know!*
* DeepSlice_PrV_Walkthrough.ipynb: JupyterNotebook containing
  * all of the python functions (custom and adapted from [DeepSlice](https://github.com/PolarBean/DeepSlice)) used to adjust DeepSlice alignments to the vole atlas
  * a walkthrough of this process, starting from DeepSlice alignment to adjustments in QuickNII
* PrV_Reference/
  * DeepSlice_PrV_Reference.csv: a combined csv file of DeepSlice's alignment of the vole reference brain ([Gustison et al., 2024](https://doi.org/10.7554/eLife.87029.3)) to the Allen Brain Atlas
  * PrV_DeepSlice_vole.xml: an xml version of the above csv (this is redundant to the csv but thought I'd include it anyway)
  * PrV_DeepSlice_vole.json: a json version of the above csv (this is redundant to the csv but thought I'd include it anyway)
* Sample_Images/
  * Sample_10x_CH2-downsized.png: png of the example section autofluorescence, downsampled to 300px wide for ease in DeepSlice
  * Sample_10x_CH2_DeepSlice_raw.csv: DeepSlice output csv of the example image
  * DeepSlice_PrV_Sample_10x_CH2_adjusted.xml: resulting xml of the example csv after it has been adjusted for vole atlas space. This can be used in QuickNII
  * DeepSlice_PrV_Sample_10x_CH2_adjusted.json: resulting xml of the example csv after it has been adjusted for vole atlas space. This can be used in QuickNII

### Sources
**PrV Reference brain and associated atlas:**<br>
* Gustison, M.L., Muñoz-Castañeda, R., Osten,P., and Phelps, S.M. (2024) **Sexual coordination in a whole-brain map of prairie vole pair bonding** *eLife* **12**:RP87029.
https://doi.org/10.7554/eLife.87029.3<br>
* Data repository: https://figshare.com/articles/journal_contribution/Data_and_code_for_Sexual_coordination_in_a_whole-brain_map_of_prairie_vole_pair_bonding_/21375666/2<br>

**DeepSlice:**<br>
* Carey, H., Pegios, M., Martin, L. et al. (2023) **DeepSlice: rapid fully automatic registration of mouse brain imaging to a volumetric atlas.** *Nat Commun* **14**, 5884 (2023). https://doi.org/10.1038/s41467-023-41645-4.<br>
* GitHub: https://github.com/PolarBean/DeepSlice<br>
* DeepSlice Web: https://www.deepslice.com.au/<br>



