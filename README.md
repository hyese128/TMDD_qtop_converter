# TMDD converter : qTMDD to pTMDD
- pTMDD is a novel approximation obtained by Taylor series first-order approximation from qTMDD (QSS approximation model assuming total receptor concentration is constant).
- This new approximation model is accurate compared to the Michaelis-menten model and computationally faster than the Quasi-steady-state model.
- For a detailed description of pTMDD, refer to the following DOI: 10.1371/journal.pcbi.1012066

## Features
1. Web-based (R-shiny)
2. Translates qTMDD models written in NONMEM language to pTMDD models.
3. Needed libraries
   - dplyr
   - stringr
   - tidyselect
   - purrr
   - shinythemes

## Instruction
1. Upload file tab : can uplod '.mod' files written NONMEM language
<div align="center">
    <img src="https://github.com/hyese128/TMDD_qtop_converter/blob/main/images/image1.png" width="400" height="300">
</div>

2. Select the appropriate option for the uploaded model (Refer to the instruction on the first tab of the right panel)
<div align="center">
    <img src="https://github.com/hyese128/TMDD_qtop_converter/blob/main/images/image2.png" width="400" height="300">
</div>
   - Whether TMDD was applied only to drug-receptor interactions or also to drug-FcRn interactions
   - Whether the observed blood concentration is non-binding or total concentration
   - Compartment number where the receptor (drug receptor, FcRn respectively) exists
   
3. Press the "process file" button to execute the conversion.
<div align="center">
    <img src="https://github.com/hyese128/TMDD_qtop_converter/blob/main/images/image3.png" width="400" height="300">
</div>

4. Click the second tab of the right panel to check the converted pTMDD file.
<div align="center">
    <img src="https://github.com/hyese128/TMDD_qtop_converter/blob/main/images/image4.png" width="400" height="300">
</div>

5. Press the "Download Processed File" button to download the pTMDD file.
<div align="center">
    <img src="https://github.com/hyese128/TMDD_qtop_converter/blob/main/images/image5.png" width="400" height="300">
</div>
