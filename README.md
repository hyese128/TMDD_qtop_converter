# TMDD converter : qTMDD to pTMDD
- pTMDD is a novel approximation obtained by Taylor series first-order approximation from qTMDD (QSS approximation model assuming total receptor concentration is constant).
- This new approximation model is accurate compared to the Michaelis-menten model and computationally faster than the Quasi-steady-state model.

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
2. Select the appropriate option for the uploaded model (Refer to the instruction on the first tab of the right panel)
   - Whether TMDD was applied only to drug-receptor interactions or also to drug-FcRn interactions
   - Whether the observed blood concentration is non-binding or total concentration
   - Compartment number where the receptor (drug receptor, FcRn respectively) exists
3. Press the "process file" button to execute the conversion.
4. Click the second tab of the right panel to check the converted pTMDD file.
5. Press the "Download Processed File" button to download the pTMDD file.
