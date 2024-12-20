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
1. The interface consists of a left sidebar and a right main panel.
   - Left sidebar
     - Upload file button : can uplod '.mod' files written NONMEM language
    

3. Select the desired ID and click "Load" at the top to fetch the associated data.
