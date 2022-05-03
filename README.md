# File Structure
- Sim.R - Runs basic simulation, no need to edit
- preamble.tex - Basic Latex script for helping display, again do not edit
- multiSim.Rmd - File to run simulations, this is where one can run simulations </br>
Note, all files can be opened with a basic text editor e.g. notepad

# multiSim.Rmd
- R sometimes is strange about where it sources files from, if the line 'source("sim.R")$value' throws an error, try changing it to 'source("sim.R", local=TRUE)$value'
- Do not change the default parameter list, instead set custom parameters at the bottom of the file
- Each instance of the combineFinal function being called is an example of choices that can be made, and meant to demonstrate how one can redefine each parameter
- For running your own simulations, I would recommend deleting the example ones for the sake of brevity
- Note R code must be inside an 'R chunk', which are indicated by the triple backticks and {r}, I would recommend re-using one of the example chunks
- In case it is not clear, a simulation takes the basic format of:
</br> combineFinal(list( 
    </br>$\qquad$ "name for first simulation" = findAvg(REDEFINE PARAMETERS HERE),
    </br>$\qquad$ "name for second simulation" = findAvg(REDEFINE PARAMETERS HERE),
    </br>$\qquad$ ...
</br>))
- One only needs to redefine the parameters that one wants to change
- For the sake of brevity, one can also lower the number of times a simulation is re-ran to find the average, the default being 10
