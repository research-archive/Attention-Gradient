Data_Codes – EJN Manuscript (Revision 1)

Overview

The figures and the statistical results reported in the manuscript (EJN) are organized into separate subfolders within the main
folder named:

    Data_Codes

Each figure and the analyses associated with it has its own dedicated subfolder. These subfolders contain:

    MATLAB scripts for generating the figures
    Associated .mat data files
    Statistical analyses corresponding to the results reported in the manuscript

How to Run the Code

To generate the figures and reproduce the analyses:

1.  Open MATLAB

2.  Navigate to the specific figure subfolder inside Data_Codes using
    the MATLAB command window. For example:

        cd('path_to_repository/Data_Codes/FigureX')

    Replace “FigureX” with the appropriate figure folder (e.g., Figure2,
    Figure3, Figure4).

3.  Run the main script within that folder:

        run('FigureX_script.m')

    The script will: Generate the respective figure panels; Reproduce
    the statistical analyses; Print relevant statistical outputs in the
    MATLAB Command Window

Notes
For any questions regarding data organization or code execution, please
refer to the corresponding author of the manuscript.
