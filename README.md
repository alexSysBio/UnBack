# UnBack
This repository includes methodologies for subtracting the background in fluorescence microscopy images.

 Uneven_background_correction:
    This script includes the functions required for a cell-free background estimation and subtraction. This is a fast implementation where a the cell-free background is calculated within square image sectors (tiles). These squared regions have a side which is a perfect divisor of the image dimensions. An alternative implementation, which is however slower, would involve a rolling window with a step of 1-pixel. An example for implementing background correction is provided in the test_uneven_background_correction.ipynb notebook.

 Cite:
    
    Nonequilibrium polysome dynamics promote chromosome segregation and its coupling to cell growth in Escherichia coli
    
    Alexandros Papagiannakis, Qiwei Yu, Sander K. Govers, Wei-Hsiang Lin,  Ned S. Wingreen, Christine Jacobs-Wagner
    
    eLife2025;14:RP104276 DOI: https://doi.org/10.7554/eLife.104276.3



