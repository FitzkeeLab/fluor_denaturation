# fluor_denaturation
Python Script to Fit Fluorescence Denaturation Data

This project consits of a Python 3 script (fitting.py) and an associated Excel Spreadsheet (Protein Denaturation Worksheet.xlsx). The worksheet is used to take fluorescence spectrum from Horiba and prepare them for denaturation measurements. This includes substracting baseline and normalizing the spectra. The worksheet will also calculate denaturant concentrations using a dilution equation and it will correct the spectra for dilution, assuming that the Fluoromax titration software is used in constant volume mode. The first workbook contains a table of two values: denaturant concentration and fluorescence. This table is used as input for the Python script; however, the text file of concentration vs fluorescence can be generated using any means.

The Python script requires the following libraries:

1. Numeric Python (numpy)
2. Scientific Python (scipy)
3. Matplotlib

These dependencies can be easily installed using the Package Installer for Python (pip):

```
pip install --user numpy
pip install --user scipy
pip install --user matplotlib
```

## Fitting Script

This program is fairly basic and requires some degree of user interaction to run. The name of the titration curve must be supplied to the program in the main function:

```
def main():
    data_file = "R2ab-Y722A-Den.txt"
    n_bootstrap = 1000
    verbose = False

    pnames     = ['af',      'bf',  'au',    'bu', 'dG',  'm']
    parameters = [-0.3, 1.05, -0.001, 0.3, 3.0,  4.0]

...
```

The titration file should be set to the data_file variable. The initial guesses for fit parameters (slope and intercept for the folded and unfolded baselines, delta G, and m-value) are listed under parameters.

## Workflow

I recommend creating a new folder and using a new fitting.py file for every fit that you do. This will keep the code together and you can go back to the initial guesses that led to the resulting fit.

To run the code, copy your data file into the folder and set the variable in fitting.py to read this file. Then run the code. You will be presented a graph where *no fitting is done*. This is simply to show you how good your initial guesses are. At this point, you can edit the initial guesses (the values in the parameters list above) and try to get your fit as close to the data points as possible. When you make the changes, quit the python code (ctrl-C) and restart it to check the new fit parameters.

Once your fit is close to the datapoints, you can close the graph window and the program will do an initial fit.

Then, it will perform a series of bootstrap calculations, resampling your original data points (with replacement) and redoing the fit, according to the number of loops in n_bootstrap. Some of these fits will fail; that's okay.

## Output

When the script completes, you will be presented with a table showing the best fit values as well as the bootstrap results. Bootstrap values are displayed as the average and then positive and negative confidence intervals. This table is stored as par_summary.txt. A histogram of bootstrap values is provided for each parameter, both as a PDF and a text file of all the values.

The program also constructs a normalized curve showing the fraction folded at each denaturant value. The residuals are used to transform the raw data onto the fraction folded curve. This transformation allows many folding curves to be plotted on the same graph. The transformed curve and transformed data are in trans_curve.txt and trans_data, respectively, and a visual graph is stored in transform.pdf.

The fit.pdf and fit_curve.txt files contain the best fit curve through the original input data in data_file. The residuals.txt file allows you to check the goodness of fit.

Note that, when reporting best fit values, you should report the overall best fit - not the average of the bootstrap fits, which is just supplied for validation and may be subject to skews from leaving data points out. While bootstrap error bars are useful, it is always better to present the error bars based on independently prepared experiments. However, in many cases we find the bootstrap uncertainties and the experimental uncertainties to be comparable.

## Citing this Work

If you publish a result using this tool, please cite the following article:

> Kariyawasam C.S., Somarathne R.P., Mayatt R.S., Conner R.A., Fitzkee N.C.. "Thermodynamic Analysis of the Nanoparticle Corona Reveals a Correlation Between Binding Affinity and Structural Stability." *Submitted.*

