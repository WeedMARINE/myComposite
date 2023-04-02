# myComposite
Laminated composites engineering simulations in Python. Forked from lamipy (https://github.com/joaopbernhardt/lamipy) for personal use in ME227 Class

<!-- <img src="https://upload.wikimedia.org/wikipedia/commons/1/13/Composite_3d.png" data-canonical-src="https://upload.wikimedia.org/wikipedia/commons/1/13/Composite_3d.png" width="300" height="200" /> -->

## Brief description:

Forked from lamipy (https://github.com/joaopbernhardt/lamipy) for personal use in ME227 Class

This project's purpose is to provide simple computations for engineering simulations of **fiber-reinforced composites**. This kind of advanced material is used in many areas of engineering: spacecrafts, pressure vessels, risers etc. Basically, the idea behind composite laminates is to produce a component which has the required engineering properties (high elastic modulus, low self-weight and others) through stacking many layers of laminae, in such way that each one can be layed in different fiber angles and can consist of different materials.

Therefore, a composite laminate consists of a material resultant of layers (laminae) bonded together. In this program, the layers are considered continuous and *orthotropic*. It is possible to input different properties for each used material in the layers.

Currently, implementation of micromechanics and macromechanics for lamina from Robert M. Jones's Mechanics of Composite Material Book is finished.

I planned to incoporate implementation uses the **Classical Laminate Theory (CLT)** for computations. 
The summary of this theory can be found in: [NASA's Basic Mechanics of Laminated Composite Plates](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19950009349.pdf), [Prof. W. Stein's document](http://wstein.org/edu/2010/480b/projects/05-lamination_theory/A%20summary%20of%20Classical%20Lamination%20Theory.pdf) and other sources.

<!-- ## Project goals:

- Provide a simple interface for laminate testing;
- Display accurate information of results, including graphical analysis;

## Technical features:
- [x] Calculates CLT stresses & strains (lamina & laminate coord. systems);
- [x] Failure tests for individual laminas (Tsai-Wu, Max. Stress, Max Strain and Hashin criteria);
- [x] Progressive failure analysis using *Ply Discount* method;
- [x] Plotting of progressive failure analysis;
- [ ] Easy entry of laminate data;
- [ ] Puck failure criteria;
- [x] Thermal & moisture effects on CLT calculations;
- [ ] Laminate materials simple database
- [ ] ...

## How to use:

Currently, *lamipy* is **not ready for general usage**. If you want to test the code, here are the steps:

*Note: using lamipy requires a python3 environment for executing the **runfailuretest.py** file. It is also required **numpy** and **matplotlib** libraries (both of which are easily encountered in many python distributions).*
1. Clone this repository;
1. Edit runfailuretest.py to input the laminate data;
1. Execute.

## Understanding lamipy:

As previously stated, *lamipy* works through the CLT computations. This way, the user input consists of the material properties, layup configurations, force vectors and other settings.

Below is a summary of the data flow inside lamipy for a *simple analysis* (i.e. not a *progressive failure analysis*).
<img src="docs/dataflow.PNG" data-canonical-src="docs/dataflow.PNG" />

### Example results:
Through the obtained results from an analysis, it is possible to plot charts for better visualisation:
<img src="docs/example_plot1.png" data-canonical-src="docs/example_plot1.png" />
<img src="docs/example_plot2.png" data-canonical-src="docs/example_plot2.png" />
<img src="docs/example_plot3.png" data-canonical-src="docs/example_plot3.png" />

Also, since lamipy is capable of running a *progressive failure analysis*, it is possible to automatically plot the Load Factor vs. Average Strain of the laminate while directly pointing *First Ply Failure* and *Last Ply Failure*:
<img src="docs/example_plot4.png" data-canonical-src="docs/example_plot4.png" /> -->