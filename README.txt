# Fick's Second Law Calculator

**Author:** Pranav Roy  
**Created:** August 2025

© 2025 Pranav Roy. All rights reserved.

This project and its code are original works by Pranav Roy.  
If you use or share this code, please credit the author.

---

## Contact

- **Email:** pranavroy282@gmail.com
- **LinkedIn:** www.linkedin.com/in/pranav-roy-64314b184

---

## Overview

This project provides two Python-based calculators for Fick's Second Law of Diffusion:
- A command-line calculator (`Fick_s 2nd Law Nonsteady-State Diffusion Calculator.py`)
- A graphical user interface (GUI) calculator (`Fick_s 2nd Law Nonsteady-Diffusion Calculator_gui.pyw`)

Both tools help you solve for:
- Cxt: Concentration at position x and time t [wt%]
- C0: Initial concentration [wt%]
- Cs: Surface concentration [wt%]
- x: Distance [m]
- D: Diffusion coefficient [m^2/s]
- t: Time [s]

with built-in error handling, logging, and unit conversion.

---

## Scientific Background

**Fick's 2nd Law of Diffusion:**
Describes nonsteady-state diffusion, where flux and concentration gradient change with time.

**General Equation:**
	∂C/∂t = ∂/∂x ( D ∂C/∂x )
If D is constant:
	∂C/∂t = D ∂²C/∂x²

**Common Solution (semi-infinite solid, constant surface concentration):**
Assumptions:
1. Before diffusion, solute atoms are uniformly distributed with concentration C₀.
2. x = 0 at the surface, increases into the solid.
3. t = 0 is the instant before diffusion begins.

Initial condition:   For t = 0, C = C₀ for 0 ≤ x < ∞
Boundary conditions: For t > 0, C = Cs at x = 0 (surface)
					For t > 0, C = C₀ at x = ∞

**Solution:**
	(Cx - C₀) / (Cs - C₀) = 1 - erf( x / (2√(Dt)) )

**Where:**
- Cx = Concentration at position x and time t
- C₀ = Initial concentration
- Cs = Surface concentration
- x  = Distance from surface
- D  = Diffusion coefficient
- t  = Time
- erf = Error function

---

## Features

### Command-line calculator (`Fick_s 2nd Law Nonsteady-State Diffusion Calculator.py`)
- Interactive prompts for variable selection and input
- Error handling for invalid inputs and divide-by-zero cases
- Logs calculations with date and time in `ficks_2nd_law_log.txt`
- Uses tabulated error function (erf) values for accurate interpolation

### GUI calculator (`Fick_s 2nd Law Nonsteady-Diffusion Calculator_gui.pyw`)
- Intuitive interface using Tkinter
- Select which variable to solve for
- Only relevant input fields are shown
- Buttons for calculation, clearing, and quitting
- Unit converter for x, D, t, Qd, T, R (supports multiple units)
- Arrhenius equation calculator for temperature-dependent diffusion
- Error messages for missing/invalid input and divide-by-zero
- Scientific explanation panel for Fick's Law and variables

---

## Getting Started

### Prerequisites
- Python 3.x installed on your system
- Tkinter (included with standard Python installations)

### Running the Command-Line Calculator
1. Open a terminal and navigate to the project folder.
2. Run:
   ```
   python "Flick_s Law Calculators/2nd Law/Fick_s 2nd Law Nonsteady-State Diffusion Calculator.py"
   ```
3. Follow the prompts to select a variable and enter values.

### Running the GUI Calculator
1. Open a terminal and navigate to the project folder.
2. Run:
   ```
   python "Flick_s Law Calculators/2nd Law/Fick_s 2nd Law Nonsteady-Diffusion Calculator_gui.pyw"
   ```
3. The GUI window will appear. Select the variable to solve for, enter the required values, and click **Calculate**. Use the unit converter and Arrhenius calculator as needed.

---

## Logging
- All calculations in the command-line version are logged to `ficks_2nd_law_log.txt` with date and time.
- You can use this log for record-keeping or analysis.

---

## License
This project is provided for educational and research purposes.

---

**Feel free to modify and expand the calculator for your needs!**
