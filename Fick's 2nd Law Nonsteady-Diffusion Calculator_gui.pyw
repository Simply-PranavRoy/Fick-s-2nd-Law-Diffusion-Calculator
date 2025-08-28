# Fick's 2nd Law Nonsteady-State Diffusion Calculator (GUI)
# ----------------------------------------------------------
# GUI tool for calculating nonsteady-state diffusion using Fick's Second Law.
# Features: variable solving, dynamic input fields, error handling.

import os
os.environ['TK_SILENCE_DEPRECATION'] = '1'
import tkinter as tk
from tkinter import ttk, messagebox
import math

# Tabulated z and erf(z) values from research paper
ERF_TABLE = [
    [0.00,  0.0000], [0.025, 0.0282], [0.05,  0.0564], [0.10,  0.1125],
    [0.15,  0.1680], [0.20,  0.2227], [0.25,  0.2763], [0.30,  0.3286],
    [0.35,  0.3794], [0.40,  0.4284], [0.45,  0.4755], [0.50,  0.5205],
    [0.55,  0.5633], [0.60,  0.6039], [0.65,  0.6420], [0.70,  0.6778],
    [0.75,  0.7112], [0.80,  0.7421], [0.85,  0.7707], [0.90,  0.7970],
    [0.95,  0.8209], [1.00,  0.8427], [1.10,  0.8802], [1.20,  0.9103],
    [1.30,  0.9340], [1.40,  0.9523], [1.50,  0.9661], [1.60,  0.9763],
    [1.70,  0.9838], [1.80,  0.9891], [1.90,  0.9928], [2.00,  0.9953],
    [2.20,  0.9981], [2.40,  0.9993], [2.60,  0.9998], [2.80,  0.9999]
]

# --- GUI Setup ---
root = tk.Tk()
root.title("Fick's 2nd Law Nonsteady-State Diffusion Calculator")
root.geometry("900x700")
mainframe = ttk.Frame(root, padding="24 24 24 24")
mainframe.grid(row=0, column=0, sticky='nw')

input_frame = ttk.Frame(mainframe, padding="12 12 12 12")
input_frame.grid(row=0, column=0, sticky='nw')

# --- Variable Selection ---
row = 0
ttk.Label(input_frame, text="Solve for:").grid(row=row, column=0, sticky=tk.W)
variable = tk.StringVar()
solve_options = ['Cxt', 'C0', 'Cs', 'x', 'D', 't']
solve_menu = ttk.Combobox(input_frame, textvariable=variable, values=solve_options, state='readonly', width=12)
solve_menu.grid(row=row, column=1, sticky=tk.W)
solve_menu.current(0)

# --- Input Fields ---
fields = {
    'Cxt': "Cxt [wt%]", 'C0': "C0 [wt%]", 'Cs': "Cs [wt%]", 'x': "x [m]", 'D': "D [m^2/s]", 't': "t [s]"
}
entries = {}
input_widgets = {}
row = 1
for key, label in fields.items():
    ttk.Label(input_frame, text=label).grid(row=row, column=0, sticky=tk.W)
    entry = ttk.Entry(input_frame, width=15)
    entry.grid(row=row, column=1, sticky=tk.W)
    entries[key] = entry
    input_widgets[key] = entry
    row += 1

# --- Required Fields for Each Calculation ---
required_fields = {
    'Cxt': ['C0', 'Cs', 'x', 'D', 't'],
    'C0': ['Cxt', 'Cs', 'x', 'D', 't'],
    'Cs': ['Cxt', 'C0', 'x', 'D', 't'],
    'x': ['Cxt', 'C0', 'Cs', 'D', 't'],
    'D': ['Cxt', 'C0', 'Cs', 'x', 't'],
    't': ['Cxt', 'C0', 'Cs', 'x', 'D']
}

def show_fields(*args):
    for widget in input_widgets.values():
        widget.grid_remove()
    target = variable.get()
    for field in required_fields[target]:
        input_widgets[field].grid()
solve_menu.bind("<<ComboboxSelected>>", show_fields)

# --- Calculation Logic ---
def calculate():
    target = variable.get()
    vals = {}
    try:
        for field in required_fields[target]:
            val = entries[field].get()
            if val == "":
                messagebox.showerror("Input Error", f"Please enter {fields[field]}.")
                return
            vals[field] = float(val)
        # Calculation for each variable
        if target == 'Cxt':
            C0, Cs, x, D, t = vals['C0'], vals['Cs'], vals['x'], vals['D'], vals['t']
            if Cs == C0:
                messagebox.showerror("Math Error", "(Cs - C0) cannot be zero (division by zero).")
                return
            denom = 2 * math.sqrt(D * t)
            if denom == 0:
                messagebox.showerror("Math Error", "(2 * sqrt(D * t)) cannot be zero.")
                return
            z = x / denom
            if not (0 <= z <= 2.8):
                messagebox.showerror("Range Error", "z must be between 0 and 2.8.")
                return
            erf_z = round(math.erf(z), 4)
            if not (0 <= erf_z <= 0.9999):
                messagebox.showerror("Range Error", "erf(z) must be between 0 and 0.9999.")
                return
            Cxt = C0 + ((Cs - C0) * (1 - erf_z))
            result_var.set(f"z = {z:.4e}\nerf(z) = {erf_z:.4f}\nCxt [wt%] = {Cxt:.4e}")
        elif target == 'C0':
            Cxt, Cs, x, D, t = vals['Cxt'], vals['Cs'], vals['x'], vals['D'], vals['t']
            denom = 2 * math.sqrt(D * t)
            if denom == 0:
                messagebox.showerror("Math Error", "(2 * sqrt(D * t)) cannot be zero.")
                return
            z = x / denom
            if not (0 <= z <= 2.8):
                messagebox.showerror("Range Error", "z must be between 0 and 2.8.")
                return
            erf_z = round(math.erf(z), 4)
            if not (0 <= erf_z <= 0.9999):
                messagebox.showerror("Range Error", "erf(z) must be between 0 and 0.9999.")
                return
            a = 1 - erf_z
            C0 = (Cxt - (a * Cs)) / (1 - a) if (1 - a) != 0 else 0
            if Cs == C0:
                messagebox.showerror("Math Error", "(Cs - C0) cannot be zero (division by zero).")
                return
            result_var.set(f"z = {z:.4e}\nerf(z) = {erf_z:.4f}\nC0 [wt%] = {C0:.4e}")
        elif target == 'Cs':
            Cxt, C0, x, D, t = vals['Cxt'], vals['C0'], vals['x'], vals['D'], vals['t']
            denom = 2 * math.sqrt(D * t)
            if denom == 0:
                messagebox.showerror("Math Error", "(2 * sqrt(D * t)) cannot be zero.")
                return
            z = x / denom
            if not (0 <= z <= 2.8):
                messagebox.showerror("Range Error", "z must be between 0 and 2.8.")
                return
            erf_z = round(math.erf(z), 4)
            if not (0 <= erf_z <= 0.9999):
                messagebox.showerror("Range Error", "erf(z) must be between 0 and 0.9999.")
                return
            Cs = ((Cxt - C0) / (1 - erf_z)) + C0 if z != 1 else None
            if Cs == C0:
                messagebox.showerror("Math Error", "(Cs - C0) cannot be zero (division by zero).")
                return
            result_var.set(f"z = {z:.4e}\nerf(z) = {erf_z:.4f}\nCs [wt%] = {Cs:.4e}")
        elif target == 'x':
            Cxt, C0, Cs, D, t = vals['Cxt'], vals['C0'], vals['Cs'], vals['D'], vals['t']
            if Cs == C0:
                messagebox.showerror("Math Error", "(Cs - C0) cannot be zero (division by zero).")
                return
            cerfz = 1 - ((Cxt - C0) / (Cs - C0))
            lower = None
            higher = None
            for i in range(len(ERF_TABLE)):
                if ERF_TABLE[i][1] <= cerfz:
                    lower = ERF_TABLE[i]
                if ERF_TABLE[i][1] >= cerfz and higher is None:
                    higher = ERF_TABLE[i]
            if lower is None or higher is None:
                messagebox.showerror("Range Error", "Calculated erf(z) is out of table range.")
                return
            lz, lerfz = lower
            hz, herfz = higher
            nz = (((cerfz - lerfz) / (herfz - lerfz)) * (hz - lz)) + lz if herfz != lerfz else lz
            x = nz * 2 * math.sqrt(D * t)
            result_var.set(f"Calculated z = {cerfz:.4f}\nLower: z = {lz}, erf(z) = {lerfz}\nHigher: z = {hz}, erf(z) = {herfz}\nInterpolated z = {nz:.4f}\nD [m^2/s] = {D:.4e}\nt [s] = {t:.4e}\nx [m] = {x:.4e}")
        elif target == 'D':
            Cxt, C0, Cs, x, t = vals['Cxt'], vals['C0'], vals['Cs'], vals['D'], vals['t']
            if Cs == C0:
                messagebox.showerror("Math Error", "(Cs - C0) cannot be zero (division by zero).")
                return
            cerfz = 1 - ((Cxt - C0) / (Cs - C0))
            lower = None
            higher = None
            for i in range(len(ERF_TABLE)):
                if ERF_TABLE[i][1] <= cerfz:
                    lower = ERF_TABLE[i]
                if ERF_TABLE[i][1] >= cerfz and higher is None:
                    higher = ERF_TABLE[i]
            if lower is None or higher is None:
                messagebox.showerror("Range Error", "Calculated erf(z) is out of table range.")
                return
            lz, lerfz = lower
            hz, herfz = higher
            nz = (((cerfz - lerfz) / (herfz - lerfz)) * (hz - lz)) + lz if herfz != lerfz else lz
            D = (x / (2 * nz * math.sqrt(t))) ** 2 if nz != 0 and t > 0 else 0
            if D < 0:
                messagebox.showerror("Math Error", "Calculated D is negative, check your inputs.")
                return
            result_var.set(f"Calculated erf(z) = {cerfz:.4f}\nLower: z = {lz}, erf(z) = {lerfz}\nHigher: z = {hz}, erf(z) = {herfz}\nInterpolated z = {nz:.4f}\nt [s] = {t:.4e}\nx [m] = {x:.4e}\nD [m^2/s] = {D:.4e}")
        elif target == 't':
            Cxt, C0, Cs, x, D = vals['Cxt'], vals['C0'], vals['Cs'], vals['x'], vals['D']
            if Cs == C0:
                messagebox.showerror("Math Error", "(Cs - C0) cannot be zero (division by zero).")
                return
            cerfz = 1 - ((Cxt - C0) / (Cs - C0))
            lower = None
            higher = None
            for i in range(len(ERF_TABLE)):
                if ERF_TABLE[i][1] <= cerfz:
                    lower = ERF_TABLE[i]
                if ERF_TABLE[i][1] >= cerfz and higher is None:
                    higher = ERF_TABLE[i]
            if lower is None or higher is None:
                messagebox.showerror("Range Error", "Calculated erf(z) is out of table range.")
                return
            lz, lerfz = lower
            hz, herfz = higher
            nz = (((cerfz - lerfz) / (herfz - lerfz)) * (hz - lz)) + lz if herfz != lerfz else lz
            t = (x / (2 * nz * math.sqrt(D))) ** 2 if nz != 0 and D > 0 else 0
            if t < 0:
                messagebox.showerror("Math Error", "Calculated t is negative, check your inputs.")
                return
            result_var.set(f"Calculated erf(z) = {cerfz:.4f}\nLower: z = {lz}, erf(z) = {lerfz}\nHigher: z = {hz}, erf(z) = {herfz}\nInterpolated z = {nz:.4f}\nD [m^2/s] = {D:.4e}\nx [m] = {x:.4e}\nt [s] = {t:.4e}")
        else:
            messagebox.showinfo("Info", "Select a variable to solve for.")
    except ValueError:
        messagebox.showerror("Input Error", "Please enter valid numeric values.")

# --- Result Display ---
result_var = tk.StringVar()
result_label = ttk.Label(input_frame, textvariable=result_var, foreground="blue", font=("Segoe UI", 11, "bold"))
result_label.grid(row=row, column=0, columnspan=2, pady=(10,0), sticky='w')

# --- Action Buttons ---
def clear_inputs():
    for entry in entries.values():
        entry.delete(0, tk.END)
    result_var.set("")
# ...existing code...
# --- Action Buttons ---
calc_btn = ttk.Button(input_frame, text="Calculate", command=calculate)
calc_btn.grid(row=row+1, column=0, pady=10, sticky='w')
calc_another_btn = ttk.Button(input_frame, text="Clear", command=clear_inputs)
calc_another_btn.grid(row=row+1, column=1, pady=10, sticky='w')
quit_btn = ttk.Button(input_frame, text="Quit", command=root.quit)
quit_btn.grid(row=row+1, column=2, pady=10, sticky='w')

# --- D Arrhenius Calculator ---
d_arrhenius_frame = tk.LabelFrame(input_frame, text="Arrhenius Equation Calculator", padx=10, pady=10)
d_arrhenius_frame.grid(row=row+2, column=0, columnspan=3, pady=(30,0), sticky='nw')

# Dropdown to select variable to solve for
arrhenius_vars = ['D', 'D₀', 'Qd', 'T']
solve_arrhenius_var = tk.StringVar()
solve_arrhenius_var.set('D')
tk.Label(d_arrhenius_frame, text="Solve for:").grid(row=0, column=0, sticky='w')
solve_arrhenius_menu = ttk.Combobox(d_arrhenius_frame, textvariable=solve_arrhenius_var, values=arrhenius_vars, state='readonly', width=6)
solve_arrhenius_menu.grid(row=0, column=1, sticky='w')

# Input fields
arrhenius_entries = {}
arrhenius_labels = {}
labels = {
    'D': "D (m²/s):",
    'D₀': "D₀ (m²/s):",
    'Qd': "Qd (J/mol):",
    'T': "T (K):",
    'R': "R (J/mol·K):"
}
row_offset = 1
for key in ['D', 'D₀', 'Qd', 'T', 'R']:
    lbl = tk.Label(d_arrhenius_frame, text=labels[key])
    ent = tk.Entry(d_arrhenius_frame, width=15)
    arrhenius_labels[key] = lbl
    arrhenius_entries[key] = ent
    lbl.grid(row=row_offset, column=0, sticky='w')
    ent.grid(row=row_offset, column=1)
    row_offset += 1
arrhenius_entries['R'].insert(0, "8.31")

# Show/hide fields based on variable to solve for
fields_for = {
    'D': ['D₀', 'Qd', 'T', 'R'],
    'D₀': ['D', 'Qd', 'T', 'R'],
    'Qd': ['D', 'D₀', 'T', 'R'],
    'T': ['D', 'D₀', 'Qd', 'R']
}
def update_arrhenius_fields(*args):
    target = solve_arrhenius_var.get()
    # Hide all input fields and labels
    for key in ['D', 'D₀', 'Qd', 'T']:
        arrhenius_entries[key].grid_remove()
        arrhenius_labels[key].grid_remove()
        arrhenius_entries[key].config(state='normal')
    # Show only relevant input fields and labels (except the one being solved for)
    for key in fields_for[target] + ['R']:
        arrhenius_entries[key].grid()
        arrhenius_labels[key].grid()
    # Hide the entry and label for the variable being solved
        arrhenius_entries[target].grid_remove()
        arrhenius_labels[target].grid_remove()

    # Binding and initial call moved outside the function

# Bind the dropdown and call once to set initial state
solve_arrhenius_menu.bind('<<ComboboxSelected>>', update_arrhenius_fields)
# Ensure fields are placed and visible at startup
root.after(100, update_arrhenius_fields)

# Result display
arrhenius_result_var = tk.StringVar()
arrhenius_result_label = tk.Label(d_arrhenius_frame, textvariable=arrhenius_result_var, fg="blue", font=("Segoe UI", 11, "bold"))
arrhenius_result_label.grid(row=row_offset+1, column=0, columnspan=2, pady=(10,0), sticky='w')

# Calculation logic
import math
def calculate_arrhenius():
    target = solve_arrhenius_var.get()
    try:
        vals = {k: float(arrhenius_entries[k].get()) for k in fields_for[target]+['R']}
        if target == 'D':
            D0, Qd, T, R = vals['D₀'], vals['Qd'], vals['T'], vals['R']
            D = D0 * math.exp(-Qd / (R * T))
            arrhenius_result_var.set(f"D = {D:.4e} m²/s")
        elif target == 'D₀':
            D, Qd, T, R = vals['D'], vals['Qd'], vals['T'], vals['R']
            D0 = D / math.exp(-Qd / (R * T))
            arrhenius_result_var.set(f"D₀ = {D0:.4e} m²/s")
        elif target == 'Qd':
            D, D0, T, R = vals['D'], vals['D₀'], vals['T'], vals['R']
            Qd = -(math.log(D / D0)) * R * T
            arrhenius_result_var.set(f"Qd = {Qd:.4e} J/mol")
        elif target == 'T':
            D, D0, Qd, R = vals['D'], vals['D₀'], vals['Qd'], vals['R']
            T = -Qd / (math.log(D / D0) * R)
            arrhenius_result_var.set(f"T = {T:.2f} K")
    except Exception:
        arrhenius_result_var.set("Enter valid numeric values.")

calc_arrhenius_btn = tk.Button(d_arrhenius_frame, text="Calculate", command=calculate_arrhenius)
calc_arrhenius_btn.grid(row=row_offset, column=0, columnspan=2, pady=(10,0))

# --- Info Section ---
info_frame = tk.LabelFrame(mainframe, text="About Fick's 2nd Law", padx=12, pady=12)
info_frame.grid(row=0, column=1, rowspan=10, sticky='nw', padx=(40,0))

info_text = (
    "Fick's Second Law — Nonsteady-State Diffusion\n"
    "\n"
    "Most practical diffusion situations are nonsteady-state: the flux and concentration gradient at a point in a solid vary with time.\n"
    "\n"
    "Fick's Second Law (general form):\n"
    "    ∂C/∂t = ∂/∂x ( D ∂C/∂x )\n"
    "If D is constant, this simplifies to:\n"
    "    ∂C/∂t = D ∂²C/∂x²\n"
    "\n"
    "A common solution is for a semi-infinite solid with constant surface concentration.\n"
    "Assumptions:\n"
    "  1. Before diffusion, solute atoms are uniformly distributed with concentration C₀.\n"
    "  2. x = 0 at the surface, increases into the solid.\n"
    "  3. t = 0 is the instant before diffusion begins.\n"
    "\n"
    "Initial condition:   For t = 0, C = C₀ for 0 ≤ x < ∞\n"
    "Boundary conditions: For t > 0, C = Cs at x = 0 (surface)\n"
    "                    For t > 0, C = C₀ at x = ∞\n"
    "\n"
    "Solution for these conditions:\n"
    "    (Cx - C₀) / (Cs - C₀) = 1 - erf( x / (2√(Dt)) )\n"
    "\n"
    "Where:\n"
    "  Cx = Concentration at position x and time t\n"
    "  C₀ = Initial concentration\n"
    "  Cs = Surface concentration\n"
    "  x  = Distance from surface\n"
    "  D  = Diffusion coefficient\n"
    "  t  = Time\n"
    "  erf = Error function\n"
)
info_label = tk.Label(info_frame, text=info_text, justify='left', font=("Segoe UI", 10), anchor='nw')
info_label.pack(fill='x', expand=False)

# --- Unit Converter Section (inside info box) ---
unit_frame = tk.LabelFrame(info_frame, text="Unit Converter", padx=10, pady=10)
unit_frame.pack(fill='x', pady=(16,0), anchor='w')

# Extensive unit options for D, C, x, t
unit_options = {
    "D": ["m²/s", "cm²/s", "mm²/s", "µm²/s", "ft²/s", "in²/s", "m²/hr", "cm²/hr"],
    "x": ["m", "cm", "mm", "µm", "nm", "ft", "in"],
    "t": ["s", "min", "hr", "day"],
    "Qd": ["J/mol", "kJ/mol", "eV/atom"],
    "T": ["K", "°C", "°F"],
    "R": ["J/(mol·K)", "kJ/(mol·K)", "eV/(atom·K)"]
}

# Conversion factors
length_factors = {
    'm': 1.0, 'cm': 0.01, 'mm': 0.001, 'µm': 1e-6, 'nm': 1e-9, 'ft': 0.3048, 'in': 0.0254
}
area_factors = {
    'm²/s': 1.0, 'cm²/s': 0.0001, 'mm²/s': 1e-6, 'µm²/s': 1e-12, 'ft²/s': 0.092903, 'in²/s': 0.00064516, 'm²/hr': 1.0/3600, 'cm²/hr': 0.0001/3600
}
conc_factors = {
    'mol/m³': 1.0, 'mol/L': 1000.0, 'mol/cm³': 1e6, 'mmol/L': 1.0, 'µmol/L': 1e-3, 'kg/m³': 1.0, 'g/L': 1.0, 'g/cm³': 1000.0, 'mg/mL': 1000.0, 'mg/L': 1e-3, 'µg/mL': 1e-6
}
time_factors = {
    's': 1.0, 'min': 60.0, 'hr': 3600.0, 'day': 86400.0
}
# Qd conversion factors
qd_factors = {
    'J/mol': 1.0,
    'kJ/mol': 1000.0,
    'eV/atom': 96.485e3,  # 1 eV/atom = 96,485 J/mol
}
# T conversion functions
def t_to_K(val, from_unit):
    if from_unit == 'K':
        return val
    elif from_unit == '°C':
        return val + 273.15
    elif from_unit == '°F':
        return (val - 32) * 5/9 + 273.15
    else:
        return val
def K_to_t(val_K, to_unit):
    if to_unit == 'K':
        return val_K
    elif to_unit == '°C':
        return val_K - 273.15
    elif to_unit == '°F':
        return (val_K - 273.15) * 9/5 + 32
    else:
        return val_K
# R conversion factors
r_factors = {
    'J/(mol·K)': 1.0,
    'kJ/(mol·K)': 1000.0,
    'eV/(atom·K)': 96.485e3,
}

# Widgets
tk.Label(unit_frame, text="Value:").grid(row=0, column=0)
value_entry = tk.Entry(unit_frame)
value_entry.grid(row=0, column=1)

tk.Label(unit_frame, text="Variable:").grid(row=1, column=0)
variable_combo = ttk.Combobox(unit_frame, values=list(unit_options.keys()), state='readonly')
variable_combo.grid(row=1, column=1)

from_unit_combo = ttk.Combobox(unit_frame, state='readonly')
from_unit_combo.grid(row=2, column=1)
to_unit_combo = ttk.Combobox(unit_frame, state='readonly')
to_unit_combo.grid(row=3, column=1)

def update_unit_dropdowns(event=None):
    var = variable_combo.get()
    units = unit_options.get(var, [])
    from_unit_combo['values'] = units
    to_unit_combo['values'] = units
    if units:
        from_unit_combo.current(0)
        to_unit_combo.current(1 if len(units) > 1 else 0)
variable_combo.bind('<<ComboboxSelected>>', update_unit_dropdowns)

# Conversion logic
def convert_units():
    value = value_entry.get()
    var = variable_combo.get()
    from_unit = from_unit_combo.get()
    to_unit = to_unit_combo.get()
    try:
        val = float(value)
    except Exception:
        result_label.config(text="Enter a valid numeric value.", fg="red")
        return
    if var == 'x':
        if from_unit not in length_factors or to_unit not in length_factors:
            result_label.config(text="Invalid units for length.", fg="red")
            return
        val_m = val * length_factors[from_unit]
        converted = val_m / length_factors[to_unit]
    elif var == 'D':
        if from_unit not in area_factors or to_unit not in area_factors:
            result_label.config(text="Invalid units for D.", fg="red")
            return
        val_base = val * area_factors[from_unit]
        converted = val_base / area_factors[to_unit]
    elif var == 't':
        if from_unit not in time_factors or to_unit not in time_factors:
            result_label.config(text="Invalid units for time.", fg="red")
            return
        val_base = val * time_factors[from_unit]
        converted = val_base / time_factors[to_unit]
    elif var == 'Qd':
        if from_unit not in qd_factors or to_unit not in qd_factors:
            result_label.config(text="Invalid units for Qd.", fg="red")
            return
        val_base = val * qd_factors[from_unit]
        converted = val_base / qd_factors[to_unit]
    elif var == 'T':
        converted = K_to_t(t_to_K(val, from_unit), to_unit)
    elif var == 'R':
        if from_unit not in r_factors or to_unit not in r_factors:
            result_label.config(text="Invalid units for R.", fg="red")
            return
        val_base = val * r_factors[from_unit]
        converted = val_base / r_factors[to_unit]
    else:
        result_label.config(text="Select a valid variable.", fg="red")
        return
    result_label.config(text=f"{val} {from_unit} = {converted:.6g} {to_unit}", fg="blue")

convert_btn = tk.Button(unit_frame, text="Convert", command=convert_units)
convert_btn.grid(row=4, column=0, columnspan=2)
result_label = tk.Label(unit_frame, text="")
result_label.grid(row=5, column=0, columnspan=2)

# --- Initialize GUI ---
show_fields()

# --- Start the Tkinter event loop ---
root.mainloop()
