import tkinter as tk
from tkinter import ttk

# Default values for data
data = {'avg': 1, 'std': 0.1, 'num': 1000, 'den': 0.25, 'olp': 0.0, 'dst': 'gamma', 'pbc': False, 'sar': False}


def settings_gui():
    # Function to collect values and update the data dictionary
    global data

    def apply_values():
        global data
        data = {
            "avg": float(avg_var.get()),
            "std": float(std_var.get()),
            "num": int(num_var.get()),
            "den": float(den_var.get()),
            "olp": float(olp_var.get()),
            "dst": dst_var.get().lower(),
            "pbc": pbc_var.get(),
            "sar": sar_var.get()
        }
        root.destroy()

    # Function to handle cancel
    def cancel():
        root.destroy()

    # Main window
    root = tk.Tk()
    root.title("Settings")

    # Variables for storing input
    avg_var = tk.StringVar(value=str(data['avg']))
    std_var = tk.StringVar(value=str(data['std']))
    num_var = tk.StringVar(value=str(data['num']))
    den_var = tk.DoubleVar(value=data['den'])
    olp_var = tk.DoubleVar(value=data['olp'])
    dst_var = tk.StringVar(value=data['dst'].capitalize())
    pbc_var = tk.BooleanVar(value=data['pbc'])
    sar_var = tk.BooleanVar(value=data['sar'])

    # Average Entry
    tk.Label(root, text="Avg").grid(row=0, column=0)
    tk.Entry(root, textvariable=avg_var).grid(row=0, column=1)
    tk.Label(root, text='(0.001 - inf)').grid(row=0, column=2)

    # Standard Deviation Entry
    tk.Label(root, text="CV").grid(row=1, column=0)
    tk.Entry(root, textvariable=std_var).grid(row=1, column=1)
    tk.Label(root, text='(0.001 - inf)').grid(row=1, column=2)

    # Number Entry
    tk.Label(root, text="Num").grid(row=2, column=0)
    tk.Entry(root, textvariable=num_var).grid(row=2, column=1)
    tk.Label(root, text='(1 - inf)').grid(row=2, column=2)

    # Density Entry/Slider
    tk.Label(root, text="Density").grid(row=3, column=0)
    tk.Scale(root, variable=den_var, from_=0.001, to=1, orient="horizontal", resolution=0.001).grid(row=3, column=1)
    tk.Entry(root, textvariable=den_var).grid(row=3, column=2)

    # Overlap Entry/Slider
    tk.Label(root, text="Overlap").grid(row=4, column=0)
    tk.Scale(root, variable=olp_var, from_=0.0, to=2.0, orient="horizontal", resolution=0.01).grid(row=4, column=1)
    tk.Entry(root, textvariable=olp_var).grid(row=4, column=2)

    # Distribution Dropdown
    tk.Label(root, text="Dst").grid(row=5, column=0)
    dst_options = ["Gamma", "Log-Normal", "Weibull", "Normal", "Half-Normal", "Physical 1 (Devries)",
                   "Physical 2 (Gal-Or)", "Physical 3 (Lemelich)"]
    dst_menu = ttk.Combobox(root, textvariable=dst_var, values=dst_options)
    dst_menu.grid(row=5, column=1)
    dst_menu.current(dst_options.index(data['dst'].capitalize()))

    # Periodic Boundary Condition Checkbox
    tk.Checkbutton(root, text="Periodic Boundary Conditions", variable=pbc_var).grid(row=6, column=0, columnspan=3)

    # Standardize Atomic Radii Checkbox
    tk.Checkbutton(root, text="Standardize Atomic Radii", variable=sar_var).grid(row=7, column=0, columnspan=2)

    # Buttons
    tk.Button(root, text="Apply", command=apply_values).grid(row=8, column=1)
    tk.Button(root, text="Cancel", command=cancel).grid(row=8, column=2)

    # Run the GUI
    root.mainloop()

    return data
